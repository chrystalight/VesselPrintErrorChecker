"""
Measurement engine for the Vessel Print Error Checker (no Slicer dependency: numpy, scipy, vtk).

All geometry is in one world frame (Slicer RAS, mm). The CT is given as a numpy array indexed [k, j, i]
plus its 4x4 world->IJK matrix. Centerlines are (N, 3) point arrays in the same world frame. The design
(original lumen) is a vtkPolyData in its own frame; it is registered to the print here.

Terms: a *vessel* is one connected vessel tree (one CT, one design STL); a *branch* is one centerline in it.
The print is registered to the design once, as a rigid body. Per station the engine reports the printed and design
equivalent diameters (2*sqrt(A/pi)), their difference, and the centerline error: the distance between the printed
and design centerlines under that single registration.
"""
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from scipy.ndimage import map_coordinates, gaussian_filter1d, binary_dilation, binary_closing, label
from scipy.spatial import cKDTree
from scipy.special import ndtr
from scipy.spatial.transform import Rotation

DEFAULTS = dict(
    station_spacing_mm=0.5,     # distance between cross-sections
    max_search_radius_mm=12.0,  # how far a ray may look for the wall on the first pass
    end_trim_mm=2.0,            # ignore this much at each path end
    junction_margin_mm=2.0,     # extra exclusion either side of a junction
    mismatch_mm=8.0,            # print/design centres further apart than this are treated as a failed match
    min_diameter_mm=0.0,        # design diameters below this are excluded as below CT resolution (0 = keep all)
    insert_hu=None,             # voxels >= this are an insert (e.g. a tube) and are treated as lumen (None = off)
    insert_dilate_mm=1.0,
    exclude_mm=None,            # {branch name: [[from_mm, to_mm], ...]} ranges along a centerline to leave out
)
N_RAYS = 90
ANG = np.linspace(0, 2 * np.pi, N_RAYS, endpoint=False)
DR = 0.05            # ray sampling step (mm)
MIN_CONTRAST = 250   # HU between lumen and wall for a ray to count
WALL_SEARCH = 2.5    # mm past the first rise searched for the wall peak
MAX_BAD_RAYS = 0.15  # sections with more rays finding no wall are "open"
WALL_RUN = 0.7       # mm of continuous wall needed to end a ray when an insert is present
MAX_INSERT_RAYS = 0.6  # sections where more of the wall than this is only seen through insert contact are excluded
MAX_SHIFT_MM = 3.0   # a measured centre may move at most max(this, 0.75 x lumen radius) from the drawn centerline
PASSES = 3
SOLID_HU = -450      # anything above this is solid or fluid rather than air
FLUID_AVG_MM = 1.0   # fluid-filled sections: profiles are averaged this far either side along the vessel
FLUID_MIN_CONTRAST = 60   # HU the wall must stand above the fluid in the lumen
FLUID_MIN_RADIUS = 3.0    # in blur widths: a smaller lumen never shows the fluid's own level, so its edge cannot be fitted
SMOOTH_MM = 1.5


# ----------------------------------------------------------------------------- small helpers
def nan_median_filter(x, w, wrap=False):
    """Sliding nan-median over an odd window w (wrap=True for circular data)."""
    h = w // 2
    xp = np.r_[x[-h:], x, x[:h]] if wrap else np.r_[np.full(h, np.nan), x, np.full(h, np.nan)]
    win = np.lib.stride_tricks.sliding_window_view(xp, w)
    with np.errstate(all="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            return np.nanmedian(win, axis=1)


def points_in_polygon(poly, pts):
    """Even-odd test of 2D points against a closed polygon."""
    x, y = pts[:, 0][:, None], pts[:, 1][:, None]
    x0, y0 = poly[:, 0][None], poly[:, 1][None]
    x1, y1 = np.roll(poly[:, 0], -1)[None], np.roll(poly[:, 1], -1)[None]
    cond = (y0 > y) != (y1 > y)
    with np.errstate(divide="ignore", invalid="ignore"):
        xc = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
    return (np.sum(cond & (x < xc), axis=1) % 2) == 1


def polygon(xy):
    x, y = xy[:, 0], xy[:, 1]; x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y; A = cr.sum() / 2
    if abs(A) < 1e-9: return 0.0, xy.mean(0)
    return abs(A), np.array([((x + x1) * cr).sum(), ((y + y1) * cr).sum()]) / (6 * A)


def resample(P, step):
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    q = np.arange(0, s[-1] + 1e-9, step)
    return np.c_[[np.interp(q, s, P[:, i]) for i in range(3)]].T, q


def tangents(C, step):
    T = np.gradient(gaussian_filter1d(C, SMOOTH_MM / step, axis=0, mode="nearest"), axis=0)
    return T / np.linalg.norm(T, axis=1, keepdims=True)


def frames(T):
    """Parallel-transported in-plane axes so ray k points the same way from station to station."""
    a = np.eye(3)[np.argmin(np.abs(T[0]))]
    u = [np.cross(T[0], a)]; u[0] /= np.linalg.norm(u[0])
    for t in T[1:]:
        w = u[-1] - np.dot(u[-1], t) * t
        u.append(w / np.linalg.norm(w))
    U = np.array(u); return U, np.cross(T, U)


def kabsch(A, B, w=None):
    w = np.ones(len(A)) if w is None else np.asarray(w, float)
    ca, cb = np.average(A, 0, w), np.average(B, 0, w); H = ((A - ca) * w[:, None]).T @ (B - cb); U, _, Vt = np.linalg.svd(H)
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]); R = Vt.T @ D @ U.T
    return R, cb - R @ ca


def icp(src, tgt, R=np.eye(3), t=np.zeros(3), iters=50, keep=1.0):
    tree = cKDTree(tgt)
    for _ in range(iters):
        X = src @ R.T + t; d, j = tree.query(X); m = d <= np.quantile(d, keep)
        dR, dt = kabsch(X[m], tgt[j[m]]); R, t = dR @ R, dR @ t + dt
    return R, t


def densify(P, step=0.05):
    seg = np.linalg.norm(np.diff(P, axis=0), axis=1); n = np.maximum(1, np.ceil(seg / step)).astype(int)
    return np.vstack([P[i] + np.linspace(0, 1, n[i], endpoint=False)[:, None] * (P[i + 1] - P[i]) for i in range(len(P) - 1)] + [P[-1:]])


# ----------------------------------------------------------------------------- CT
class CTVolume:
    def __init__(self, array_kji, world_to_ijk, insert_hu=None, insert_dilate_mm=1.0, spacing_mm=0.4):
        self.v = np.asarray(array_kji, dtype=np.float32).copy()
        self.M = np.asarray(world_to_ijk, float)
        self.core = self.rim = None
        if insert_hu is not None:
            # An insert (e.g. a sheath) lies inside the lumen, so it counts as lumen. `core` is the insert itself,
            # closed so that a hollow sheath's bore is included; `rim` adds its blurred edge (about one voxel).
            box = np.ones((3, 3, 3), bool)
            core = binary_closing(self.v >= insert_hu, structure=box, iterations=max(1, int(round(0.8 / spacing_mm))))
            self.core = core.astype(np.float32)
            self.rim = binary_dilation(core, structure=box, iterations=max(1, int(round(0.5 / spacing_mm)))).astype(np.float32)

    def _ijk(self, X):
        return X @ self.M[:3, :3].T + self.M[:3, 3]

    def sample(self, X):
        ijk = self._ijk(X)
        return map_coordinates(self.v, [ijk[..., 2], ijk[..., 1], ijk[..., 0]], order=1, cval=np.nan)

    def in_insert(self, X, rim=True):
        m = self.rim if rim else self.core
        if m is None: return np.zeros(X.shape[:-1], bool)
        ijk = self._ijk(X)
        return map_coordinates(m, [ijk[..., 2], ijk[..., 1], ijk[..., 0]], order=1, cval=0) > 0.5


# ----------------------------------------------------------------------------- printed lumen from CT
def lumen_section_insert(ct, c, u, v, r_max):
    """Lumen contour when an insert (e.g. a sheath) is in the lumen. Lumen = air or insert; the wall starts at the
    first stretch of at least WALL_RUN mm that is neither (the blurred edge of the insert in air is shorter than
    that and is skipped). Coming from air, the edge is the usual half-maximum between air and wall; where the insert
    touches the wall, the edge is the insert's own outer surface (half-maximum between insert and wall), so the
    search never passes through the wall behind it. Returns radii, air level, rays ending on insert contact, and
    rays that pass through the insert."""
    r = np.arange(0, r_max, DR); n_run = int(WALL_RUN / DR)
    dirs = np.cos(ANG)[:, None] * u + np.sin(ANG)[:, None] * v
    X = c + r[None, :, None] * dirs[:, None, :]
    f = ct.sample(X); rim = ct.in_insert(X); core = ct.in_insert(X, rim=False)
    L = np.nanpercentile(f[~rim], 2) if (~rim).any() else np.nanpercentile(f, 2)
    solid = f[~rim & (f > L + MIN_CONTRAST)]
    rad = np.full(N_RAYS, np.nan); contact = np.zeros(N_RAYS, bool); through = np.zeros(N_RAYS, bool)
    if solid.size < 20: return rad, L, contact, through
    Wref = np.median(solid); air_thr = (L + Wref) / 2
    lumen = core | (f < air_thr)
    for k in range(N_RAYS):
        wall = ~lumen[k] & np.isfinite(f[k])
        run = np.convolve(wall, np.ones(n_run), "valid") == n_run       # run[i]: wall from i to i + n_run
        if not run.any(): continue
        i0 = int(np.argmax(run))
        if i0 == 0: continue
        through[k] = core[k, :i0].any()
        fk = f[k]; back = max(0, i0 - int(1.0 / DR))
        if core[k, i0 - 1]:                                           # insert touching the wall
            Wk = np.median(fk[i0 + int(0.4 / DR):i0 + n_run]); Tk = np.nanmax(fk[back:i0 + 1])
            if Tk - Wk < MIN_CONTRAST: continue
            thr = (Tk + Wk) / 2; j = np.nonzero(fk[back:i0 + n_run] >= thr)[0]
            if not len(j): continue
            j = back + j[-1]
            rad[k] = r[j] + DR * (fk[j] - thr) / max(fk[j] - fk[j + 1], 1e-6) if j + 1 < len(r) else r[j]
            contact[k] = True
        else:                                                         # air to wall
            seg = fk[i0:i0 + int(WALL_SEARCH / DR)]; Wk = np.nanmax(np.where(core[k, i0:i0 + len(seg)], -np.inf, seg))
            Wk = min(Wk, Wref + 2 * MIN_CONTRAST)                     # do not let a nearby insert set the wall level
            if Wk - L < MIN_CONTRAST: continue
            thr = (L + Wk) / 2; a = np.nonzero(fk[back:i0 + n_run] < thr)[0]
            if not len(a): continue
            j = back + a[-1]                                          # last sample below the half-maximum
            if j + 1 >= len(r): continue
            rad[k] = r[j] + DR * (thr - fk[j]) / max(fk[j + 1] - fk[j], 1e-6)
    if np.isfinite(rad).sum() > N_RAYS / 2:
        med = nan_median_filter(rad, 7, wrap=True)
        rad[np.abs(rad - med) > np.maximum(0.5, 0.15 * med)] = np.nan
    # the caller treats a centre above (returned level + MIN_CONTRAST / 2) as "in the wall"; with thin slivers of air
    # beside an insert the right cut-off is the air/wall half-maximum, so return the level that gives exactly that
    return rad, air_thr - MIN_CONTRAST / 2, contact, through


def hull_radii(rad):
    """Radii of the convex outline of the valid edge points, along the ray directions. Used where an insert lies
    against the wall: the lumen is convex, but the sliver of air between insert and wall is too thin to resolve and
    would otherwise leave notches in the contour."""
    from scipy.spatial import ConvexHull
    ok = np.isfinite(rad)
    if ok.sum() < 8: return None
    pts = np.c_[rad[ok] * np.cos(ANG[ok]), rad[ok] * np.sin(ANG[ok])]
    try: H = pts[ConvexHull(pts).vertices]
    except Exception: return None
    a, b = H, np.roll(H, -1, axis=0); e = b - a
    out = np.full(N_RAYS, np.nan)
    for k in range(N_RAYS):                                    # ray from the origin against each hull edge
        d = np.array([np.cos(ANG[k]), np.sin(ANG[k])])
        den = d[0] * e[:, 1] - d[1] * e[:, 0]
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (a[:, 0] * e[:, 1] - a[:, 1] * e[:, 0]) / den
            w = (a[:, 0] * d[1] - a[:, 1] * d[0]) / den
        hit = (np.abs(den) > 1e-12) & (t > 0) & (w >= -1e-9) & (w <= 1 + 1e-9)
        if hit.any(): out[k] = t[hit].max()
    return out if np.isfinite(out).all() else None


def lumen_section(ct, c, u, v, r_max):
    if ct.rim is not None: return lumen_section_insert(ct, c, u, v, r_max)
    through = np.zeros(N_RAYS, bool)
    r = np.arange(0, r_max, DR)
    dirs = np.cos(ANG)[:, None] * u + np.sin(ANG)[:, None] * v
    f = ct.sample(c + r[None, :, None] * dirs[:, None, :])
    L = np.nanpercentile(f, 5)
    rad = np.full(N_RAYS, np.nan)
    for k in range(N_RAYS):
        fk = f[k]; up = np.nonzero(fk >= L + MIN_CONTRAST / 2)[0]
        if not len(up): continue
        i0 = up[0]; W = np.nanmax(fk[i0:i0 + int(WALL_SEARCH / DR)])
        if W - L < MIN_CONTRAST: continue
        thr = (L + W) / 2; i = np.nonzero(fk >= thr)[0][0]
        if i == 0: continue
        rad[k] = r[i - 1] + DR * (thr - fk[i - 1]) / (fk[i] - fk[i - 1])
    ins = np.zeros(N_RAYS, bool)
    if np.isfinite(rad).sum() > N_RAYS / 2:      # drop spikes against the neighbouring rays
        med = nan_median_filter(rad, 7, wrap=True)
        rad[np.abs(rad - med) > np.maximum(0.5, 0.15 * med)] = np.nan
    return rad, L, ins, through


def reseed(ct, c, u, v, r_max, px=0.1):
    """If the path point sits in the wall, move it to the nearest enclosed dark region in the plane."""
    g = np.arange(-r_max, r_max + px, px); A, B = np.meshgrid(g, g, indexing="ij")
    f = ct.sample(c + A[..., None] * u + B[..., None] * v)
    disc = A ** 2 + B ** 2 <= r_max ** 2
    L, W = np.nanpercentile(f[disc], 5), np.nanpercentile(f[disc], 90)
    if W - L < MIN_CONTRAST: return None
    lab, n = label((f < (L + W) / 2) & disc)
    rim = set(np.unique(lab[disc & ~(A ** 2 + B ** 2 <= (r_max - 2 * px) ** 2)]).tolist())
    best = None
    for k in range(1, n + 1):
        if k in rim: continue
        m = lab == k
        if m.sum() < 5: continue
        d = np.sqrt(A[m] ** 2 + B[m] ** 2).min()
        if best is None or d < best[0]: best = (d, A[m].mean(), B[m].mean())
    return None if best is None else c + best[1] * u + best[2] * v


def fill_rays(rad):
    ok = np.isfinite(rad)
    if ok.sum() < 3: return None
    a = np.r_[ANG[ok] - 2 * np.pi, ANG[ok], ANG[ok] + 2 * np.pi]
    return np.interp(ANG, a, np.tile(rad[ok], 3))


def fill_rays_circle(rad):
    ok = np.isfinite(rad)
    if ok.sum() < 8: return None
    x, y = rad[ok] * np.cos(ANG[ok]), rad[ok] * np.sin(ANG[ok])
    a, b, c = np.linalg.lstsq(np.c_[x, y, np.ones_like(x)], x ** 2 + y ** 2, rcond=None)[0]
    cc = np.array([a / 2, b / 2]); R2 = c + cc @ cc
    d = np.c_[np.cos(ANG), np.sin(ANG)]; pr = d @ cc; disc = pr ** 2 - cc @ cc + R2
    if R2 <= 0 or (disc < 0).any(): return fill_rays(rad)
    out = rad.copy(); out[~ok] = (pr + np.sqrt(disc))[~ok]; return out


# ----------------------------------------------------------------------------- fluid-filled lumen
def _halfmax_shift_table():
    """For a wall of thickness w blurred by a Gaussian of width sigma with air on both sides, the half-maximum edge
    sits slightly inside the true inner edge. Returns q = w / sigma and that offset in units of sigma (<= 0)."""
    q = np.linspace(0.3, 8, 78); lo, hi = np.full_like(q, -4.0), np.zeros_like(q); target = ndtr(q / 2) - 0.5
    for _ in range(40):
        mid = (lo + hi) / 2; below = ndtr(mid) - ndtr(mid - q) < target
        lo, hi = np.where(below, mid, lo), np.where(below, hi, mid)
    return q, (lo + hi) / 2


HALFMAX_Q, HALFMAX_SHIFT = _halfmax_shift_table()


def fluid_centre(ct, c, u, v, r_max=8.0, px=0.1):
    """Centre of the solid-or-fluid region (vessel wall plus filled lumen) that contains c, in the plane (u, v)."""
    g = np.arange(-r_max, r_max + px, px); A, B = np.meshgrid(g, g, indexing="ij")
    f = ct.sample(c + A[..., None] * u + B[..., None] * v)
    lab, _ = label(f > SOLID_HU); k = lab[len(g) // 2, len(g) // 2]
    if k == 0: return None
    m = lab == k
    if m[0].any() or m[-1].any() or m[:, 0].any() or m[:, -1].any(): return None     # joined to something bigger
    return c + A[m].mean() * u + B[m].mean() * v


def lumen_section_fluid(ct, c, t, u, v, r_max=8.0):
    """Lumen contour where the lumen holds fluid (e.g. water left in after flushing) instead of air.
    Fluid and wall differ by only about 100 HU, so a half-maximum on single rays is not usable. Instead:
      1. profiles are averaged over +-FLUID_AVG_MM along the vessel;
      2. the outer wall edge (wall to outside air, a large step) is found on every ray;
      3. all rays, aligned at that edge, give one low-noise profile, fitted with a blurred wall model
         (fluid level, wall level, wall thickness, blur width);
      4. with those levels fixed, each ray is fitted for its inner and outer edge; the wall thickness is smoothed
         around the circumference and subtracted from the outer edge;
      5. the edge is moved to where the half-maximum would lie if the lumen held air, so fluid-filled and
         air-filled stations are measured by the same convention.
    Returns (radii, info), radii being NaN where a ray has no clean outer edge, or None if this is not a
    fluid-filled lumen."""
    import warnings
    from scipy.optimize import least_squares
    r = np.arange(0, r_max, DR); dirs = np.cos(ANG)[:, None] * u + np.sin(ANG)[:, None] * v
    f = np.mean([ct.sample(c + dz * t + r[None, :, None] * dirs[:, None, :])
                 for dz in np.arange(-FLUID_AVG_MM, FLUID_AVG_MM + 1e-6, 0.25)], 0)
    if not f[0, 0] > SOLID_HU: return None
    ro = np.full(N_RAYS, np.nan); Lout = np.full(N_RAYS, np.nan); n1, n2 = int(1.0 / DR), int(2.0 / DR)
    for k in range(N_RAYS):
        fk = f[k]; below = np.nonzero(~(fk >= SOLID_HU))[0]
        if not len(below) or below[0] == 0: continue
        j = below[0]; tail = fk[j + n1:j + n2]
        if not np.isfinite(fk[j]) or len(tail) < 5 or not np.isfinite(tail).all() or tail.max() > -600: continue
        ro[k] = r[j - 1] + DR * (fk[j - 1] - SOLID_HU) / (fk[j - 1] - fk[j]); Lout[k] = np.median(tail)
    med = nan_median_filter(ro, 7, wrap=True)
    ro[np.abs(ro - med) > np.maximum(0.5, 0.15 * med)] = np.nan
    ok = np.isfinite(ro)
    if ok.mean() < 1 - MAX_BAD_RAYS: return None
    x = np.arange(-4.0, 1.5, DR)                                       # distance from the outer edge
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        g = np.nanmedian([np.interp(ro[k] + x, r, f[k], left=np.nan) for k in np.nonzero(ok)[0]], 0)
    Lo = float(np.median(Lout[ok])); m = np.isfinite(g) & (x > -np.median(ro[ok]) + 0.3)
    if m.sum() < 30: return None
    Lw0 = float(np.median(g[m][:6]))

    def model(q):
        Lw, W, w, sg, x0 = q; z = x[m] - x0
        return Lw + (W - Lw) * ndtr((z + w) / sg) - (W - Lo) * ndtr(z / sg) - g[m]
    lo_b, hi_b = np.array([-1100, -300, 0.2, 0.2, -0.5]), np.array([800, 2500, 3.5, 0.9, 0.5])
    q0 = np.clip([Lw0, max(Lw0 + 130, 200), 1.1, 0.4, 0.0], lo_b + 1e-3, hi_b - 1e-3)
    Lw, W, w, sg, x0 = least_squares(model, q0, bounds=(lo_b, hi_b)).x
    # a fluid-filled lumen: inside is much nearer the wall level than the outside air, yet darker than the wall
    if W - Lw < FLUID_MIN_CONTRAST or Lw < Lo + 0.8 * (W - Lo) or not 0.3 <= w <= 3.0: return None
    A = W - Lw; Bk = np.where(ok, W - Lout, 0.0); rq = np.where(ok, ro + x0, 1.0); ri = np.clip(rq - w, 0.2, None)
    mask = ok[:, None] & (r[None] < (rq + 1.2)[:, None]) & np.isfinite(f); f0 = np.where(mask, f, 0.0)
    pdf = lambda z: np.exp(-0.5 * z * z) / np.sqrt(2 * np.pi)
    for _ in range(15):                                                # damped Gauss-Newton on (inner, outer) per ray
        zi, zo = (r[None] - ri[:, None]) / sg, (r[None] - rq[:, None]) / sg
        res = np.where(mask, Lw + A * ndtr(zi) - Bk[:, None] * ndtr(zo) - f0, 0.0)
        Ji = np.where(mask, -A * pdf(zi) / sg, 0.0); Jo = np.where(mask, Bk[:, None] * pdf(zo) / sg, 0.0)
        a, b, d = (Ji * Ji).sum(1), (Ji * Jo).sum(1), (Jo * Jo).sum(1); gi, go = (Ji * res).sum(1), (Jo * res).sum(1)
        a = a * 1.05 + 1e-6; d = d * 1.05 + 1e-6; det = a * d - b * b
        ri = np.clip(ri + np.clip(-(d * gi - b * go) / det, -0.2, 0.2), 0.1, r_max)
        rq = np.clip(rq + np.clip(-(a * go - b * gi) / det, -0.2, 0.2), 0.3, r_max)
    wk = nan_median_filter(np.where(ok, rq - ri, np.nan), 15, wrap=True)
    rad = rq - wk + sg * np.interp(wk / sg, HALFMAX_Q, HALFMAX_SHIFT)
    rad[~ok | ~(rad > 0.2)] = np.nan
    if not np.nanmedian(rad) >= FLUID_MIN_RADIUS * sg: return None
    inner = (r[None] < (rad - 0.4)[:, None]) & np.isfinite(f)          # air inside: the lumen is only partly filled
    if inner.any() and np.mean(f[inner] < Lw - 200) > 0.05: return None
    return rad, dict(level=float(Lw), contrast=float(W - Lw), sigma=float(sg), wall=float(w))


def measure_fluid_station(ct, c0, t, u, v, lim):
    """Fluid-filled station: centre on the wall-plus-lumen region, measure, re-centre on the lumen, measure again.
    Returns (centre used for the rays, radii, info) or None."""
    c = fluid_centre(ct, c0, u, v)
    if c is None or np.linalg.norm(c - c0) > lim: return None
    for _ in range(2):
        m = lumen_section_fluid(ct, c, t, u, v)
        if m is None: return None
        rad, info = m; rf = fill_rays(rad)
        if rf is None: return None
        _, cen = polygon(np.c_[rf * np.cos(ANG), rf * np.sin(ANG)]); c2 = c + cen[0] * u + cen[1] * v
        if np.linalg.norm(c2 - c0) > lim: return None
        if np.linalg.norm(c2 - c) < 0.05: break
        c = c2
    return c, rad, info


def measure_print_path(ct, P0, step, r_search, max_diameter=None):
    """Iteratively re-centred cross-sections along one centerline.
    Ray length: max_diameter/2 + 1 mm if a caliper diameter is given; otherwise r_search on the first pass and
    the local lumen radius + 2.5 mm afterwards (about a wall plus margin; stops rays running into a neighbouring vessel).
    Re-centring may move a point at most max(3 mm, 0.75 x lumen radius) from the drawn centerline. Where no lumen is
    found that close (e.g. the lumen is filled, or the drawn path runs through the wall), the station is marked `lost`
    rather than letting the centre jump into a neighbouring vessel."""
    C, s = resample(P0, step); n = len(C); C0 = C.copy()
    rmax = np.full(n, max_diameter / 2 + 1.0 if max_diameter else r_search)
    lim = np.full(n, MAX_SHIFT_MM)
    for it in range(PASSES):
        T = tangents(C, step); U, V = frames(T)
        newC, D = C.copy(), np.full(n, np.nan); bad, insert = np.ones(n), np.zeros(n)
        rings = np.full((n, N_RAYS, 3), np.nan); lost = np.zeros(n, bool); valid = np.zeros((n, N_RAYS), bool)
        for i in range(n):
            rad, L, ins, thru = lumen_section(ct, C[i], U[i], V[i], rmax[i])
            if ct.sample(C[i][None])[0] > L + MIN_CONTRAST / 2 and not ct.in_insert(C[i][None])[0]:
                c2 = reseed(ct, C[i], U[i], V[i], min(rmax[i], lim[i] + 1.0))
                if c2 is None or np.linalg.norm(c2 - C0[i]) > lim[i]:
                    lost[i] = True; newC[i] = C0[i]; continue
                C[i] = c2; rad, _, ins, thru = lumen_section(ct, C[i], U[i], V[i], rmax[i])
            bad[i] = np.mean(~np.isfinite(rad)); insert[i] = np.mean(ins)
            valid[i] = np.isfinite(rad)
            rf = (hull_radii(rad) if thru.mean() > 0.05 else None)
            if rf is None: rf = fill_rays(rad)
            if rf is None: continue
            xy = np.c_[rf * np.cos(ANG), rf * np.sin(ANG)]; A, cen = polygon(xy)
            ci = C[i] + cen[0] * U[i] + cen[1] * V[i]
            if np.linalg.norm(ci - C0[i]) > lim[i]:                 # drifted towards another lumen
                lost[i] = True; newC[i] = C0[i]; bad[i] = 1.0; valid[i] = False; continue
            D[i] = 2 * np.sqrt(A / np.pi); newC[i] = ci
            rings[i] = C[i] + xy[:, :1] * U[i] + xy[:, 1:] * V[i]
        C = newC
        if not max_diameter and np.isfinite(D).any():
            Rl = nan_median_filter(D / 2, 2 * int(2.5 / step) + 1)
            Rl[~np.isfinite(Rl)] = np.nanmedian(D / 2)
            rmax = np.minimum(r_search, Rl + 2.5)
            lim = np.maximum(MAX_SHIFT_MM, 0.75 * Rl)
    # Stations with no air lumen near the path may hold fluid instead: measure those from the fluid/wall edge.
    fluid = np.zeros(n, bool); fluid_info = [None] * n
    if ct.rim is None and lost.any():
        T = tangents(C, step); U, V = frames(T)
        found = {i: measure_fluid_station(ct, C0[i], T[i], U[i], V[i], lim[i]) for i in np.nonzero(lost)[0]}
        for i, m in found.items():
            if m is None: continue
            c, rad, q = m; rf = fill_rays(rad); fluid_info[i] = q
            xy = np.c_[rf * np.cos(ANG), rf * np.sin(ANG)]; A, cen = polygon(xy)
            C[i] = c + cen[0] * U[i] + cen[1] * V[i]; D[i] = 2 * np.sqrt(A / np.pi)
            rings[i] = c + xy[:, :1] * U[i] + xy[:, 1:] * V[i]; valid[i] = np.isfinite(rad)
            bad[i] = np.mean(~np.isfinite(rad)); lost[i] = False; fluid[i] = True
    return dict(s=s, C=C, T=tangents(C, step), D=D, bad=bad, insert=insert, rings=rings, lost=lost, valid=valid, fluid=fluid, fluid_info=fluid_info)


# ----------------------------------------------------------------------------- design (original) sections
class DesignSlicer:
    def __init__(self, poly):
        self.plane = vtk.vtkPlane()
        self.cut = vtk.vtkCutter(); self.cut.SetCutFunction(self.plane); self.cut.SetInputData(poly)
        self.strip = vtk.vtkStripper(); self.strip.SetInputConnection(self.cut.GetOutputPort()); self.strip.JoinContiguousSegmentsOn()

    def section(self, c, n, d_max):
        """Design lumen loop in the plane (c, n): the loop enclosing c, else the nearest one."""
        self.plane.SetOrigin(*c); self.plane.SetNormal(*n); self.strip.Update()
        out = self.strip.GetOutput()
        if out.GetNumberOfPoints() < 3: return None
        P = vtk_to_numpy(out.GetPoints().GetData()); lines = out.GetLines(); lines.InitTraversal()
        u = np.cross(n, np.eye(3)[np.argmin(np.abs(n))]); u /= np.linalg.norm(u); v = np.cross(n, u)
        best, ids = None, vtk.vtkIdList()
        while lines.GetNextCell(ids):
            idx = [ids.GetId(k) for k in range(ids.GetNumberOfIds())]
            if len(idx) < 6: continue
            pts = P[idx]; xy = np.c_[(pts - c) @ u, (pts - c) @ v]
            A, cen = polygon(xy)
            inside = bool(points_in_polygon(xy, np.zeros((1, 2)))[0]); dist = np.linalg.norm(cen)
            key = (0 if inside else 1, dist)
            if dist < d_max and (best is None or key < best[0]): best = (key, A, c + cen[0] * u + cen[1] * v, pts)
        return None if best is None else (2 * np.sqrt(best[1] / np.pi), best[2], best[3])


def crop(poly, C, margin):
    lo, hi = np.nanmin(C, 0) - margin, np.nanmax(C, 0) + margin
    box = vtk.vtkBox(); box.SetBounds(lo[0], hi[0], lo[1], hi[1], lo[2], hi[2])
    ex = vtk.vtkExtractPolyDataGeometry(); ex.SetInputData(poly); ex.SetImplicitFunction(box)
    ex.ExtractInsideOn(); ex.ExtractBoundaryCellsOn(); ex.Update()
    out = vtk.vtkPolyData(); out.DeepCopy(ex.GetOutput()); return out


def measure_design_path(poly, C0, T0, d_max, step, passes=PASSES):
    sl = DesignSlicer(poly); C, T = C0.copy(), T0.copy(); n = len(C)
    D = np.full(n, np.nan); loops = [None] * n
    for _ in range(passes):
        newC = C.copy()
        for i in range(n):
            r = sl.section(C[i], T[i], d_max)
            if r is None: D[i] = np.nan; loops[i] = None; continue
            D[i], newC[i], loops[i] = r
        C = newC; T = tangents(C, step)
    return C, D, loops


# ----------------------------------------------------------------------------- registration
def register(wall_pts, design_pts, progress=None):
    """Print wall points -> design surface. 24 axis-aligned orientations x 9 placements, annealed ICP
    (only pairs closer than a shrinking threshold), scored by the fraction of points within 1 mm."""
    tree = cKDTree(design_pts); Y = design_pts
    X = wall_pts[np.random.default_rng(0).choice(len(wall_pts), min(10000, len(wall_pts)), replace=False)]; Xc = X[:1500]
    anchors = [Y.mean(0)]
    for _ in range(8): anchors.append(Y[np.argmax(np.min([np.linalg.norm(Y - a, axis=1) for a in anchors], 0))])

    def anneal(P, R, t, sched):
        for thr in sched:
            Z = P @ R.T + t; d, j = tree.query(Z, distance_upper_bound=thr); m = np.isfinite(d)
            if m.sum() < 30: break
            dR, dt = kabsch(Z[m], Y[j[m]]); R, t = dR @ R, dR @ t + dt
        return R, t

    inl = lambda R, t, P: float(np.mean(tree.query(P @ R.T + t, distance_upper_bound=1.0)[0] < 1.0))
    starts, rots = [], Rotation.create_group("O").as_matrix()
    for k, R0 in enumerate(rots):
        if progress: progress(k / len(rots))
        for a in anchors:
            R, t = anneal(Xc, R0, a - R0 @ Xc.mean(0), list(np.geomspace(20, 2, 20))); starts.append((inl(R, t, Xc), R, t))
    best = None
    for _, R, t in sorted(starts, key=lambda z: -z[0])[:3]:
        R, t = anneal(X, R, t, list(np.geomspace(3, 1.5, 10)) + [1.5] * 30); f = inl(R, t, X)
        if best is None or f > best[2]: best = (R, t, f)
    R, t, f = best
    return R, t, dict(rotation_deg=float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))),
                      median_wall_distance_mm=float(np.median(tree.query(X @ R.T + t)[0])), wall_within_1mm=f)


def wall_offset(P, valid, u, v, tree, Y, N, shift0=None):
    """Fit  signed distance to design = dR + shift . normal  over the wall points P of one station (design frame).
    `valid` marks which of the N_RAYS rays the points belong to. Returns (dR, |shift|, rms) or None if the wall
    seen is too little or too one-sided to separate a radius error from a sideways shift."""
    if len(P) < 0.4 * N_RAYS: return None
    shift = np.zeros(3) if shift0 is None else np.asarray(shift0, float).copy()
    for _ in range(4):                                    # re-pair after removing the sideways shift found so far
        Q = P - shift; d, j = tree.query(Q); keep = d < 2.5
        if keep.sum() < 0.4 * N_RAYS: return None
        n = N[j[keep]]; sd = ((Q[keep] - Y[j[keep]]) * n).sum(1)
        A = np.c_[np.ones(keep.sum()), n @ u, n @ v]
        x, *_ = np.linalg.lstsq(A, sd, rcond=None)
        shift = shift + x[1] * u + x[2] * v
    seen = np.zeros(N_RAYS, bool); seen[np.nonzero(valid)[0][keep]] = True
    gap = 0; run = 0
    for b in np.r_[seen, seen]:                           # longest stretch of rays without wall (circular)
        run = 0 if b else run + 1; gap = max(gap, run)
    if gap > N_RAYS * 200 / 360: return None
    rms = float(np.sqrt(np.mean((A @ x - sd) ** 2)))
    if rms > 0.5: return None
    return float(x[0]), float(np.linalg.norm(shift)), rms


# ----------------------------------------------------------------------------- whole module
def run(ct, centerlines, design_poly, params=None, progress=None):
    """
    ct: CTVolume. centerlines: {name: (N,3) points, or (points, max_diameter_mm or None)}. design_poly: vtkPolyData.
    Returns dict(stations=list of dicts, paths={name: inspect data}, R, t (print->design), registration=dict).
    """
    p = dict(DEFAULTS); p.update(params or {}); step = p["station_spacing_mm"]
    say = progress or (lambda frac, msg: None)
    names = list(centerlines)

    pr = {}
    for k, name in enumerate(names):
        say(0.6 * k / len(names), f"Measuring printed lumen: {name}")
        pts, dmax = centerlines[name] if isinstance(centerlines[name], tuple) else (centerlines[name], None)
        pr[name] = measure_print_path(ct, np.asarray(pts, float), step, p["max_search_radius_mm"], dmax)
        pr[name]["dmax"] = dmax or np.nanpercentile(pr[name]["D"], 95) + 3.0

    for name, r in pr.items():                     # exclusions known from the print alone
        flag = np.array([""] * len(r["s"]), dtype=object)
        flag[r["bad"] > MAX_BAD_RAYS] = "open"; flag[~np.isfinite(r["D"])] = "open"
        flag[r["insert"] > MAX_INSERT_RAYS] = "insert contact"
        flag[r["lost"]] = "no lumen near path"
        # sections next to a stretch with no measurable lumen are often partly filled (fluid surface) -> also excluded
        nolum = np.isin(flag, ["open", "no lumen near path"]); k = int(p["junction_margin_mm"] / step)
        edge = (np.convolve(nolum, np.ones(2 * k + 1), "same") > 0) & (flag == "")
        flag[edge] = "next to unmeasurable section"
        fl = np.convolve(r["fluid"], np.ones(2 * k + 1), "same"); n_meas = np.convolve(np.isfinite(r["D"]) & ~nolum, np.ones(2 * k + 1), "same")
        flag[(fl > 0) & (fl < n_meas) & (flag == "")] = "next to fluid surface"
        flag[(r["s"] < p["end_trim_mm"]) | (r["s"] > r["s"][-1] - p["end_trim_mm"])] = "end"
        for lo, hi in (p.get("exclude_mm") or {}).get(name, []):     # user-chosen ranges, e.g. a water-filled stretch
            flag[(r["s"] >= lo) & (r["s"] <= hi)] = "user excluded"
        # Another branch's centerline passing through this cross-section means either
        #  - a shared segment: both centerlines run along the same lumen for a while (e.g. ICA and ECA paths both
        #    starting in the common carotid). It is measured once, in the branch listed first, and marked
        #    "shared with X" in the others. Only stretches of at least 5 mm and 1.5 diameters count as shared;
        #  - a junction: the other branch leaves at an angle, or only briefly runs alongside where it starts inside
        #    its parent. The section straddles both lumens, so it is excluded (plus a margin either side).
        n_st = len(r["s"]); junc = np.zeros(n_st, bool); par = np.full(n_st, "", dtype=object)
        for other, o in pr.items():
            if other == name: continue
            for i in range(n_st):
                if not np.isfinite(r["rings"][i, 0, 0]): continue
                rel = o["C"] - r["C"][i]; near = np.nonzero(np.abs(rel @ r["T"][i]) < 1.0)[0]
                if not len(near): continue
                U, V = frames(r["T"][i:i + 1]); ring = r["rings"][i] - r["C"][i]
                inside = near[points_in_polygon(np.c_[ring @ U[0], ring @ V[0]], np.c_[rel[near] @ U[0], rel[near] @ V[0]])]
                if not len(inside): continue
                if np.max(np.abs(o["T"][inside] @ r["T"][i])) > np.cos(np.radians(30)):   # running alongside
                    if not par[i]: par[i] = other
                else:
                    junc[i] = True
        shared = np.full(n_st, "", dtype=object); i = 0
        while i < n_st:                                   # runs of "running alongside the same other branch"
            if not par[i]: i += 1; continue
            j = i
            while j + 1 < n_st and par[j + 1] == par[i]: j += 1
            run_mm = (j - i + 1) * step; d_run = np.nanmedian(r["D"][i:j + 1]) if np.isfinite(r["D"][i:j + 1]).any() else 0
            if run_mm < max(5.0, 1.5 * d_run): junc[i:j + 1] = True                  # brief: a branch origin
            elif names.index(par[i]) < names.index(name): shared[i:j + 1] = f"shared with {par[i]}"
            i = j + 1
        if junc.any():
            junc = np.convolve(junc, np.ones(2 * int(p["junction_margin_mm"] / step) + 1), "same") > 0
            flag[junc & (flag == "")] = "junction"
        sh = (shared != "") & (flag == ""); flag[sh] = shared[sh]
        r["flag"] = flag; r["junc"] = junc

    # Junction stations that can still be measured: the branch passes a side branch's opening (or has just left its
    # parent). Stations at the very start or end of a branch that still lie inside another branch's lumen are not
    # used: that wall belongs to the other branch and is measured there.
    for name, r in pr.items():
        n_st = len(r["s"]); inside = np.zeros(n_st, bool)
        for other, o in pr.items():
            if other == name: continue
            good = np.isfinite(o["D"]) & ~o["junc"] & np.isin(o["flag"], ["", "junction"])
            if good.sum() < 3: continue
            Ro = np.interp(np.arange(len(o["s"])), np.nonzero(good)[0], o["D"][good] / 2)
            d, j = cKDTree(o["C"]).query(r["C"])
            from_end = np.minimum(r["s"], r["s"][-1] - r["s"])
            inside |= (d < Ro[j]) & (from_end <= Ro[j] + 2.0)
        k = int(p["junction_margin_mm"] / step)        # not next to a filled or partly filled stretch
        wet = np.convolve(r["lost"] | r["fluid"], np.ones(2 * k + 1), "same") > 0
        r["jcand"] = (r["junc"] & ~inside & ~wet & (r["valid"].mean(1) >= 0.4)
                      & np.isin(r["flag"], ["junction", "open", "next to unmeasurable section"]))

    say(0.6, "Registering print to design")
    wall = np.vstack([r["rings"][np.isin(r["flag"], ["", "junction", "insert contact"])].reshape(-1, 3) for r in pr.values()])
    wall = wall[np.isfinite(wall).all(1)]
    design_pts = vtk_to_numpy(design_poly.GetPoints().GetData()).astype(float)
    R, t, reg = register(wall, design_pts, lambda f: say(0.6 + 0.2 * f, "Registering print to design"))

    nf = vtk.vtkPolyDataNormals(); nf.SetInputData(design_poly); nf.SplittingOff(); nf.ConsistencyOn()
    nf.AutoOrientNormalsOn(); nf.ComputePointNormalsOn(); nf.Update()
    d_pts = vtk_to_numpy(nf.GetOutput().GetPoints().GetData()).astype(float)
    d_nrm = vtk_to_numpy(nf.GetOutput().GetPointData().GetNormals()).astype(float); d_tree = cKDTree(d_pts)

    stations, paths = [], {}
    for k, name in enumerate(names):
        say(0.8 + 0.2 * k / len(names), f"Sectioning design: {name}")
        r = pr[name]; use = r["flag"] == ""
        Cp, Tp = r["C"] @ R.T + t, r["T"] @ R.T               # the one rigid-body registration of the whole vessel
        sub = crop(design_poly, Cp, 2 * r["dmax"])
        Co, Do, loops = measure_design_path(sub, Cp, Tp, r["dmax"], step)
        mismatch = ~(np.linalg.norm(Co - Cp, axis=1) <= p["mismatch_mm"])
        ok = use & np.isfinite(Do) & ~mismatch
        dev = np.full(len(Cp), np.nan)
        if ok.sum() > 5: dev[ok] = cKDTree(densify(Co[ok])).query(Cp[ok])[0]
        small = Do < (p["min_diameter_mm"] or 0)
        trend = nan_median_filter(Do, 2 * int(3 / step) + 1)
        artefact = np.abs(Do - trend) > np.maximum(0.3, 0.1 * trend)
        ok &= ~small & ~artefact
        # Wall offset: signed distance from the printed wall points to the design surface, fitted per station as
        # offset = dR + shift . normal. dR is the radius error (2 dR is a diameter error) and |shift| the distance
        # between the printed and design centres there. It needs only part of the wall, so it also works where the
        # cross-section is open towards a side branch.
        # The two estimates are independent, so ordinary stations where they disagree by more than 0.3 mm are
        # treated as matched to the wrong place in the design and left out.
        Uc, Vc = frames(Tp); wall = [None] * len(Cp)
        off = Cp - Co; idx = np.arange(len(Cp))                       # starting guess for the sideways shift
        if ok.sum() >= 2: off = np.c_[[np.interp(idx, idx[ok], (Cp - Co)[ok, c]) for c in range(3)]].T
        for i in np.nonzero(ok | r["jcand"])[0]:
            wall[i] = wall_offset(r["rings"][i][r["valid"][i]] @ R.T + t, r["valid"][i], Uc[i], Vc[i], d_tree, d_pts, d_nrm,
                                  off[i] if np.isfinite(off[i]).all() else None)
        agree = np.array([w is not None and abs(2 * w[0] - (r["D"][i] - Do[i])) <= 0.3 for i, w in enumerate(wall)])
        mismatch |= ok & ~agree; ok &= agree
        Dprint, Dref, at_j = r["D"].copy(), Do.copy(), np.zeros(len(Cp), bool)
        if ok.sum() >= 2:
            Dint = np.interp(np.arange(len(Cp)), np.nonzero(ok)[0], Do[ok])
            dev_int = np.interp(np.arange(len(Cp)), np.nonzero(ok)[0], dev[ok])
            for i in np.nonzero(r["jcand"] & ~ok)[0]:
                if wall[i] is None or wall[i][1] > p["mismatch_mm"] or wall[i][2] > 0.25: continue   # poor fit: skip
                if abs(wall[i][1] - dev_int[i]) > 0.5: continue    # centre jumps away from its neighbours: wrong match
                Dref[i] = Dint[i]; Dprint[i] = Dint[i] + 2 * wall[i][0]; dev[i] = wall[i][1]
                at_j[i] = True; ok[i] = True
        for i in range(len(Cp)):
            why = "" if at_j[i] else (r["flag"][i] or ("no design section" if not np.isfinite(Do[i]) else "design mismatch" if mismatch[i]
                                      else "below CT resolution" if small[i] else "design artefact" if artefact[i] else ""))
            stations.append(dict(branch=name, station=i, s_mm=float(r["s"][i]), excluded=why, at_junction=int(at_j[i]), edited=0,
                                 D_print_mm=float(Dprint[i]), D_design_mm=float(Dref[i]), dD_mm=float(Dprint[i] - Dref[i]),
                                 dD_pct=float(100 * (Dprint[i] - Dref[i]) / Dref[i]) if Dref[i] > 0 else np.nan,
                                 cl_dev_mm=float(dev[i]),
                                 insert_hidden_frac=float(r["insert"][i]), fluid_filled=int(r["fluid"][i]),
                                 x=float(r["C"][i, 0]), y=float(r["C"][i, 1]), z=float(r["C"][i, 2])))
        paths[name] = dict(s=r["s"], C=r["C"], T=r["T"], rings=r["rings"], ok=ok, valid=r["valid"], at_junction=at_j,
                           D=r["D"], wall=wall, Co=Co, Do=Dref, fluid=r["fluid"], fluid_info=r["fluid_info"],
                           design_loops_in_ct=[None if L is None else (L - t) @ R for L in loops])
    say(1.0, "Done")
    return dict(stations=stations, paths=paths, R=R, t=t, registration=reg)


def lumen_surface(ct, res):
    """Printed lumen as one continuous surface (junctions included), segmented in 3D from the CT.
    Air (or insert) voxels below the lumen/wall half-maximum are kept if they lie within the local lumen radius
    + 1 mm of a measured centerline and are connected to it; the surface is the half-maximum iso-surface, smoothed.
    Each vertex carries the diameter error and centerline error of the nearest included station (NaN if none is
    near), so the surface can be coloured by either. Returns vtkPolyData in the CT's world frame."""
    from vtk.util.numpy_support import numpy_to_vtk
    C, Rl, Tout, isend, okc, vd, vc = [], [], [], [], [], [], []
    for name, p in res["paths"].items():
        rows = [st for st in res["stations"] if st["branch"] == name]; n = len(p["s"])
        g = p["ok"] & np.isfinite(p["D"])           # radius from included stations only: elsewhere (junctions, open
        if g.sum() >= 2:                            # sections) the measured size can be inflated and let in outside air
            r = np.interp(np.arange(n), np.nonzero(g)[0], nan_median_filter(np.where(g, p["D"] / 2, np.nan), 11)[g])
        else:
            r = np.full(n, np.nanmedian(p["D"] / 2) if np.isfinite(p["D"]).any() else 2.0)
        e = np.zeros(n, bool); e[[0, -1]] = True; to = p["T"].copy(); to[0] = -to[0]
        C.append(p["C"]); Rl.append(r); Tout.append(to); isend.append(e)
        okc.append(np.array([st["excluded"] == "" for st in rows]))
        vd.append(np.array([st["dD_mm"] for st in rows])); vc.append(np.array([st["cl_dev_mm"] for st in rows]))
    C, Rl, Tout, isend, okc, vd, vc = (np.concatenate(a) for a in (C, Rl, Tout, isend, okc, vd, vc))
    Minv = np.linalg.inv(ct.M); sp = np.linalg.norm(Minv[:3, :3], axis=0).min()
    ijk = C @ ct.M[:3, :3].T + ct.M[:3, 3]; pad = int(np.ceil((Rl.max() + 3.0) / sp))
    shape = np.array(ct.v.shape)[::-1]                                  # (ni, nj, nk)
    lo = np.clip(np.floor(ijk.min(0)).astype(int) - pad, 0, shape - 1); hi = np.clip(np.ceil(ijk.max(0)).astype(int) + pad + 1, 1, shape)
    sub = ct.v[lo[2]:hi[2], lo[1]:hi[1], lo[0]:hi[0]].astype(np.float32).copy()
    L = float(np.percentile(sub, 2))
    if ct.rim is not None: sub[ct.rim[lo[2]:hi[2], lo[1]:hi[1], lo[0]:hi[0]] > 0.5] = L
    rings = np.vstack([p["rings"][p["ok"]] for p in res["paths"].values()]); cen = np.vstack([np.repeat(p["C"][p["ok"]], N_RAYS, 0) for p in res["paths"].values()])
    out = rings.reshape(-1, 3) - cen; nrm = np.linalg.norm(out, axis=1, keepdims=True)
    good = np.isfinite(out).all(1) & (nrm[:, 0] > 0)
    W = float(np.nanmedian(ct.sample((rings.reshape(-1, 3) + 0.7 * out / nrm)[good]))) if good.sum() > 50 else L + 800
    thr = (L + W) / 2
    # Fluid-filled stations have no air to segment: draw their measured contours into the volume instead.
    for p in res["paths"].values():
        for i in np.nonzero(p["fluid"] & np.isfinite(p["rings"][:, 0, 0]))[0]:
            c, t = p["C"][i], p["T"][i]; q = p["rings"][i] - c; rr = np.linalg.norm(q, axis=1)
            u = q[0] - (q[0] @ t) * t; u /= np.linalg.norm(u); v = np.cross(t, u)
            a0 = (np.round(c @ ct.M[:3, :3].T + ct.M[:3, 3]).astype(int) - lo); h = int(np.ceil((rr.max() + 1.0) / sp))
            b0, b1 = np.maximum(a0 - h, 0), np.minimum(a0 + h + 1, np.array(sub.shape)[::-1])
            if (b1 <= b0).any(): continue
            I, J, K = np.meshgrid(*[np.arange(b0[d], b1[d]) for d in range(3)], indexing="ij")
            X = (np.stack([I, J, K], -1) + lo) @ Minv[:3, :3].T + Minv[:3, 3] - c
            ang = np.arctan2(X @ v, X @ u); ring_ang = np.arctan2(q @ v, q @ u); o = np.argsort(ring_ang)
            edge = np.interp(ang, np.r_[ring_ang[o] - 2 * np.pi, ring_ang[o], ring_ang[o] + 2 * np.pi], np.tile(rr[o], 3))
            val = np.clip(thr + (np.hypot(X @ u, X @ v) - edge) * 800.0, L, None)
            slab = np.abs(X @ t) <= 0.3
            blk = sub[b0[2]:b1[2], b0[1]:b1[1], b0[0]:b1[0]].transpose(2, 1, 0)    # a view of sub, indexed [i, j, k]
            blk[slab] = np.minimum(blk[slab], val[slab])
    air = sub < thr; kji = np.argwhere(air)
    X = (kji[:, ::-1] + lo) @ Minv[:3, :3].T + Minv[:3, 3]
    d, j = cKDTree(C).query(X); keep = d < Rl[j] + 1.0
    beyond = isend[j] & (((X - C[j]) * Tout[j]).sum(1) > 0.5)           # past the end of a centerline
    keep &= ~beyond
    mask = np.zeros(sub.shape, bool); mask[tuple(kji[keep].T)] = True
    lab, _ = label(mask)
    seeds = np.round(ijk).astype(int) - lo; inb = ((seeds >= 0) & (seeds < np.array(sub.shape)[::-1])).all(1)
    ids = np.unique(lab[seeds[inb, 2], seeds[inb, 1], seeds[inb, 0]]); ids = ids[ids > 0]
    comp = np.isin(lab, ids); shell = binary_dilation(comp, iterations=2)
    vol = sub.copy(); high = thr + 400.0
    vol[~shell] = high; vol[shell & air & ~comp] = high
    img = vtk.vtkImageData(); img.SetDimensions(*vol.shape[::-1]); img.SetOrigin(*lo.astype(float)); img.SetSpacing(1, 1, 1)
    arr = numpy_to_vtk(vol.ravel(), deep=True); arr.SetName("HU"); img.GetPointData().SetScalars(arr)
    fe = vtk.vtkFlyingEdges3D(); fe.SetInputData(img); fe.SetValue(0, thr); fe.ComputeNormalsOff()
    tr = vtk.vtkTransform(); tr.SetMatrix(Minv.ravel().tolist())
    tf = vtk.vtkTransformPolyDataFilter(); tf.SetTransform(tr); tf.SetInputConnection(fe.GetOutputPort())
    sm = vtk.vtkWindowedSincPolyDataFilter(); sm.SetInputConnection(tf.GetOutputPort()); sm.SetNumberOfIterations(20)
    sm.SetPassBand(0.1); sm.BoundarySmoothingOff(); sm.NonManifoldSmoothingOn(); sm.NormalizeCoordinatesOn()
    nf = vtk.vtkPolyDataNormals(); nf.SetInputConnection(sm.GetOutputPort()); nf.SplittingOff(); nf.AutoOrientNormalsOn(); nf.Update()
    poly = vtk.vtkPolyData(); poly.DeepCopy(nf.GetOutput())
    if poly.GetNumberOfPoints() and okc.any():
        V = vtk_to_numpy(poly.GetPoints().GetData()).astype(float)
        d, j = cKDTree(C[okc]).query(V); near = d < 1.5 * Rl[okc][j] + 1.5
        for key, vals in (("dD_mm", vd), ("cl_dev_mm", vc)):
            a = numpy_to_vtk(np.where(near, vals[okc][j], np.nan), deep=True); a.SetName(key); poly.GetPointData().AddArray(a)
    return poly


def bland_altman_rows(stations):
    """One row per included station, ready for a Bland-Altman plot: x = mean_D_mm, y = diff_mm."""
    return [dict(branch=s["branch"], station=s["station"], s_mm=s["s_mm"], D_design_mm=s["D_design_mm"], D_print_mm=s["D_print_mm"],
                 mean_D_mm=(s["D_print_mm"] + s["D_design_mm"]) / 2, diff_mm=s["dD_mm"], diff_pct=s["dD_pct"],
                 at_junction=s["at_junction"], fluid_filled=s["fluid_filled"], edited=s["edited"])
            for s in stations if s["excluded"] == ""]


def write_csv(rows, path):
    import csv
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader()
        for r in rows: w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})


def summarise(stations, step):
    """Per branch and for the whole vessel (ALL): n, length, diameter error and centerline error stats."""
    import math
    names = list(dict.fromkeys(s["branch"] for s in stations)) + ["ALL"]
    out = []
    for nm in names:
        v = [s for s in stations if s["excluded"] == "" and (nm == "ALL" or s["branch"] == nm)]
        a = lambda k: np.array([s[k] for s in v], float)
        dD, pct, cl = a("dD_mm"), a("dD_pct"), a("cl_dev_mm")
        f = lambda fn, x: float(fn(x[np.isfinite(x)])) if np.isfinite(x).any() else math.nan
        sd = f(lambda x: np.std(x, ddof=1) if len(x) > 1 else math.nan, dD)
        out.append(dict(branch=nm, n_stations=len(v), length_mm=len(v) * step,
                        dD_mean_mm=f(np.mean, dD), dD_sd_mm=sd, dD_LoA_low_mm=f(np.mean, dD) - 1.96 * sd,
                        dD_LoA_high_mm=f(np.mean, dD) + 1.96 * sd, dD_mean_pct=f(np.mean, pct),
                        cl_dev_mean_mm=f(np.mean, cl), cl_dev_p95_mm=f(lambda x: np.percentile(x, 95), cl),
                        cl_dev_max_mm=f(np.max, cl)))
    return out


def parse_exclude(text):
    """'Aorta: 0-60, 100-110; ICA: 5-12'  ->  {'Aorta': [[0, 60], [100, 110]], 'ICA': [[5, 12]]}"""
    out = {}
    for part in (text or "").split(";"):
        if ":" not in part: continue
        name, ranges = part.split(":", 1)
        for rg in ranges.split(","):
            if "-" not in rg: continue
            lo, hi = rg.strip().split("-", 1)
            out.setdefault(name.strip(), []).append([float(lo), float(hi)])
    return out


def format_exclude(d):
    return "; ".join(f"{k}: " + ", ".join(f"{lo:g}-{hi:g}" for lo, hi in v) for k, v in (d or {}).items())


def read_pth(path, scale=10.0, lps_to_ras=True):
    """SimVascular path file -> (N,3) points in mm, in Slicer RAS by default (SimVascular uses LPS)."""
    import xml.etree.ElementTree as ET
    P = np.array([[float(q.find("pos").get(k)) for k in "xyz"] for q in ET.parse(path).getroot().iter("path_point")]) * scale
    if lps_to_ras: P[:, :2] *= -1
    return P
