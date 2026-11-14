"""
Measurement engine for the Vessel Print Error Checker (no Slicer dependency: numpy, scipy, vtk).

All geometry is in one world frame (Slicer RAS, mm). The CT is given as a numpy array indexed [k, j, i]
plus its 4x4 world->IJK matrix. Centerlines are (N, 3) point arrays in the same world frame. The design
(original lumen) is a vtkPolyData in its own frame; it is registered to the print here.

Per station the engine reports the printed and design equivalent diameters (2*sqrt(A/pi)), their difference,
and the centerline deviation after (a) one rigid fit for the whole module and (b) one rigid fit per vessel.
"""
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from scipy.ndimage import map_coordinates, gaussian_filter1d, binary_dilation, label
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

DEFAULTS = dict(
    station_spacing_mm=0.5,     # distance between cross-sections
    max_search_radius_mm=12.0,  # how far a ray may look for the wall on the first pass
    end_trim_mm=2.0,            # ignore this much at each path end
    junction_margin_mm=2.0,     # extra exclusion either side of a junction
    min_path_fit_mm=15.0,       # shortest vessel that gets its own rigid fit
    mismatch_mm=2.0,            # print/design centres further apart than this (after per-vessel fit) are not compared
    min_diameter_mm=0.0,        # design diameters below this are excluded as below CT resolution (0 = keep all)
    insert_hu=None,             # voxels >= this are an insert (e.g. a tube) and are treated as lumen (None = off)
    insert_dilate_mm=1.0,
    exclude_mm=None,            # {vessel name: [[from_mm, to_mm], ...]} ranges along a centerline to leave out
)
N_RAYS = 90
ANG = np.linspace(0, 2 * np.pi, N_RAYS, endpoint=False)
DR = 0.05            # ray sampling step (mm)
MIN_CONTRAST = 250   # HU between lumen and wall for a ray to count
WALL_SEARCH = 2.5    # mm past the first rise searched for the wall peak
MAX_BAD_RAYS = 0.15  # sections with more rays finding no wall are "open"
MAX_INSERT_RAYS = 0.25
PASSES = 3
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


def kabsch(A, B):
    ca, cb = A.mean(0), B.mean(0); H = (A - ca).T @ (B - cb); U, _, Vt = np.linalg.svd(H)
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
        self.mask = None
        if insert_hu is not None:
            it = max(1, int(round(insert_dilate_mm / spacing_mm)))
            m = binary_dilation(self.v >= insert_hu, iterations=it)
            self.v[m] = np.percentile(self.v[::4, ::4, ::4], 2)
            self.mask = binary_dilation(m, iterations=max(1, int(round(0.5 / spacing_mm)))).astype(np.uint8)

    def _ijk(self, X):
        return X @ self.M[:3, :3].T + self.M[:3, 3]

    def sample(self, X):
        ijk = self._ijk(X)
        return map_coordinates(self.v, [ijk[..., 2], ijk[..., 1], ijk[..., 0]], order=1, cval=np.nan)

    def in_insert(self, X):
        if self.mask is None: return np.zeros(len(X), bool)
        ijk = self._ijk(X)
        return map_coordinates(self.mask, [ijk[..., 2], ijk[..., 1], ijk[..., 0]], order=0, cval=0) > 0


# ----------------------------------------------------------------------------- printed lumen from CT
def lumen_section(ct, c, u, v, r_max):
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
    ins = np.zeros(N_RAYS, bool); hit = np.isfinite(rad)
    if hit.any(): ins[np.nonzero(hit)[0][ct.in_insert(c + rad[hit, None] * dirs[hit])]] = True
    rad[ins] = np.nan
    if np.isfinite(rad).sum() > N_RAYS / 2:      # drop spikes against the neighbouring rays
        med = nan_median_filter(rad, 7, wrap=True)
        rad[np.abs(rad - med) > np.maximum(0.5, 0.15 * med)] = np.nan
    return rad, L, ins


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


def measure_print_path(ct, P0, step, r_search, max_diameter=None):
    """Iteratively re-centred cross-sections along one centerline.
    Ray length: max_diameter/2 + 1 mm if a caliper diameter is given; otherwise r_search on the first pass and
    the local lumen radius + 2.5 mm afterwards (about a wall plus margin; stops rays running into a neighbouring vessel)."""
    C, s = resample(P0, step); n = len(C)
    rmax = np.full(n, max_diameter / 2 + 1.0 if max_diameter else r_search)
    for it in range(PASSES):
        T = tangents(C, step); U, V = frames(T)
        newC, D = C.copy(), np.full(n, np.nan); bad, insert = np.ones(n), np.zeros(n)
        rings = np.full((n, N_RAYS, 3), np.nan)
        for i in range(n):
            rad, L, ins = lumen_section(ct, C[i], U[i], V[i], rmax[i])
            if ct.sample(C[i][None])[0] > L + MIN_CONTRAST / 2:
                c2 = reseed(ct, C[i], U[i], V[i], rmax[i])
                if c2 is None: continue
                C[i] = c2; rad, _, ins = lumen_section(ct, C[i], U[i], V[i], rmax[i])
            bad[i] = np.mean(~np.isfinite(rad) & ~ins); insert[i] = np.mean(ins)
            rf = fill_rays_circle(rad) if ins.any() else fill_rays(rad)
            if rf is None: continue
            xy = np.c_[rf * np.cos(ANG), rf * np.sin(ANG)]; A, cen = polygon(xy)
            D[i] = 2 * np.sqrt(A / np.pi); newC[i] = C[i] + cen[0] * U[i] + cen[1] * V[i]
            rings[i] = C[i] + xy[:, :1] * U[i] + xy[:, 1:] * V[i]
        C = newC
        if not max_diameter and np.isfinite(D).any():
            Rl = nan_median_filter(D / 2, 2 * int(2.5 / step) + 1)
            Rl[~np.isfinite(Rl)] = np.nanmedian(D / 2)
            rmax = np.minimum(r_search, Rl + 2.5)
    return dict(s=s, C=C, T=tangents(C, step), D=D, bad=bad, insert=insert, rings=rings)


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
        flag[(r["s"] < p["end_trim_mm"]) | (r["s"] > r["s"][-1] - p["end_trim_mm"])] = "end"
        for lo, hi in (p.get("exclude_mm") or {}).get(name, []):     # user-chosen ranges, e.g. a water-filled stretch
            flag[(r["s"] >= lo) & (r["s"] <= hi)] = "user excluded"
        junc = np.zeros(len(r["s"]), bool)
        for other, o in pr.items():
            if other == name: continue
            for i in range(len(r["s"])):
                if not np.isfinite(r["rings"][i, 0, 0]): continue
                rel = o["C"] - r["C"][i]; near = np.abs(rel @ r["T"][i]) < 1.0
                if not near.any(): continue
                U, V = frames(r["T"][i:i + 1]); ring = r["rings"][i] - r["C"][i]
                if points_in_polygon(np.c_[ring @ U[0], ring @ V[0]], np.c_[rel[near] @ U[0], rel[near] @ V[0]]).any(): junc[i] = True
        if junc.any():
            junc = np.convolve(junc, np.ones(2 * int(p["junction_margin_mm"] / step) + 1), "same") > 0
            flag[junc & (flag == "")] = "junction"
        r["flag"] = flag

    say(0.6, "Registering print to design")
    wall = np.vstack([r["rings"][np.isin(r["flag"], ["", "junction", "insert contact"])].reshape(-1, 3) for r in pr.values()])
    wall = wall[np.isfinite(wall).all(1)]
    design_pts = vtk_to_numpy(design_poly.GetPoints().GetData()).astype(float)
    R, t, reg = register(wall, design_pts, lambda f: say(0.6 + 0.2 * f, "Registering print to design"))

    stations, paths = [], {}
    for k, name in enumerate(names):
        say(0.8 + 0.2 * k / len(names), f"Sectioning design: {name}")
        r = pr[name]; use = r["flag"] == ""
        Cp, Tp = r["C"] @ R.T + t, r["T"] @ R.T
        sub = crop(design_poly, Cp, 2 * r["dmax"])
        Co, _, _ = measure_design_path(sub, Cp, Tp, r["dmax"], step, passes=1)
        Rp, tp = np.eye(3), np.zeros(3); okp = use & np.isfinite(Co).all(1)
        if (okp.sum() - 1) * step >= p["min_path_fit_mm"]:
            Rp, tp = icp(Cp[okp], densify(Co[okp]), iters=50, keep=0.8)
        Cq = Cp @ Rp.T + tp
        Co, Do, loops = measure_design_path(sub, Cq, Tp @ Rp.T, r["dmax"], step)
        mismatch = ~(np.linalg.norm(Co - Cq, axis=1) <= p["mismatch_mm"])
        ok = use & np.isfinite(Do) & ~mismatch
        dev_mod, dev_path = np.full(len(Cp), np.nan), np.full(len(Cp), np.nan)
        if ok.sum() > 5:
            dense = densify(Co[ok]); dev_mod[ok] = cKDTree(dense).query(Cp[ok])[0]
            if (ok.sum() - 1) * step >= p["min_path_fit_mm"]:
                Rp, tp = icp(Cp[ok], dense, Rp, tp, iters=50, keep=1.0)
                dev_path[ok] = cKDTree(dense).query(Cp[ok] @ Rp.T + tp)[0]
        mismatch |= dev_path > p["mismatch_mm"]
        small = Do < (p["min_diameter_mm"] or 0)
        trend = nan_median_filter(Do, 2 * int(3 / step) + 1)
        artefact = np.abs(Do - trend) > np.maximum(0.3, 0.1 * trend)
        ok &= ~mismatch & ~small & ~artefact
        Rc, tc = Rp @ R, Rp @ t + tp                      # print -> design for this vessel
        for i in range(len(Cp)):
            why = r["flag"][i] or ("no design section" if not np.isfinite(Do[i]) else "design mismatch" if mismatch[i]
                                   else "below CT resolution" if small[i] else "design artefact" if artefact[i] else "")
            stations.append(dict(path=name, station=i, s_mm=float(r["s"][i]), excluded=why,
                                 D_print_mm=float(r["D"][i]), D_design_mm=float(Do[i]), dD_mm=float(r["D"][i] - Do[i]),
                                 dD_pct=float(100 * (r["D"][i] - Do[i]) / Do[i]) if Do[i] > 0 else np.nan,
                                 cl_dev_module_mm=float(dev_mod[i]), cl_dev_vessel_mm=float(dev_path[i]),
                                 insert_hidden_frac=float(r["insert"][i]),
                                 x=float(r["C"][i, 0]), y=float(r["C"][i, 1]), z=float(r["C"][i, 2])))
        paths[name] = dict(s=r["s"], C=r["C"], T=r["T"], rings=r["rings"], ok=ok,
                           design_loops_in_ct=[None if L is None else (L - tc) @ Rc for L in loops])
    say(1.0, "Done")
    return dict(stations=stations, paths=paths, R=R, t=t, registration=reg)


def summarise(stations, step):
    """Per vessel and for all vessels: n, length, diameter error stats, centerline error stats (valid stations)."""
    import math
    names = list(dict.fromkeys(s["path"] for s in stations)) + ["ALL"]
    out = []
    for nm in names:
        v = [s for s in stations if s["excluded"] == "" and (nm == "ALL" or s["path"] == nm)]
        a = lambda k: np.array([s[k] for s in v], float)
        dD, pct, cm, cv = a("dD_mm"), a("dD_pct"), a("cl_dev_module_mm"), a("cl_dev_vessel_mm")
        f = lambda fn, x: float(fn(x[np.isfinite(x)])) if np.isfinite(x).any() else math.nan
        sd = f(lambda x: np.std(x, ddof=1) if len(x) > 1 else math.nan, dD)
        out.append(dict(vessel=nm, n_stations=len(v), length_mm=len(v) * step,
                        dD_mean_mm=f(np.mean, dD), dD_sd_mm=sd, dD_LoA_low_mm=f(np.mean, dD) - 1.96 * sd,
                        dD_LoA_high_mm=f(np.mean, dD) + 1.96 * sd, dD_mean_pct=f(np.mean, pct),
                        cl_module_mean_mm=f(np.mean, cm), cl_module_p95_mm=f(lambda x: np.percentile(x, 95), cm),
                        cl_module_max_mm=f(np.max, cm), cl_vessel_mean_mm=f(np.mean, cv),
                        cl_vessel_p95_mm=f(lambda x: np.percentile(x, 95), cv), cl_vessel_max_mm=f(np.max, cv)))
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
