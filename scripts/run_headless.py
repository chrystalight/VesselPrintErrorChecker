"""
Run the Vessel Print Error Checker without the Slicer GUI and write CSV tables.
Uses Slicer's bundled Python (numpy, scipy, pydicom, vtk); any Python with those packages also works.

  PythonSlicer scripts/run_headless.py                       # all SIMONE samples -> SIMONE_validation_files/expected_results
  PythonSlicer scripts/run_headless.py Femoral               # one SIMONE sample
  PythonSlicer scripts/run_headless.py --ct DICOM_DIR --design design.stl --paths PTH_DIR [--insert-hu 600] --out OUT_DIR

PythonSlicer is in the Slicer install folder: <Slicer>/bin/PythonSlicer(.exe).
"""
import argparse, csv, glob, json, os, sys
import numpy as np, pydicom, vtk

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "VesselPrintErrorChecker"))
from VesselPrintErrorCheckerLib import core  # noqa: E402


def load_dicom(folder):
    """DICOM series -> (array [k,j,i], world(LPS)->IJK matrix, voxel size). Works in the scanner (LPS) frame,
    which is also the frame of SimVascular paths, so no conversion is needed here."""
    ds = [pydicom.dcmread(f) for f in sorted(glob.glob(os.path.join(folder, "*"))) if os.path.isfile(f)]
    iop = np.array(ds[0].ImageOrientationPatient, float); nrm = np.cross(iop[:3], iop[3:])
    ds.sort(key=lambda d: float(np.dot(d.ImagePositionPatient, nrm)))
    P = np.array([d.ImagePositionPatient for d in ds], float); px = np.array(ds[0].PixelSpacing, float)
    vol = np.stack([d.pixel_array * float(getattr(d, "RescaleSlope", 1)) + float(getattr(d, "RescaleIntercept", 0)) for d in ds])
    A = np.eye(4); A[:3, :3] = np.c_[iop[:3] * px[1], iop[3:] * px[0], (P[-1] - P[0]) / (len(P) - 1)]; A[:3, 3] = P[0]
    return vol, np.linalg.inv(A), float(min(px))


def write_csv(rows, path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader()
        for r in rows: w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})


def run(name, ct_dir, design, paths_dir, mm_per_unit, insert_hu, out, params, exclude=None):
    print(f"== {name}")
    vol, M, px = load_dicom(ct_dir)
    ct = core.CTVolume(vol, M, insert_hu, 1.0, px)
    rd = vtk.vtkSTLReader(); rd.SetFileName(design); rd.Update()
    curves = {os.path.splitext(os.path.basename(f))[0]: core.read_pth(f, mm_per_unit, lps_to_ras=False)
              for f in sorted(glob.glob(os.path.join(paths_dir, "*.pth")))}
    if exclude: print("   excluding:", core.format_exclude(exclude))
    res = core.run(ct, curves, rd.GetOutput(), dict(params, exclude_mm=exclude), lambda f, m: print(f"   {100 * f:3.0f}%  {m}", end="\r"))
    summ = core.summarise(res["stations"], params.get("station_spacing_mm", 0.5))
    os.makedirs(out, exist_ok=True)
    write_csv(res["stations"], os.path.join(out, f"{name}_stations.csv"))
    write_csv(summ, os.path.join(out, f"{name}_summary.csv"))
    print(f"\n   registration: " + ", ".join(f"{k} {v:.2f}" for k, v in res["registration"].items()))
    for r in summ:
        print(f"   {r['vessel']:20s} n={r['n_stations']:4d}  dD {r['dD_mean_mm']:+.2f} ± {r['dD_sd_mm']:.2f} mm ({r['dD_mean_pct']:+.1f}%)  "
              f"centerline {r['cl_module_mean_mm']:.2f} (module) / {r['cl_vessel_mean_mm']:.2f} (vessel) mm")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("samples", nargs="*", help="SIMONE sample names (default: all)")
    ap.add_argument("--ct"); ap.add_argument("--design"); ap.add_argument("--paths")
    ap.add_argument("--mm-per-unit", type=float, default=10.0); ap.add_argument("--insert-hu", type=float)
    ap.add_argument("--out"); ap.add_argument("--name", default="custom")
    ap.add_argument("--exclude", help='ranges to leave out, e.g. "Aorta: 0-60; ICA: 5-12"')
    for k, v in core.DEFAULTS.items():
        if k not in ("insert_hu", "insert_dilate_mm", "exclude_mm"): ap.add_argument("--" + k.replace("_", "-"), type=float, default=v)
    a = ap.parse_args()
    params = {k: getattr(a, k) for k in core.DEFAULTS if k not in ("insert_hu", "insert_dilate_mm", "exclude_mm")}
    if a.ct:
        run(a.name, a.ct, a.design, a.paths, a.mm_per_unit, a.insert_hu, a.out or "results", params, core.parse_exclude(a.exclude))
    else:
        S = os.path.join(ROOT, "SIMONE_validation_files"); cfg = json.load(open(os.path.join(S, "samples.json")))
        for name in a.samples or cfg:
            c = cfg[name]
            run(name, os.path.join(S, c["ct_folder"]), os.path.join(S, c["design"]), os.path.join(S, c["centerlines"]),
                c["mm_per_unit"], c["insert_hu"], a.out or os.path.join(S, "expected_results"), params, c.get("exclude_mm"))
