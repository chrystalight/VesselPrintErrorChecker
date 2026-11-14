# Vessel Print Error Checker

This is a 3D Slicer module that measures how closely a 3D-printed vessel model matches its design. You give it:
- a CT scan of the printed model;
- the original lumen surface (STL);
- a centerline for each vessel.

It reports two errors along every vessel:
- **Diameter error:** printed minus design lumen diameter, as the equivalent diameter 2√(A/π).
- **Centerline error:** how far the printed vessel's path is from the design. It is reported twice: after one rigid fit for the whole model (this includes any bending of the model), and after a separate fit for each vessel (this leaves only the vessel's own shape).

It was written to validate the SIMONE simulator's printed vessel modules. Their scans and designs are included, so you can run it straight away.

![screenshot](docs/screenshot.png)

---

## Requirements
- [3D Slicer](https://download.slicer.org) 5.x (tested on 5.12.4). Nothing else needs installing.
- [Git LFS](https://git-lfs.com) is not needed. The sample data is stored as ordinary files (about 700 MB).

## Install
1. Clone this repository:
   ```
   git clone <repo-url>
   ```
2. In Slicer, open **Edit → Application Settings → Modules**. Under **Additional module paths**, click **Add** and choose the `VesselPrintErrorChecker/VesselPrintErrorChecker` folder inside the clone (the folder that contains `VesselPrintErrorChecker.py`).
3. Restart Slicer. The module is under **Modules → Quantification → Vessel Print Error Checker**.

## Quick start with the SIMONE data
1. Open the module.
2. Under **SIMONE sample data**, pick *Torso*, *Cranial* or *Femoral* and click **Load sample**. This loads the CT, the design STL and the centerlines, and fills in the inputs.
3. Click **Measure print error**. It takes about 20–60 s.
4. Compare with `SIMONE_validation_files/expected_results/`.

## Using your own data
| Input | How to provide it |
|---|---|
| **Printed CT** | Load it with Slicer's DICOM browser (or any volume format). |
| **Design (STL)** | Drag the STL into Slicer. It must be the lumen surface (inner wall) of the design, in mm. |
| **Centerlines** | Either click **Import SimVascular paths (.pth)...** and pick your path files (set *mm per unit* to 10 if SimVascular worked in cm), or draw one open curve per vessel with Markups → Curve. A curve only needs to run inside the printed lumen; it doesn't have to be centred. |

The centerlines must be drawn on the printed CT, not on the design. The design can be in any position or orientation, because the module registers it automatically.

Tick **Treat bright insert as lumen** under *Parameters* if something bright, such as a tube, sits inside a vessel. Stations where it covers more than 25% of the wall are excluded.

## Outputs
- **Summary table** (per vessel and overall) and **per-station table**. Switch between them with the drop-down, and save both with **Save tables as CSV...**.
- **Diameter error plot** and **centerline error plot** along each vessel. Excluded stations are left as gaps.
- **Inspect cross-section:** choose a vessel and station, or click a point in either plot. The left view reslices the CT onto that cross-section and shows:
  - the printed lumen (cyan);
  - the design lumen mapped into the CT (orange);
  - that station's numbers.

| Column (per-station table) | Meaning |
|---|---|
| `s_mm` | distance along the centerline |
| `D_print_mm`, `D_design_mm` | equivalent lumen diameters |
| `dD_mm`, `dD_pct` | diameter error, print − design |
| `cl_dev_module_mm` | centerline error, whole-model rigid fit |
| `cl_dev_vessel_mm` | centerline error, per-vessel rigid fit (vessels of at least 15 mm) |
| `insert_hidden_frac` | fraction of the wall hidden by an insert |
| `excluded` | empty if the station is included, otherwise the reason (see below) |
| `x, y, z` | printed centerline point (CT world coordinates) |

## Method
1. **Stations:** each centerline is resampled every 0.5 mm.
2. **Printed lumen:** the CT is resliced in a plane normal to the centerline. The wall is found on 90 rays (every 4°) at the **half-maximum** between the lumen level and the wall peak. A ray counts only if the wall is at least 250 HU above the lumen, and rays that disagree with their neighbours are dropped. Each ray is limited to the local lumen radius + 2.5 mm, so it cannot run into a neighbouring vessel.
3. **Re-centring:** each section is re-centred on its lumen centroid and the tangent is re-estimated; this is done 3 times. The printed centerline is therefore measured from the CT, not taken from where the curve was drawn.
4. **Registration:** the printed lumen is rigidly registered to the design with ICP. The search starts from 24 orientations × 9 positions. Each iteration only uses point pairs closer than a threshold that shrinks from 20 mm to 1.5 mm, which makes it robust to diameter bias and unmatched regions.
5. **Design sections:** the design is sectioned normal to its own centerline at the matching stations, after a per-vessel fit, so bending of an elastic model cannot shift which part of the design is compared.
6. **Exclusions:** each excluded station gets one of these reasons in the `excluded` column.

| Reason | Meaning |
|---|---|
| `open` | more than 15% of rays found no wall: an outlet, a junction, or a lumen with no visible contrast |
| `end` | within 2 mm of a centerline end |
| `junction` | another vessel's centerline passes through the section (plus 2 mm either side) |
| `insert contact` | an insert hides more than 25% of the wall |
| `no design section` / `design mismatch` | the design vessel is not found where expected |
| `design artefact` | the design section jumps away from its local trend, as at a capped outlet or a branch ostium |
| `below CT resolution` | the design diameter is under the *Min design diameter* setting (0 = off) |

All parameters are under *Parameters* in the module, with tooltips.

## Command line (no GUI)
`scripts/run_headless.py` runs the same engine with Slicer's Python and writes the CSV tables:
```
<Slicer>/bin/PythonSlicer scripts/run_headless.py                  # all SIMONE samples
<Slicer>/bin/PythonSlicer scripts/run_headless.py --ct DICOM_DIR --design design.stl --paths PTH_DIR --out results
```

## Repository layout
```
VesselPrintErrorChecker/            Slicer module (add this folder to Slicer's module paths)
  VesselPrintErrorChecker.py        user interface, tables, plots, inspect view
  VesselPrintErrorCheckerLib/core.py  measurement engine (numpy/scipy/vtk only; no Slicer dependency)
SIMONE_validation_files/            sample data (see its README)
scripts/run_headless.py             command-line runner
```

## Known limitations
- Lumens below about 2 mm (around 5 voxels at 0.4 mm) are measured narrow because of CT blur, so treat them as approximate.
- Fluid left in a lumen during scanning can hide the wall. Fully filled sections are excluded as `open`, but partly filled sections can read narrow.
- Registration is least certain when the scanned vessels cover only part of the design, as with the SIMONE Cranial module. There the result can shift slightly (a few stations, ~0.03 mm) depending on centerline order.
