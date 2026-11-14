# Vessel Print Error Checker

This is a 3D Slicer module that measures how closely a 3D-printed vessel model matches its design. You give it:
- a CT scan of the printed model;
- the original lumen surface (STL);
- a centerline for each branch.

A **vessel** here is one connected vessel tree (one CT and one design STL), for example the SIMONE Torso: the aorta plus its four branches. A **branch** is one centerline within it.

It reports two errors along every branch, plus a summary for the whole vessel:
- **Diameter error:** printed minus design lumen diameter, as the equivalent diameter 2√(A/π).
- **Centerline error:** the distance between the printed and design centerlines after the branch is rigidly aligned on its own. This measures the branch's printed shape. It is not affected by how a soft, elastic print happened to lie in the scanner.

It was written to validate the SIMONE simulator's printed vessel modules. Their scans and designs are included, so you can run it straight away.

![screenshot](docs/screenshot.png)

---

## Install (about 5 minutes, once)

You need 3D Slicer and a copy of this repository. There is nothing else to install: no Python packages, no Slicer extensions, and no build step.

### Step 1: Install 3D Slicer
Download and install **3D Slicer 5.x** from <https://download.slicer.org>. It was tested on 5.12.4; any recent stable release should work.

### Step 2: Download this repository
Pick **one** of these:
- **Without git:** on <https://github.com/chrystalight/VesselPrintErrorChecker>, click the green **Code** button, then **Download ZIP**. Unzip it somewhere permanent, for example `Documents\VesselPrintErrorChecker`, not your Downloads folder, because Slicer will load the module from there every time it starts.
- **With git:**
  ```
  git clone https://github.com/chrystalight/VesselPrintErrorChecker.git
  ```

The download is about 700 MB, almost all of it the SIMONE sample CT scans.

You should now have this structure:
```
VesselPrintErrorChecker\                     <- repository folder
├── VesselPrintErrorChecker\                 <- THE MODULE FOLDER (you'll point Slicer here)
│   ├── VesselPrintErrorChecker.py
│   └── VesselPrintErrorCheckerLib\
├── SIMONE_validation_files\
├── scripts\
└── README.md
```

### Step 3: Tell Slicer where the module is
1. Open 3D Slicer.
2. Go to **Edit → Application Settings** (on macOS: **Slicer → Preferences**).
3. Click **Modules** in the list on the left.
4. Find the box **Additional module paths** and click **Add** next to it.
5. Browse to the **inner** `VesselPrintErrorChecker` folder, the one that contains `VesselPrintErrorChecker.py` (marked *THE MODULE FOLDER* above), and select it.
6. Click **OK**. When Slicer asks to restart, click **Yes**, or close and reopen Slicer yourself.

> **Shortcut:** instead of steps 2–6, drag `VesselPrintErrorChecker.py` from your file browser onto the Slicer window. Choose **Add Python scripted modules to the application**, tick the module, click **Yes** to add it to the module paths, then restart Slicer.

### Step 4: Check it worked
1. Click the **magnifying glass** next to the *Modules* drop-down on the toolbar (or press **Ctrl+F**).
2. Type `Vessel Print` and open **Vessel Print Error Checker**. It also appears under **Modules → Quantification**.
3. The left panel should show these sections: **SIMONE sample data**, **Inputs**, **Parameters**, **Measure print error** and **Inspect cross-section**.

### Troubleshooting
| Problem | Fix |
|---|---|
| The module isn't in the list | You probably added the outer repository folder. Remove it from *Additional module paths*, add the inner `VesselPrintErrorChecker` folder (the one with the `.py` file) and restart. |
| It's listed but shows an error | Open **View → Python Console** to see the message. The folder must contain `VesselPrintErrorCheckerLib\` next to the `.py` file; don't move the `.py` on its own. |
| The SIMONE sample list is empty | `SIMONE_validation_files` must sit next to the module folder, as in the repository. Keep the repository layout as it is. |
| You moved or renamed the folder | Update the path in *Additional module paths* and restart. |

### Updating
Run `git pull`, or download the ZIP again into the same place, then restart Slicer.

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
| **Centerlines** | Either click **Import SimVascular paths (.pth)...** and pick your path files (set *mm per unit* to 10 if SimVascular worked in cm), or draw one open curve per branch with Markups → Curve. A curve only needs to run inside the printed lumen; it doesn't have to be centred. |

The centerlines must be drawn on the printed CT, not on the design. The design can be in any position or orientation, because the module registers it automatically.

Tick **Treat bright insert as lumen** under *Parameters* if something bright, such as a tube, sits inside a vessel. Stations where it covers more than 25% of the wall are excluded.

To leave out a stretch you know is unreliable, such as a section that held water during scanning, fill in **Exclude ranges (mm)**. Give distances along the centerline, for example `Aorta: 0-60; ICA: 5-12`; several ranges per vessel are separated by commas. Excluded stretches are left out of the results and also out of the registration. They show as `user excluded` in the per-station table. The SIMONE Torso sample fills this in for you (`Aorta: 0-60`, where water pooled in the lower aorta).

## Outputs
- **Summary table** (each branch, plus `ALL` for the whole vessel) and **per-station table**. Switch between them with the drop-down, and save both with **Save tables as CSV...**.
- **Diameter error plot** and **centerline error plot** for the branch you are inspecting. Pick another branch in *Inspect cross-section* to switch.
  - Both plots run along the branch's centerline, and each title gives the branch mean.
  - A grey vertical line marks the inspected station. Excluded stations are left as gaps.
  - Use the summary table to compare branches.
- **Inspect cross-section:** choose a branch and station, or click a point in either plot. The left view reslices the CT onto that cross-section and shows:
  - the printed lumen (cyan);
  - the design lumen mapped into the CT (orange);
  - that station's numbers.

  The view re-centres and zooms onto the contours every time you change station, even if you have panned or zoomed in the meantime.
- **3D view:** tick **Show 3D view** under *3D view*. A 3D panel opens next to the cross-section. It shows the measured printed lumen over a translucent grey copy of the design, registered into the CT.
  - **Colour printed lumen by:** choose *Diameter error* (blue = narrower than designed, red = wider), *Centerline error*, or *None*. A colour legend appears in the 3D view.
  - **Colour scale limit:** sets the end of the colour scale. Diameter error runs from −limit to +limit; centerline error runs from 0 to the limit.
  - **Centre 3D view on the inspected station:** keeps the 3D camera on the station you are inspecting, without changing your zoom or angle.
  - Only included stations are drawn. Gaps in the surface are excluded stations.
  - The original design model is hidden in 3D because it sits in its own coordinate frame, away from the CT. The registered copy is called *VPE design (registered to CT)*.

| Column (per-station table) | Meaning |
|---|---|
| `s_mm` | distance along the centerline |
| `D_print_mm`, `D_design_mm` | equivalent lumen diameters |
| `dD_mm`, `dD_pct` | diameter error, print − design |
| `branch` | which centerline the station belongs to |
| `cl_dev_mm` | centerline error: distance to the design centerline after aligning this branch on its own (branches with at least 15 mm included) |
| `cl_dev_whole_vessel_mm` | for reference only: the same distance after one alignment of the whole vessel. For elastic prints it includes how the vessel lay during the scan. |
| `insert_hidden_frac` | fraction of the wall hidden by an insert |
| `excluded` | empty if the station is included, otherwise the reason (see below) |
| `x, y, z` | printed centerline point (CT world coordinates) |

## Method
1. **Stations:** each centerline is resampled every 0.5 mm.
2. **Printed lumen:** the CT is resliced in a plane normal to the centerline. The wall is found on 90 rays (every 4°) at the **half-maximum** between the lumen level and the wall peak. A ray counts only if the wall is at least 250 HU above the lumen, and rays that disagree with their neighbours are dropped. Each ray is limited to the local lumen radius + 2.5 mm, so it cannot run into a neighbouring vessel.
3. **Re-centring:** each section is re-centred on its lumen centroid and the tangent is re-estimated; this is done 3 times. The printed centerline is therefore measured from the CT, not taken from where the curve was drawn.
4. **Registration:** the printed lumen is rigidly registered to the design with ICP. The search starts from 24 orientations × 9 positions. Each iteration only uses point pairs closer than a threshold that shrinks from 20 mm to 1.5 mm, which makes it robust to diameter bias and unmatched regions.
5. **Design sections:** each branch is then aligned on its own, and the design is sectioned normal to its own centerline at the matching stations. This way, bending of an elastic print cannot shift which part of the design is compared.
6. **Errors:** diameter error is printed minus design equivalent diameter. Centerline error is the distance from the printed to the design centerline after aligning the branch on its own (rigid ICP).
7. **Exclusions:** each excluded station gets one of these reasons in the `excluded` column.

| Reason | Meaning |
|---|---|
| `open` | more than 15% of rays found no wall: an outlet, a junction, or a lumen with no visible contrast |
| `end` | within 2 mm of a centerline end |
| `junction` | another branch leaves here: its centerline passes through the section at an angle, or only briefly alongside where it starts inside this one (plus 2 mm either side) |
| `shared with X` | this stretch runs along the same lumen as branch X (for example, two branch paths that both start in a common trunk). It is measured once, in the branch listed first (X), so the trunk isn't counted twice. |
| `insert contact` | an insert hides more than 25% of the wall |
| `no design section` / `design mismatch` | the design branch is not found where expected |
| `design artefact` | the design section jumps away from its local trend, as at a capped outlet or a branch ostium |
| `below CT resolution` | the design diameter is under the *Min design diameter* setting (0 = off) |
| `user excluded` | inside a range given in *Exclude ranges (mm)* |

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
- Fluid left in a lumen during scanning can hide the wall. Fully filled sections are excluded as `open`, but partly filled sections read narrow, because the contour follows the air–fluid surface instead of the wall. In the inspect view the contour then has a flat side. Scan the prints dry, or use *Treat bright insert as lumen* only for solid inserts, not fluid.
- Registration is least certain when the scanned vessels cover only part of the design, as with the SIMONE Cranial module. There the result can shift slightly (a few stations, ~0.03 mm) depending on centerline order.
