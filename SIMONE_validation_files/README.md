# SIMONE validation files

These are CT scans, designs and centerlines for the three printed vessel modules of the SIMONE simulator. Load any of them from the module's **SIMONE sample data** panel, or point the inputs at them yourself.

| Module | Material | CT (`CT/`) | Design (`design.stl`) | Centerlines (`centerlines/*.pth`) | Notes |
|---|---|---|---|---|---|
| Torso | elastic | 396 slices, 0.465 mm | aortic arch lumen | Aorta, Left/Right Common Carotid, Left/Right Subclavian | water was pooled in the lower aorta during scanning, so ~15–55 mm along *Aorta* reads 1–3 mm narrow (the contour follows the water surface); this stretch is excluded by default (`exclude_mm` in `samples.json`, shown in the module's *Exclude ranges* field) |
| Cranial | rigid | 508 slices, 0.362 mm | carotid lumen | ICA, ECA, OPTH | some sections held water during scanning (no lumen contrast) and are excluded; OPTH is ~1–1.5 mm, near the resolution limit |
| Femoral | elastic | 396 slices, 0.465 mm | aorto-iliac lumen | DescendingAorta, LeftIliac, RightIliac | a tube sits in both iliacs, so use *Treat bright insert as lumen* (HU ≥ 600); the iliacs end up fully excluded because the tube covers about half the wall |

## What's in the files
- **CT:** Siemens cone-beam CT (DynaCT), 109 kVp, 512 × 512, in HU. Each module was scanned separately; the simulator was moved between scans but not rotated.
- **Design STL:** the original lumen segmentation, in mm, in its own design coordinate frame.
- **Centerlines:** SimVascular paths drawn on each CT, in cm, in the scanner's (LPS) coordinates. `samples.json` lists each module's files and settings.

## Expected results
`expected_results/` has the summary and per-station CSVs from `scripts/run_headless.py` with default settings. Running a module in Slicer should give the same numbers, to within a station or two.

| Module | Diameter error, mean ± SD | Centerline error, whole module | Centerline error, per vessel |
|---|---|---|---|
| Torso (Aorta 0–60 mm excluded: water) | −0.51 ± 0.21 mm (−8.3%) | 0.51 mm | 0.18 mm |
| Cranial | −0.53 ± 0.24 mm (−17.8%) | 0.74 mm | 0.22 mm |
| Femoral (descending aorta) | −0.25 ± 0.19 mm (−3.0%) | 0.51 mm | 0.17 mm |

The DICOM headers come from a phantom scan (patient name `TEST^SIMONE`), not a person. They still carry the scanning site's name and the scanner serial number, so remove these before making the repository public if that matters.
