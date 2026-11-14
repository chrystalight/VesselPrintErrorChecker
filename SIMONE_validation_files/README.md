# SIMONE validation files

These are CT scans, designs and centerlines for the three printed vessel modules of the SIMONE simulator. Load any of them from the module's **SIMONE sample data** panel, or point the inputs at them yourself.

| Vessel | Material | CT (`CT/`) | Design (`design.stl`) | Branch centerlines (`centerlines/*.pth`) | Notes |
|---|---|---|---|---|---|
| Torso | elastic | 396 slices, 0.465 mm | aortic arch lumen | Aorta, Left/Right Common Carotid, Left/Right Subclavian | water was pooled in the lower aorta during scanning, so ~15–55 mm along *Aorta* reads 1–3 mm narrow (the contour follows the water surface); this stretch is excluded by default (`exclude_mm` in `samples.json`, shown in the module's *Exclude ranges* field) |
| Cranial | rigid | 508 slices, 0.362 mm | carotid lumen | ICA, ECA, OPTH | some sections held water during scanning (no lumen contrast) and are excluded; OPTH is ~1–1.5 mm, near the resolution limit |
| Femoral | elastic | 396 slices, 0.465 mm | aorto-iliac lumen | DescendingAorta, LeftIliac, RightIliac | a catheter lies in both iliacs, so *Treat bright insert as lumen* (HU ≥ 600) is on; the iliacs are measured with the catheter counted as lumen and are noisier (±0.3 mm) than the aorta |

## What's in the files
- **CT:** Siemens cone-beam CT (DynaCT), 109 kVp, 512 × 512, in HU. Each module was scanned separately; the simulator was moved between scans but not rotated.
- **Design STL:** the original lumen segmentation, in mm, in its own design coordinate frame.
- **Centerlines:** SimVascular paths drawn on each CT, in cm, in the scanner's (LPS) coordinates. `samples.json` lists each module's files and settings.

## Shared trunks
Some branch paths start in a common trunk and run along it together. The trunk is measured once, in the branch listed first:
- **Cranial:** ICA and ECA both start in the common carotid, which is counted under ECA.
- **Torso:** RightCommonCarotid and RightSubclavian share the brachiocephalic trunk, which is counted under RightCommonCarotid.

## Where Cranial loses stations
- **No air in the lumen** (`no lumen near path`, plus 2 mm either side): the common carotid at 8–24 mm, the ECA root at about 33–42 mm, and the ICA at about 61–104 mm (apart from two short air pockets). These read +10 to +190 HU instead of about −950, so they were filled when scanned. Only a dry rescan would recover them.
- **The print extends past the design** (`no design section`): the ICA beyond 146 mm, the ECA beyond 68 mm (where the lumen widens to 5.5 mm, an outlet connector) and the OPTH before 27 mm. There is nothing in the design to compare against there.
- **Junctions:** the ICA/ECA bifurcation and the OPTH origin.

The ICA and ECA cross over each other near the bifurcation. Where the drawn ECA path runs through filled lumen there, the module deliberately does not let its centre jump across into the ICA.

## Expected results
`expected_results/` has the summary and per-station CSVs from `scripts/run_headless.py` with default settings. Running a module in Slicer should give the same numbers, to within a station or two.

| Vessel | Diameter error, mean ± SD | Centerline error (rigid-body registration) | Stations (of which at junctions) |
|---|---|---|---|
| Torso (Aorta 0–60 mm excluded: water) | −0.49 ± 0.20 mm (−8.0%) | 0.51 mm | 467 (69) |
| Cranial | −0.46 ± 0.14 mm (−15.4%) | 0.57 mm | 183 (41) |
| Femoral | −0.07 ± 0.30 mm (−0.6%) | 1.43 mm | 234 (7) |

Per branch in the Femoral: descending aorta −0.24 ± 0.15 mm (0.36 mm), left iliac +0.04 ± 0.32 mm (0.41 mm), right iliac +0.03 ± 0.33 mm (4.73 mm). The right iliac's large centerline error is the elastic vessel lying bent relative to the design; its diameter is unaffected.

The DICOM headers come from a phantom scan (patient name `TEST^SIMONE`), not a person. They still carry the scanning site's name and the scanner serial number, so remove these before making the repository public if that matters.
