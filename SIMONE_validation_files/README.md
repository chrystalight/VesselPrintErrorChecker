# SIMONE validation files

These are CT scans, designs and centerlines for the three printed vessel modules of the SIMONE simulator. Load any of them from the module's **SIMONE sample data** panel, or point the inputs at them yourself.

| Vessel | Material | CT (`CT/`) | Design (`design.stl`) | Branch centerlines (`centerlines/*.pth`) | Notes |
|---|---|---|---|---|---|
| Torso | elastic | 396 slices, 0.465 mm | aortic arch lumen | Aorta, Left/Right Common Carotid, Left/Right Subclavian | water was pooled in the lower aorta during scanning, so ~15–55 mm along *Aorta* reads 1–3 mm narrow (the contour follows the water surface); this stretch is excluded by default (`exclude_mm` in `samples.json`, shown in the module's *Exclude ranges* field) |
| Cranial | rigid | 508 slices, 0.362 mm | carotid lumen | ICA, ECA, OPTH | parts of the lumen held water during scanning; fully filled stretches are measured from the water/wall edge (`fluid_filled` = 1); OPTH is ~1–1.5 mm, near the resolution limit |
| Femoral | elastic | 396 slices, 0.465 mm | aorto-iliac lumen | DescendingAorta, LeftIliac, RightIliac | an introducer sheath lies in each iliac and is too bright for a reliable lumen edge, so both iliacs are left out of the results by default (`LeftIliac: all; RightIliac: all` in *Exclude ranges*). They are still measured, with *Treat bright insert as lumen* (HU ≥ 600) on, because the registration needs them |

## What's in the files
- **CT:** Siemens cone-beam CT (DynaCT), 109 kVp, 512 × 512, in HU. Each module was scanned separately; the simulator was moved between scans but not rotated.
- **Design STL:** the original lumen segmentation, in mm, in its own design coordinate frame.
- **Centerlines:** SimVascular paths drawn on each CT, in cm, in the scanner's (LPS) coordinates. `samples.json` lists each module's files and settings.

## Shared trunks
Some branch paths start in a common trunk and run along it together. The trunk is measured once, in the branch listed first:
- **Cranial:** ICA and ECA both start in the common carotid, which is counted under ECA.
- **Torso:** RightCommonCarotid and RightSubclavian share the brachiocephalic trunk, which is counted under RightCommonCarotid.

## Where Cranial loses stations
- **Water in the lumen:** the common carotid at 8–24 mm, the ECA root at about 33–42 mm and the ICA at about 61–104 mm (apart from two short air pockets) read +10 to +190 HU instead of about −950, so they held water when scanned. The fully filled parts are measured from the water/wall edge: the common carotid at 11.5–21.5 mm and the ICA at 64–83.5 and 98.5–101 mm (67 stations, 33.5 mm). Still lost are the 2 mm either side of each water surface, and the ECA root, where the lumen is too small (about 2 mm) for the water/wall edge to be resolved.
- **The print extends past the design** (`no design section`): the ICA beyond 146 mm, the ECA beyond 68 mm (where the lumen widens to 5.5 mm, an outlet connector) and the OPTH before 27 mm. There is nothing in the design to compare against there.
- **Junctions:** the ICA/ECA bifurcation and the OPTH origin.

The ICA and ECA cross over each other near the bifurcation. Where the drawn ECA path runs through filled lumen there, the module deliberately does not let its centre jump across into the ICA.

## Expected results
`expected_results/` has the summary, per-station and Bland–Altman CSVs from `scripts/run_headless.py` with default settings. Running a module in Slicer should give the same numbers, to within a station or two.

| Vessel | Diameter error, mean ± SD | Centerline error (rigid-body registration) | Stations (of which at junctions) |
|---|---|---|---|
| Torso (Aorta 0–60 mm excluded: water) | −0.49 ± 0.20 mm (−8.1%) | 0.48 mm | 456 (58) |
| Cranial | −0.44 ± 0.14 mm (−14.3%) | 0.61 mm | 239 (37) |
| Femoral (descending aorta only; iliacs not reported) | −0.24 ± 0.17 mm (−2.8%) | 0.45 mm | 87 (0) |

In the Cranial, 67 of the 239 stations are water-filled. They read −0.39 ± 0.14 mm, against −0.46 ± 0.13 mm for the 172 air-filled stations.

If the Femoral iliacs are reported as well (clear the *Exclude ranges* field): left iliac 0.00 ± 0.28 mm (0.38 mm), right iliac +0.01 ± 0.32 mm (4.37 mm), whole vessel −0.09 ± 0.28 mm (1.38 mm). The right iliac's large centerline error is the elastic vessel lying bent relative to the design.

The DICOM headers come from a phantom scan (patient name `TEST^SIMONE`), not a person. They still carry the scanning site's name and the scanner serial number, so remove these before making the repository public if that matters.
