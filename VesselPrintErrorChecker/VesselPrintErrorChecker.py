import os
import sys

import numpy as np
import qt
import ctk
import vtk
import slicer
from slicer.ScriptedLoadableModule import ScriptedLoadableModule, ScriptedLoadableModuleWidget, ScriptedLoadableModuleLogic

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from VesselPrintErrorCheckerLib import core  # noqa: E402

SAMPLES_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "SIMONE_validation_files")
LAYOUT_ID = 50731
LAYOUT_XML = """
<layout type="horizontal" split="true">
 <item splitSize="450">
  <layout type="vertical" split="true">
   <item splitSize="600">
    <view class="vtkMRMLSliceNode" singletontag="Red">
     <property name="orientation" action="default">Axial</property>
     <property name="viewlabel" action="default">R</property>
     <property name="viewcolor" action="default">#F34A33</property>
    </view>
   </item>
   <item splitSize="300"><view class="vtkMRMLTableViewNode" singletontag="VPETable"><property name="viewlabel" action="default">T</property></view></item>
  </layout>
 </item>
 <item splitSize="650">
  <layout type="vertical" split="true">
   <item><view class="vtkMRMLPlotViewNode" singletontag="VPEPlotDiameter"><property name="viewlabel" action="default">D</property></view></item>
   <item><view class="vtkMRMLPlotViewNode" singletontag="VPEPlotCenterline"><property name="viewlabel" action="default">C</property></view></item>
  </layout>
 </item>
</layout>"""
LAYOUT_3D_ID = 50732
LAYOUT_3D_XML = LAYOUT_XML.replace(
    '<view class="vtkMRMLSliceNode" singletontag="Red">',
    '<layout type="horizontal" split="true"><item><view class="vtkMRMLSliceNode" singletontag="Red">', 1).replace(
    '    </view>\n   </item>\n   <item splitSize="300">',
    '    </view></item><item><view class="vtkMRMLViewNode" singletontag="1"><property name="viewlabel" action="default">1</property></view></item></layout>\n   </item>\n   <item splitSize="300">', 1)
COLOR_BY = [("None", None), ("Diameter error", "dD_mm"), ("Centerline error", "cl_dev_mm")]
LEGEND_TITLE = {"dD_mm": "Diameter error (mm)", "cl_dev_mm": "Centerline error (mm)"}
COLORS = [(0.16, 0.47, 0.84), (0.77, 0.43, 0.11), (0.09, 0.57, 0.50), (0.42, 0.31, 0.77), (0.80, 0.25, 0.45),
          (0.45, 0.60, 0.15), (0.20, 0.65, 0.85), (0.55, 0.40, 0.25)]


class VesselPrintErrorChecker(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        parent.title = "Vessel Print Error Checker"
        parent.categories = ["Quantification"]
        parent.contributors = ["Capstone DMD team"]
        parent.helpText = (
            "Compares a 3D-printed vessel model, from its CT scan, with the original design (STL).<br>"
            "For each centerline the CT is resliced in planes normal to the centerline every 0.5 mm; the lumen wall is found "
            "at the half-maximum between lumen and wall on 90 rays and the section is re-centred (3 passes). The print is "
            "rigidly registered to the design, which is sectioned at the matching locations.<br>"
            "Outputs: per-station table, diameter error and centerline error plots, and a click-to-inspect cross-section view.")
        parent.acknowledgementText = ""


# ----------------------------------------------------------------------------------------------------------- widget
class VesselPrintErrorCheckerWidget(ScriptedLoadableModuleWidget):
    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)
        self.logic = VesselPrintErrorCheckerLogic()

        # --- SIMONE sample data
        sbox = ctk.ctkCollapsibleButton(); sbox.text = "SIMONE sample data"; self.layout.addWidget(sbox)
        srow = qt.QHBoxLayout(sbox)
        self.sampleCombo = qt.QComboBox(); self.sampleCombo.addItems(self.logic.sampleNames())
        self.sampleButton = qt.QPushButton("Load sample"); self.sampleButton.enabled = self.sampleCombo.count > 0
        self.sampleButton.setToolTip("Loads the CT, design STL and centerlines of one SIMONE module and fills in the inputs")
        srow.addWidget(self.sampleCombo); srow.addWidget(self.sampleButton)

        # --- inputs
        box = ctk.ctkCollapsibleButton(); box.text = "Inputs"; self.layout.addWidget(box)
        form = qt.QFormLayout(box)
        self.ctSelector = self._combo(["vtkMRMLScalarVolumeNode"], "CT of the printed model")
        form.addRow("Printed CT:", self.ctSelector)
        self.designSelector = self._combo(["vtkMRMLModelNode"], "Original lumen surface (e.g. the STL)")
        form.addRow("Design (STL):", self.designSelector)
        self.curvesSelector = slicer.qMRMLCheckableNodeComboBox()
        self.curvesSelector.nodeTypes = ["vtkMRMLMarkupsCurveNode"]
        self.curvesSelector.setMRMLScene(slicer.mrmlScene)
        self.curvesSelector.setToolTip("One open curve per vessel, drawn in the lumen of the printed CT")
        form.addRow("Centerlines:", self.curvesSelector)
        row = qt.QHBoxLayout()
        self.importButton = qt.QPushButton("Import SimVascular paths (.pth)...")
        self.unitsSpin = qt.QDoubleSpinBox(); self.unitsSpin.setRange(0.001, 1000); self.unitsSpin.value = 10.0
        self.unitsSpin.setToolTip("mm per SimVascular unit (10 if SimVascular worked in cm)")
        row.addWidget(self.importButton); row.addWidget(qt.QLabel("mm per unit:")); row.addWidget(self.unitsSpin)
        form.addRow(row)
        self.excludeEdit = qt.QLineEdit(); self.excludeEdit.placeholderText = "optional, e.g.  Aorta: 0-60; ICA: 5-12"
        self.excludeEdit.setToolTip("Stretches of a vessel to leave out, in mm along its centerline (e.g. a section that held "
                                    "water during scanning). Format: Branch: from-to, from-to; Branch2: from-to")
        form.addRow("Exclude ranges (mm):", self.excludeEdit)

        # --- parameters
        pbox = ctk.ctkCollapsibleButton(); pbox.text = "Parameters"; pbox.collapsed = True; self.layout.addWidget(pbox)
        pform = qt.QFormLayout(pbox)
        d = core.DEFAULTS
        self.spacing = self._spin(0.1, 5, d["station_spacing_mm"], " mm", "Distance between cross-sections")
        self.search = self._spin(2, 50, d["max_search_radius_mm"], " mm", "Largest distance from the centerline searched for the wall (first pass)")
        self.minDiam = self._spin(0, 10, d["min_diameter_mm"], " mm", "Design diameters below this are excluded as below CT resolution (0 = keep all)")
        self.junction = self._spin(0, 20, d["junction_margin_mm"], " mm", "Extra exclusion either side of a vessel junction")
        self.trim = self._spin(0, 20, d["end_trim_mm"], " mm", "Ignore this much at each end of a centerline")
        self.mismatch = self._spin(0.5, 20, d["mismatch_mm"], " mm", "Print and design centres further apart than this are treated as a failed match and not compared")
        self.insertCheck = qt.QCheckBox("Treat bright insert as lumen, HU ≥")
        self.insertHU = self._spin(100, 5000, 600, "", "Voxels at or above this value are an insert (e.g. a tube) inside the lumen")
        irow = qt.QHBoxLayout(); irow.addWidget(self.insertCheck); irow.addWidget(self.insertHU)
        for lab, w in [("Station spacing:", self.spacing), ("Max search radius:", self.search), ("Min design diameter:", self.minDiam),
                       ("Junction margin:", self.junction), ("End trim:", self.trim), ("Mismatch limit:", self.mismatch)]:
            pform.addRow(lab, w)
        pform.addRow(irow)

        # --- run
        self.applyButton = qt.QPushButton("Measure print error"); self.applyButton.setMinimumHeight(30)
        self.layout.addWidget(self.applyButton)
        self.progress = qt.QProgressBar(); self.progress.setRange(0, 100); self.progress.visible = False
        self.layout.addWidget(self.progress)
        self.status = qt.QLabel(""); self.status.wordWrap = True; self.layout.addWidget(self.status)

        # --- inspect
        ibox = ctk.ctkCollapsibleButton(); ibox.text = "Inspect cross-section"; self.layout.addWidget(ibox)
        iform = qt.QFormLayout(ibox)
        self.vesselCombo = qt.QComboBox(); iform.addRow("Branch:", self.vesselCombo)
        self.stationSlider = ctk.ctkSliderWidget(); self.stationSlider.decimals = 0; self.stationSlider.singleStep = 1
        self.stationSlider.minimum = 0; self.stationSlider.maximum = 0
        iform.addRow("Station:", self.stationSlider)
        self.info = qt.QLabel("Run a measurement, then pick a branch and station here, or click a point in either plot. "
                              "The plots show the selected branch; the grey line marks the inspected station.")
        self.info.wordWrap = True; iform.addRow(self.info)
        legend = qt.QLabel("<span style='color:#22d3ee'>■</span> printed lumen &nbsp; "
                           "<span style='color:#f59e0b'>■</span> design lumen (mapped into the CT)")
        iform.addRow(legend)
        self.editCheck = qt.QCheckBox("Edit printed contour (drag its points in the cross-section view)")
        iform.addRow(self.editCheck)
        erow = qt.QHBoxLayout()
        self.applyEditButton = qt.QPushButton("Apply edited contour"); self.applyEditButton.enabled = False
        self.applyEditButton.setToolTip("Use the edited outline for this station: its diameter, error and centre are recomputed")
        self.toggleButton = qt.QPushButton("Exclude / include station")
        self.toggleButton.setToolTip("Leave this station out of the results, or bring back an excluded one that has both contours")
        erow.addWidget(self.applyEditButton); erow.addWidget(self.toggleButton); iform.addRow(erow)

        # --- 3D view
        vbox = ctk.ctkCollapsibleButton(); vbox.text = "3D view"; self.layout.addWidget(vbox)
        vform = qt.QFormLayout(vbox)
        self.show3D = qt.QCheckBox("Show 3D view (printed lumen over the registered design)")
        vform.addRow(self.show3D)
        self.colorBy = qt.QComboBox(); self.colorBy.addItems([c[0] for c in COLOR_BY]); self.colorBy.currentIndex = 1
        vform.addRow("Colour printed lumen by:", self.colorBy)
        self.colorMax = self._spin(0.1, 10, 1.0, " mm", "Colour scale limit: diameter error -limit..+limit, centerline error 0..limit")
        vform.addRow("Colour scale limit:", self.colorMax)
        self.follow3D = qt.QCheckBox("Centre 3D view on the inspected station"); self.follow3D.checked = True
        vform.addRow(self.follow3D)

        rrow = qt.QHBoxLayout()
        self.tableChoice = qt.QComboBox(); self.tableChoice.addItems(["Summary table", "Per-station table"])
        self.saveButton = qt.QPushButton("Export CSV..."); self.saveButton.enabled = False
        rrow.addWidget(self.tableChoice); rrow.addWidget(self.saveButton); self.layout.addLayout(rrow)
        self.layout.addStretch(1)

        self.importButton.connect("clicked()", self.onImport)
        self.sampleButton.connect("clicked()", self.onLoadSample)
        self.applyButton.connect("clicked()", self.onApply)
        self.vesselCombo.connect("currentIndexChanged(int)", self.onVesselChanged)
        self.stationSlider.connect("valueChanged(double)", self.onStationChanged)
        self.saveButton.connect("clicked()", self.onSave)
        self.tableChoice.connect("currentIndexChanged(int)", self.onTableChoice)
        self.show3D.connect("toggled(bool)", self.on3DOptions)
        self.editCheck.connect("toggled(bool)", self.onEditToggled)
        self.applyEditButton.connect("clicked()", self.onApplyEdit)
        self.toggleButton.connect("clicked()", self.onToggleStation)
        self.colorBy.connect("currentIndexChanged(int)", self.on3DOptions)
        self.colorMax.connect("valueChanged(double)", self.on3DOptions)
        self._plotConnected = []

    def _combo(self, types, tip):
        c = slicer.qMRMLNodeComboBox(); c.nodeTypes = types; c.selectNodeUponCreation = True
        c.addEnabled = False; c.removeEnabled = False; c.noneEnabled = True; c.showHidden = False
        c.setMRMLScene(slicer.mrmlScene); c.setToolTip(tip); return c

    def _spin(self, lo, hi, val, suffix, tip):
        s = qt.QDoubleSpinBox(); s.setRange(lo, hi); s.setDecimals(2); s.value = val; s.suffix = suffix; s.setToolTip(tip)
        return s

    # --- actions
    def onImport(self):
        files = qt.QFileDialog.getOpenFileNames(None, "SimVascular path files", "", "SimVascular path (*.pth)")
        for f in files:
            node = self.logic.importPth(f, self.unitsSpin.value)
            self.curvesSelector.setCheckState(node, qt.Qt.Checked)

    def onLoadSample(self):
        name = self.sampleCombo.currentText
        with slicer.util.tryWithErrorDisplay(f"Could not load sample {name}.", waitCursor=True):
            self.status.text = f"Loading SIMONE {name} (CT, design, centerlines)..."; slicer.app.processEvents()
            vol, model, curves, cfg = self.logic.loadSample(name, self.unitsSpin.value)
            self.ctSelector.setCurrentNode(vol); self.designSelector.setCurrentNode(model)
            for i in range(self.curvesSelector.nodeCount()):
                self.curvesSelector.setCheckState(self.curvesSelector.nodeFromIndex(i), qt.Qt.Unchecked)
            for c in curves: self.curvesSelector.setCheckState(c, qt.Qt.Checked)
            self.insertCheck.checked = cfg.get("insert_hu") is not None
            if cfg.get("insert_hu") is not None: self.insertHU.value = cfg["insert_hu"]
            self.excludeEdit.text = core.format_exclude(cfg.get("exclude_mm"))
            if cfg.get("_exclude_reason"): self.excludeEdit.setToolTip(f"SIMONE {name}: {cfg['_exclude_reason']}")
            self.status.text = f"Loaded SIMONE {name}. Click 'Measure print error'."

    def params(self):
        return dict(station_spacing_mm=self.spacing.value, max_search_radius_mm=self.search.value,
                    min_diameter_mm=self.minDiam.value, junction_margin_mm=self.junction.value,
                    end_trim_mm=self.trim.value, mismatch_mm=self.mismatch.value,
                    insert_hu=self.insertHU.value if self.insertCheck.checked else None,
                    exclude_mm=core.parse_exclude(self.excludeEdit.text))

    def onApply(self):
        ct, design, curves = self.ctSelector.currentNode(), self.designSelector.currentNode(), self.curvesSelector.checkedNodes()
        if not ct or not design or not curves:
            slicer.util.errorDisplay("Select the printed CT, the design model and at least one centerline."); return
        self.progress.visible = True; self.applyButton.enabled = False

        def progress(frac, msg):
            self.progress.value = int(100 * frac); self.status.text = msg; slicer.app.processEvents()
        try:
            with slicer.util.tryWithErrorDisplay("Measurement failed.", waitCursor=True):
                self.logic.run(ct, design, curves, self.params(), progress)
                self.showResults()
        finally:
            self.progress.visible = False; self.applyButton.enabled = True

    def showResults(self):
        lg = self.logic; reg = lg.result["registration"]
        all_ = [r for r in lg.summary if r["branch"] == "ALL"][0]
        self.status.text = (f"Registration: {reg['rotation_deg']:.1f}°, median wall distance {reg['median_wall_distance_mm']:.2f} mm, "
                            f"{100 * reg['wall_within_1mm']:.0f}% of wall within 1 mm.<br>"
                            f"Whole vessel (all branches): diameter error {all_['dD_mean_mm']:+.2f} ± {all_['dD_sd_mm']:.2f} mm "
                            f"({all_['dD_mean_pct']:+.1f}%), centerline error {all_['cl_dev_mean_mm']:.2f} mm, "
                            f"n = {all_['n_stations']} stations.")
        lg.showLayout(self.show3D.checked)
        lg.build3D(); lg.fit3D()
        lg.colour3D(COLOR_BY[self.colorBy.currentIndex][1], self.colorMax.value)
        self.connectPlots()
        self.vesselCombo.blockSignals(True); self.vesselCombo.clear()
        for name in lg.result["paths"]: self.vesselCombo.addItem(name)
        self.vesselCombo.blockSignals(False); self.saveButton.enabled = True
        self.onVesselChanged(0)

    def connectPlots(self):
        lm = slicer.app.layoutManager()
        for i in range(lm.plotViewCount):
            pw = lm.plotWidget(i)
            if pw in self._plotConnected: continue
            pw.mrmlPlotViewNode().SetInteractionMode(slicer.vtkMRMLPlotViewNode.InteractionModeSelectPoints)
            pw.plotView().connect("dataSelected(vtkStringArray*, vtkCollection*)", self.onPlotSelection)
            self._plotConnected.append(pw)

    def onPlotSelection(self, ids, selection):
        for k in range(ids.GetNumberOfValues()):
            sel = selection.GetItemAsObject(k)
            if sel is None or sel.GetNumberOfTuples() == 0: continue
            name = self.logic.seriesToPath.get(ids.GetValue(k))
            if name is None: continue
            self.vesselCombo.setCurrentText(name); self.stationSlider.value = int(sel.GetValue(0)); return

    def onVesselChanged(self, idx):
        name = self.vesselCombo.currentText
        if not name or self.logic.result is None: return
        n = len(self.logic.result["paths"][name]["s"])
        self.stationSlider.maximum = n - 1
        ok = np.nonzero(self.logic.result["paths"][name]["ok"])[0]
        self.stationSlider.value = int(ok[len(ok) // 2]) if len(ok) else 0
        self.onStationChanged(self.stationSlider.value)

    def onStationChanged(self, value):
        name = self.vesselCombo.currentText
        if not name or self.logic.result is None: return
        st = self.logic.inspect(name, int(value), follow3D=self.show3D.checked and self.follow3D.checked,
                                editable=self.editCheck.checked)
        self.logic.updatePlots(name, int(value))
        f = lambda x, fmt: "–" if x is None or not np.isfinite(x) else format(x, fmt)
        tag = (st['excluded'] or 'included') + (", at a junction" if st.get('at_junction') else "") + (", fluid-filled" if st.get('fluid_filled') else "") + (", edited" if st.get('edited') else "")
        self.info.text = (f"<b>{name}</b>, s = {st['s_mm']:.1f} mm &nbsp; ({tag})<br>"
                          f"Diameter: printed {f(st['D_print_mm'], '.2f')} mm, design {f(st['D_design_mm'], '.2f')} mm, "
                          f"error {f(st['dD_mm'], '+.2f')} mm ({f(st['dD_pct'], '+.1f')}%)<br>"
                          f"Centerline error: {f(st['cl_dev_mm'], '.2f')} mm")

    def onEditToggled(self, on):
        self.applyEditButton.enabled = on
        self.onStationChanged(self.stationSlider.value)

    def _afterEdit(self, msg):
        lg = self.logic; all_ = [r for r in lg.summary if r["branch"] == "ALL"][0]
        self.status.text = (f"{msg}<br>Whole vessel (all branches): diameter error {all_['dD_mean_mm']:+.2f} \u00b1 {all_['dD_sd_mm']:.2f} mm "
                            f"({all_['dD_mean_pct']:+.1f}%), centerline error {all_['cl_dev_mean_mm']:.2f} mm, n = {all_['n_stations']} stations.")
        lg.colour3D(COLOR_BY[self.colorBy.currentIndex][1], self.colorMax.value)
        self.onTableChoice(self.tableChoice.currentIndex)
        self.onStationChanged(self.stationSlider.value)

    def onApplyEdit(self):
        if self.logic.result is None or not self.vesselCombo.currentText: return
        with slicer.util.tryWithErrorDisplay("Could not apply the contour.", waitCursor=True):
            self._afterEdit(self.logic.applyContour(self.vesselCombo.currentText, int(self.stationSlider.value)))

    def onToggleStation(self):
        if self.logic.result is None or not self.vesselCombo.currentText: return
        with slicer.util.tryWithErrorDisplay("Could not change the station.", waitCursor=True):
            self._afterEdit(self.logic.toggleStation(self.vesselCombo.currentText, int(self.stationSlider.value)))

    def on3DOptions(self, *args):
        if self.logic.result is None: return
        self.logic.showLayout(self.show3D.checked)
        self.logic.colour3D(COLOR_BY[self.colorBy.currentIndex][1], self.colorMax.value)
        self.connectPlots(); self.onTableChoice(self.tableChoice.currentIndex)
        self.onStationChanged(self.stationSlider.value)

    def onTableChoice(self, idx):
        tv = slicer.mrmlScene.GetSingletonNode("VPETable", "vtkMRMLTableViewNode")
        node = self.logic.summaryTable if idx == 0 else self.logic.stationsTable
        if tv and node: tv.SetTableNodeID(node.GetID())

    def onSave(self):
        folder = qt.QFileDialog.getExistingDirectory(None, "Save tables to folder")
        if folder:
            lg = self.logic
            core.write_csv(lg.result["stations"], os.path.join(folder, "VPE_stations.csv"))
            core.write_csv(lg.summary, os.path.join(folder, "VPE_summary.csv"))
            ba = core.bland_altman_rows(lg.result["stations"])
            if ba: core.write_csv(ba, os.path.join(folder, "VPE_bland_altman.csv"))
            self.status.text = f"Saved VPE_stations.csv, VPE_summary.csv and VPE_bland_altman.csv to {folder}"

    def cleanup(self):
        pass


# ----------------------------------------------------------------------------------------------------------- logic
class VesselPrintErrorCheckerLogic(ScriptedLoadableModuleLogic):
    def __init__(self):
        ScriptedLoadableModuleLogic.__init__(self)
        self.result = None; self.summary = None; self.seriesToPath = {}
        self.stationsTable = self.summaryTable = None

    # --- inputs
    @staticmethod
    def sampleNames():
        f = os.path.join(SAMPLES_DIR, "samples.json")
        if not os.path.exists(f): return []
        import json
        return list(json.load(open(f)))

    def loadSample(self, name, mm_per_unit=10.0):
        """Load one SIMONE module from SIMONE_validation_files: CT (DICOM), design STL and .pth centerlines."""
        import glob, json
        from DICOMLib import DICOMUtils
        cfg = json.load(open(os.path.join(SAMPLES_DIR, "samples.json")))[name]
        with DICOMUtils.TemporaryDICOMDatabase() as db:
            DICOMUtils.importDicom(os.path.join(SAMPLES_DIR, cfg["ct_folder"]), db)
            uids = [se for pa in db.patients() for st in db.studiesForPatient(pa) for se in db.seriesForStudy(st)]
            ids = DICOMUtils.loadSeriesByUID(uids)
        vol = slicer.mrmlScene.GetNodeByID(ids[0]); vol.SetName(f"SIMONE {name} CT")
        model = slicer.util.loadModel(os.path.join(SAMPLES_DIR, cfg["design"])); model.SetName(f"SIMONE {name} design")
        curves = [self.importPth(f, cfg.get("mm_per_unit", mm_per_unit))
                  for f in sorted(glob.glob(os.path.join(SAMPLES_DIR, cfg["centerlines"], "*.pth")))]
        return vol, model, curves, cfg

    @staticmethod
    def importPth(path, mm_per_unit=10.0):
        """SimVascular .pth (LPS) -> open markups curve in RAS."""
        P = core.read_pth(path, mm_per_unit, lps_to_ras=True)
        node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsCurveNode", os.path.splitext(os.path.basename(path))[0])
        node.SetCurveTypeToLinear()
        slicer.util.updateMarkupsControlPointsFromArray(node, P)
        dn = node.GetDisplayNode(); dn.SetPointLabelsVisibility(False); dn.SetPropertiesLabelVisibility(False)
        dn.SetGlyphScale(0.5); node.SetLocked(True)
        return node

    @staticmethod
    def _world_matrix(node):
        m = vtk.vtkMatrix4x4(); M = np.eye(4)
        if node.GetParentTransformNode():
            if not node.GetParentTransformNode().IsTransformToWorldLinear():
                raise ValueError(f"{node.GetName()} is under a non-linear transform; harden it first.")
            node.GetParentTransformNode().GetMatrixTransformToWorld(m); M = slicer.util.arrayFromVTKMatrix(m)
        return M

    def _ct(self, volNode, insert_hu):
        m = vtk.vtkMatrix4x4(); volNode.GetIJKToRASMatrix(m)
        ijk_to_world = self._world_matrix(volNode) @ slicer.util.arrayFromVTKMatrix(m)
        return core.CTVolume(slicer.util.arrayFromVolume(volNode), np.linalg.inv(ijk_to_world), insert_hu, 1.0, min(volNode.GetSpacing()))

    def _design(self, modelNode):
        poly = vtk.vtkPolyData(); poly.DeepCopy(modelNode.GetPolyData())
        if modelNode.GetParentTransformNode():
            t = vtk.vtkGeneralTransform(); modelNode.GetParentTransformNode().GetTransformToWorld(t)
            f = vtk.vtkTransformPolyDataFilter(); f.SetTransform(t); f.SetInputData(poly); f.Update(); poly = f.GetOutput()
        tri = vtk.vtkTriangleFilter(); tri.SetInputData(poly); tri.Update()
        return tri.GetOutput()

    # --- run
    def run(self, volNode, modelNode, curveNodes, params, progress=None):
        ct = self._ct(volNode, params.get("insert_hu")); self.ct = ct; self.step = params.get("station_spacing_mm", 0.5)
        curves = {}
        for c in curveNodes:
            P = slicer.util.arrayFromMarkupsCurvePoints(c, world=True)
            if len(P) < 2: continue
            curves[c.GetName()] = P
        self.volumeNode = volNode
        self.designPoly = self._design(modelNode); self.inputDesignNode = modelNode
        self.result = core.run(ct, curves, self.designPoly, params, progress)
        self.refresh(rebuild3D=False)

    def refresh(self, rebuild3D=True):
        """Recompute the summary, tables and plots (and the 3D surface) from the per-station results."""
        self.summary = core.summarise(self.result["stations"], self.step)
        self.stationsTable = self._table("VPE stations", self.result["stations"])
        self.summaryTable = self._table("VPE summary", self.summary)
        self._plots()
        if rebuild3D: self.build3D()

    # --- manual corrections
    def applyContour(self, name, i):
        """Replace station i's printed contour with the (edited) curve shown in the cross-section view."""
        node = slicer.mrmlScene.GetFirstNodeByName("VPE printed contour")
        pts = slicer.util.arrayFromMarkupsCurvePoints(node, world=True) if node else np.zeros((0, 3))
        if len(pts) < 8: return "No contour to apply at this station."
        res = self.result; p = res["paths"][name]; row = [s for s in res["stations"] if s["branch"] == name][i]
        c0, n = p["C"][i], p["T"][i]; U, _ = core.frames(p["T"][max(0, i - 1):i + 1]); u = U[-1]; v = np.cross(n, u)
        xy = np.c_[(pts - c0) @ u, (pts - c0) @ v]; A, cen = core.polygon(xy)
        if A <= 0: return "The contour is not a valid outline."
        c = c0 + cen[0] * u + cen[1] * v; D = 2 * np.sqrt(A / np.pi)
        ang = np.arctan2(xy[:, 1] - cen[1], xy[:, 0] - cen[0]) % (2 * np.pi); rad = np.hypot(xy[:, 0] - cen[0], xy[:, 1] - cen[1])
        o = np.argsort(ang); rr = np.interp(core.ANG, np.r_[ang[o] - 2 * np.pi, ang[o], ang[o] + 2 * np.pi], np.tile(rad[o], 3))
        p["rings"][i] = c + (rr * np.cos(core.ANG))[:, None] * u + (rr * np.sin(core.ANG))[:, None] * v
        p["C"][i] = c; p["D"][i] = D; p["valid"][i] = True
        row.update(D_print_mm=float(D), x=float(c[0]), y=float(c[1]), z=float(c[2]), edited=1)
        Dd = row["D_design_mm"]
        if not np.isfinite(Dd):
            self.refresh(); return "Contour saved, but there is no design cross-section here to compare with, so the station stays excluded."
        ref = p["Co"][np.isfinite(p["Co"]).all(1) & p["ok"]]
        cl = float(core.cKDTree(core.densify(ref)).query(c @ res["R"].T + res["t"])[0]) if len(ref) >= 2 else float("nan")
        row.update(dD_mm=float(D - Dd), dD_pct=float(100 * (D - Dd) / Dd), cl_dev_mm=cl, excluded="")
        p["ok"][i] = True
        self.refresh(); return f"Contour applied: printed diameter {D:.2f} mm, error {D - Dd:+.2f} mm."

    def toggleStation(self, name, i):
        """Exclude an included station, or bring back one that has both diameters."""
        p = self.result["paths"][name]; row = [s for s in self.result["stations"] if s["branch"] == name][i]
        if row["excluded"] == "":
            row["excluded"] = "user excluded"; p["ok"][i] = False; msg = "Station excluded."
        elif np.isfinite(row["D_print_mm"]) and np.isfinite(row["D_design_mm"]) and np.isfinite(row["dD_mm"]):
            msg = f"Station included (it was: {row['excluded']})."; row["excluded"] = ""; row["edited"] = 1; p["ok"][i] = True
        else:
            return "This station has no printed or no design cross-section, so it cannot be included."
        self.refresh(); return msg

    # --- outputs
    @staticmethod
    def _node(cls, name, hidden=False):
        n = slicer.mrmlScene.GetFirstNodeByName(name)
        if n is None or not n.IsA(cls):
            n = slicer.mrmlScene.AddNewNodeByClass(cls, name)
        n.SetHideFromEditors(hidden)
        return n

    def _table(self, name, rows):
        node = self._node("vtkMRMLTableNode", name); node.RemoveAllColumns()
        keys = list(rows[0].keys())
        for k in keys:
            vals = [r[k] for r in rows]
            if isinstance(vals[0], str):
                arr = vtk.vtkStringArray()
                for v in vals: arr.InsertNextValue(v)
            else:
                arr = vtk.vtkDoubleArray()
                for v in vals: arr.InsertNextValue(np.nan if v is None else float(v))
            arr.SetName(k); node.AddColumn(arr)
        node.Modified(); return node

    def _plots(self):
        """One table per vessel (all stations; excluded ones blank so they leave gaps) and one set of series per
        vessel. The two charts only ever show the vessel being inspected (see updatePlots)."""
        for n in getattr(self, "_plotNodes", []):
            if slicer.mrmlScene.IsNodePresent(n): slicer.mrmlScene.RemoveNode(n)
        self._plotNodes, self.seriesToPath, self.vesselSeries, self._plotted = [], {}, {}, None
        chartD = self._node("vtkMRMLPlotChartNode", "VPE diameter error")
        chartC = self._node("vtkMRMLPlotChartNode", "VPE centerline error")
        for ch in (chartD, chartC):
            ch.RemoveAllPlotSeriesNodeIDs(); ch.SetXAxisTitle("distance along centerline (mm)")
        chartD.SetYAxisTitle("print \u2212 design (mm)"); chartD.SetLegendVisibility(False)
        chartC.SetYAxisTitle("print to design (mm)"); chartC.SetLegendVisibility(False)

        def series(name, table, x, y, color, dashed=False, marker=True, width=2):
            sn = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotSeriesNode", name); sn.SetHideFromEditors(True)
            sn.SetAndObserveTableNodeID(table.GetID()); sn.SetXColumnName(x); sn.SetYColumnName(y)
            sn.SetPlotType(slicer.vtkMRMLPlotSeriesNode.PlotTypeScatter)
            sn.SetMarkerStyle(slicer.vtkMRMLPlotSeriesNode.MarkerStyleCircle if marker else slicer.vtkMRMLPlotSeriesNode.MarkerStyleNone)
            sn.SetMarkerSize(3)
            sn.SetLineStyle(slicer.vtkMRMLPlotSeriesNode.LineStyleDash if dashed else slicer.vtkMRMLPlotSeriesNode.LineStyleSolid)
            sn.SetLineWidth(width); sn.SetColor(*color); self._plotNodes.append(sn); return sn

        for name, p in self.result["paths"].items():
            rows = [st for st in self.result["stations"] if st["branch"] == name]; ok = p["ok"]
            masked = lambda key: [r[key] if ok[i] else np.nan for i, r in enumerate(rows)]
            t = self._table(f"VPE plot {name}", [dict(s_mm=r["s_mm"], dD_mm=a, cl_mm=b) for r, a, b in
                                                  zip(rows, masked("dD_mm"), masked("cl_dev_mm"))])
            t.SetHideFromEditors(True)
            sD = series("diameter error", t, "s_mm", "dD_mm", (0.16, 0.47, 0.84))
            sC = series("centerline error", t, "s_mm", "cl_mm", (0.80, 0.45, 0.12))
            for sn in (sD, sC): self.seriesToPath[sn.GetID()] = name
            self.vesselSeries[name] = (sD, sC, rows, ok)
        # zero line and station marker, shared by both charts
        self.guideTable = self._table("VPE plot guides", [dict(x=0.0, zero=0.0, yD=0.0, yC=0.0), dict(x=1.0, zero=0.0, yD=0.0, yC=0.0)])
        self.markerTable = self._table("VPE plot marker", [dict(x=0.0, yD=0.0, yC=0.0), dict(x=0.0, yD=0.0, yC=0.0)])
        for t in (self.guideTable, self.markerTable): t.SetHideFromEditors(True)
        self.zeroSeries = series("zero", self.guideTable, "x", "zero", (0.55, 0.55, 0.55), marker=False, width=1)
        self.markD = series("inspected station", self.markerTable, "x", "yD", (0.35, 0.35, 0.35), marker=False, width=1)
        self.markC = series("inspected station", self.markerTable, "x", "yC", (0.35, 0.35, 0.35), marker=False, width=1)
        self.chartD, self.chartC = chartD, chartC

    def updatePlots(self, name, i):
        """Show only branch `name` in both charts, with a vertical line at station i."""
        if name not in getattr(self, "vesselSeries", {}): return
        sD, sC, rows, ok = self.vesselSeries[name]
        s_all = np.array([r["s_mm"] for r in rows])
        val = lambda k: np.array([r[k] if ok[j] else np.nan for j, r in enumerate(rows)], float)
        dD, cl = val("dD_mm"), val("cl_dev_mm")
        fin = lambda a: a[np.isfinite(a)]
        lo = min(0.0, *fin(dD)) if fin(dD).size else -1.0; hi = max(0.0, *fin(dD)) if fin(dD).size else 1.0
        top = max([1e-3] + list(fin(cl)))
        if self._plotted != name:
            for ch in (self.chartD, self.chartC): ch.RemoveAllPlotSeriesNodeIDs()
            for sn in (self.zeroSeries, sD, self.markD): self.chartD.AddAndObservePlotSeriesNodeID(sn.GetID())
            for sn in (self.zeroSeries, sC, self.markC): self.chartC.AddAndObservePlotSeriesNodeID(sn.GetID())
            n_ok = int(ok.sum()); m = fin(dD)
            self.chartD.SetTitle(f"{name}: diameter error" + (f"  (mean {m.mean():+.2f} mm, n = {n_ok})" if m.size else "  (no included stations)"))
            c = fin(cl)
            self.chartC.SetTitle(f"{name}: centerline error" + (f"  (mean {c.mean():.2f} mm)" if c.size else
                                                                 "  (no included stations)"))
            g = self.guideTable.GetTable(); g.SetValue(0, 0, float(s_all[0])); g.SetValue(1, 0, float(s_all[-1])); g.Modified()
            self._plotted = name
        t = self.markerTable.GetTable(); x = float(s_all[i])
        for r, (yd, yc) in enumerate([(lo - 0.1 * (hi - lo + 1e-6), 0.0), (hi + 0.1 * (hi - lo + 1e-6), top * 1.1)]):
            t.SetValue(r, 0, x); t.SetValue(r, 1, float(yd)); t.SetValue(r, 2, float(yc))
        t.Modified(); self.markerTable.Modified()

    def showLayout(self, with3D=False):
        lm = slicer.app.layoutManager(); ln = lm.layoutLogic().GetLayoutNode()
        if not ln.IsLayoutDescription(LAYOUT_ID): ln.AddLayoutDescription(LAYOUT_ID, LAYOUT_XML)
        if not ln.IsLayoutDescription(LAYOUT_3D_ID): ln.AddLayoutDescription(LAYOUT_3D_ID, LAYOUT_3D_XML)
        target = LAYOUT_3D_ID if with3D else LAYOUT_ID
        if lm.layout != target:
            lm.setLayout(target); slicer.app.processEvents()
            if with3D: self.fit3D()
        for tag, chart in [("VPEPlotDiameter", self.chartD), ("VPEPlotCenterline", self.chartC)]:
            v = slicer.mrmlScene.GetSingletonNode(tag, "vtkMRMLPlotViewNode")
            if v: v.SetPlotChartNodeID(chart.GetID())
        tv = slicer.mrmlScene.GetSingletonNode("VPETable", "vtkMRMLTableViewNode")
        if tv: tv.SetTableNodeID(self.summaryTable.GetID())
        slicer.util.setSliceViewerLayers(background=self.volumeNode, fit=False)
        vd = self.volumeNode.GetDisplayNode()
        if vd and not self.volumeNode.GetAttribute("VPE.WindowSet"):          # once per volume: air -> resin contrast
            vd.AutoWindowLevelOff(); vd.SetWindowLevel(1800, -100); self.volumeNode.SetAttribute("VPE.WindowSet", "1")

    def inspect(self, name, i, follow3D=False, editable=False):
        """Reorient the Red slice view onto station i of vessel `name` and draw both contours."""
        p = self.result["paths"][name]
        row = [s for s in self.result["stations"] if s["branch"] == name][i]
        c, n = p["C"][i], p["T"][i]
        U, _ = core.frames(p["T"][max(0, i - 1):i + 1]); u = U[-1]
        sn = slicer.mrmlScene.GetNodeByID("vtkMRMLSliceNodeRed")
        pts = [q for q in (p["rings"][i], p["design_loops_in_ct"][i]) if q is not None and np.isfinite(q).all()]
        if pts:                                  # frame both contours: centre between them, wide enough for both
            allp = np.vstack(pts); c = (allp.min(0) + allp.max(0)) / 2
            c = c - np.dot(c - p["C"][i], n) * n                        # stay in the station's plane
            r = float(np.max(np.linalg.norm(allp - c, axis=1)))
        else:
            r = 5.0
        if sn is not None:
            sn.SetSliceToRASByNTP(*n, *u, *c, 0)
            sn.SetXYZOrigin(0, 0, 0)             # undo any pan left over from the user, so the view is centred on c
            fov = 2 * r + 6.0; sn.SetFieldOfView(fov, fov, sn.GetFieldOfView()[2])
            sn.UpdateMatrices()
        if follow3D and slicer.app.layoutManager().threeDViewCount:
            viewNode = slicer.app.layoutManager().threeDWidget(0).mrmlViewNode()
            cam = slicer.modules.cameras.logic().GetViewActiveCameraNode(viewNode)
            if cam:
                fp, pos = np.array(cam.GetFocalPoint()), np.array(cam.GetPosition())
                cam.SetFocalPoint(*c); cam.SetPosition(*(c + pos - fp))
        self._contour("VPE printed contour", p["rings"][i] if np.isfinite(p["rings"][i]).all() else None, (0.13, 0.83, 0.93), editable)
        self._contour("VPE design contour", p["design_loops_in_ct"][i], (0.96, 0.62, 0.04))
        return row

    def build3D(self):
        """Printed lumen as one continuous surface segmented from the CT (junctions included), carrying the errors
        of the nearest included station, plus the design surface mapped into CT space."""
        res = self.result
        poly = core.lumen_surface(self.ct, res)
        m = self._node("vtkMRMLModelNode", "VPE printed lumen"); m.SetAndObservePolyData(poly)
        if m.GetDisplayNode() is None: m.CreateDefaultDisplayNodes()
        m.GetDisplayNode().SetVisibility2D(False); m.GetDisplayNode().SetBackfaceCulling(False)
        R, t = res["R"], res["t"]; M = np.eye(4); M[:3, :3] = R.T; M[:3, 3] = -R.T @ t   # design -> CT
        tr = vtk.vtkTransform(); tr.SetMatrix(slicer.util.vtkMatrixFromArray(M))
        f = vtk.vtkTransformPolyDataFilter(); f.SetTransform(tr); f.SetInputData(self.designPoly); f.Update()
        d = self._node("vtkMRMLModelNode", "VPE design (registered to CT)"); d.SetAndObservePolyData(f.GetOutput())
        if d.GetDisplayNode() is None: d.CreateDefaultDisplayNodes()
        dd = d.GetDisplayNode(); dd.SetColor(0.75, 0.75, 0.75); dd.SetOpacity(0.25); dd.SetVisibility2D(False)
        self.lumenModel, self.designModel = m, d
        # the input design sits in its own (design) frame, away from the CT; hide it in 3D so the camera frames the
        # printed lumen and the registered copy of the design
        idn = self.inputDesignNode.GetDisplayNode() if getattr(self, "inputDesignNode", None) else None
        if idn: idn.SetVisibility3D(False)

    def fit3D(self):
        """Point the 3D camera at the printed lumen (not the whole CT box) and hide the box and axis labels."""
        lm = slicer.app.layoutManager(); m = getattr(self, "lumenModel", None)
        if not lm.threeDViewCount or m is None or m.GetPolyData() is None or m.GetPolyData().GetNumberOfPoints() == 0: return
        w = lm.threeDWidget(0); vn = w.mrmlViewNode(); vn.SetBoxVisible(False); vn.SetAxisLabelsVisible(False)
        ren = w.threeDView().renderWindow().GetRenderers().GetFirstRenderer()
        ren.ResetCamera(m.GetPolyData().GetBounds()); w.threeDView().forceRender()

    def _colorTable(self, kind):
        """Blue-white-red for diameter error (signed), cream-to-dark-red for centerline error; grey = no value."""
        name = "VPE diverging" if kind == "div" else "VPE sequential"
        node = slicer.mrmlScene.GetFirstNodeByName(name)
        if node is None:
            node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLColorTableNode", name)
            node.SetTypeToUser(); node.SetNumberOfColors(256); node.SetHideFromEditors(True)
            stops = ([(0, (0.13, 0.30, 0.67)), (0.5, (0.97, 0.97, 0.97)), (1, (0.70, 0.09, 0.17))] if kind == "div" else
                     [(0, (1.0, 0.97, 0.88)), (0.4, (0.99, 0.70, 0.30)), (1, (0.50, 0.0, 0.15))])
            for k in range(256):
                x = k / 255.0
                for (x0, c0), (x1, c1) in zip(stops[:-1], stops[1:]):
                    if x0 <= x <= x1:
                        w = (x - x0) / (x1 - x0); node.SetColor(k, *[a + w * (b - a) for a, b in zip(c0, c1)], 1.0); break
            node.GetLookupTable().SetNanColor(0.6, 0.6, 0.6, 1.0)
        return node

    def colour3D(self, key, limit):
        m = getattr(self, "lumenModel", None)
        if m is None: return
        dn = m.GetDisplayNode()
        if key is None:
            dn.SetScalarVisibility(False); dn.SetColor(0.13, 0.83, 0.93)
        else:
            signed = key == "dD_mm"; ct = self._colorTable("div" if signed else "seq")
            dn.SetActiveScalarName(key); dn.SetAndObserveColorNodeID(ct.GetID())
            dn.SetScalarRangeFlag(dn.UseManualScalarRange); dn.SetScalarRange(-limit if signed else 0.0, limit)
            dn.SetScalarVisibility(True)
        try:                                     # colour legend in the 3D view, rebuilt so it follows the colouring
            cl = slicer.modules.colors.logic()
            old = cl.GetColorLegendDisplayNode(m)
            if old: slicer.mrmlScene.RemoveNode(old)
            if key is not None:
                legend = cl.AddDefaultColorLegendDisplayNode(m)
                legend.SetTitleText(LEGEND_TITLE[key])
                legend.SetSize(0.12, 0.45); legend.SetPosition(0.98, 0.5)
                views = slicer.util.getNodesByClass("vtkMRMLViewNode")
                if views: legend.SetViewNodeIDs([v.GetID() for v in views])
        except Exception:
            pass

    def _contour(self, name, pts, color, editable=False):
        node = self._node("vtkMRMLMarkupsClosedCurveNode", name, hidden=True)
        if node.GetDisplayNode() is None: node.CreateDefaultDisplayNodes()
        dn = node.GetDisplayNode()
        dn.SetPointLabelsVisibility(False); dn.SetPropertiesLabelVisibility(False)
        dn.SetUseGlyphScale(False)
        if editable:                              # a few draggable points on a smooth closed curve
            node.SetCurveTypeToKochanekSpline(); node.SetLocked(False); dn.SetGlyphSize(0.35)
            if pts is not None and len(pts) >= 3: pts = np.asarray(pts)[np.linspace(0, len(pts), 16, endpoint=False).astype(int)]
        else:
            node.SetCurveTypeToLinear(); node.SetLocked(True); dn.SetGlyphSize(0.05)   # hide points, keep the line
        dn.SetCurveLineSizeMode(dn.UseLineDiameter); dn.SetLineDiameter(0.12)
        dn.SetColor(*color); dn.SetSelectedColor(*color)
        dn.SetSliceProjection(True); dn.SetSliceProjectionUseFiducialColor(True); dn.SetSliceProjectionOpacity(1.0)
        if pts is None or len(pts) < 3:
            node.RemoveAllControlPoints(); return
        pts = np.asarray(pts)
        if len(pts) > 200: pts = pts[np.linspace(0, len(pts) - 1, 200).astype(int)]
        slicer.util.updateMarkupsControlPointsFromArray(node, pts)
