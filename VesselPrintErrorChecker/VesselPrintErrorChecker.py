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

        # --- parameters
        pbox = ctk.ctkCollapsibleButton(); pbox.text = "Parameters"; pbox.collapsed = True; self.layout.addWidget(pbox)
        pform = qt.QFormLayout(pbox)
        d = core.DEFAULTS
        self.spacing = self._spin(0.1, 5, d["station_spacing_mm"], " mm", "Distance between cross-sections")
        self.search = self._spin(2, 50, d["max_search_radius_mm"], " mm", "Largest distance from the centerline searched for the wall (first pass)")
        self.minDiam = self._spin(0, 10, d["min_diameter_mm"], " mm", "Design diameters below this are excluded as below CT resolution (0 = keep all)")
        self.junction = self._spin(0, 20, d["junction_margin_mm"], " mm", "Extra exclusion either side of a vessel junction")
        self.trim = self._spin(0, 20, d["end_trim_mm"], " mm", "Ignore this much at each end of a centerline")
        self.mismatch = self._spin(0.5, 20, d["mismatch_mm"], " mm", "Print and design centres further apart than this (after the per-vessel fit) are not compared")
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
        self.vesselCombo = qt.QComboBox(); iform.addRow("Vessel:", self.vesselCombo)
        self.stationSlider = ctk.ctkSliderWidget(); self.stationSlider.decimals = 0; self.stationSlider.singleStep = 1
        self.stationSlider.minimum = 0; self.stationSlider.maximum = 0
        iform.addRow("Station:", self.stationSlider)
        self.info = qt.QLabel("Run a measurement, then pick a station here or click a point in either plot.")
        self.info.wordWrap = True; iform.addRow(self.info)
        legend = qt.QLabel("<span style='color:#22d3ee'>■</span> printed lumen &nbsp; "
                           "<span style='color:#f59e0b'>■</span> design lumen (mapped into the CT)")
        iform.addRow(legend)

        rrow = qt.QHBoxLayout()
        self.tableChoice = qt.QComboBox(); self.tableChoice.addItems(["Summary table", "Per-station table"])
        self.saveButton = qt.QPushButton("Save tables as CSV..."); self.saveButton.enabled = False
        rrow.addWidget(self.tableChoice); rrow.addWidget(self.saveButton); self.layout.addLayout(rrow)
        self.layout.addStretch(1)

        self.importButton.connect("clicked()", self.onImport)
        self.sampleButton.connect("clicked()", self.onLoadSample)
        self.applyButton.connect("clicked()", self.onApply)
        self.vesselCombo.connect("currentIndexChanged(int)", self.onVesselChanged)
        self.stationSlider.connect("valueChanged(double)", self.onStationChanged)
        self.saveButton.connect("clicked()", self.onSave)
        self.tableChoice.connect("currentIndexChanged(int)", self.onTableChoice)
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
            self.status.text = f"Loaded SIMONE {name}. Click 'Measure print error'."

    def params(self):
        return dict(station_spacing_mm=self.spacing.value, max_search_radius_mm=self.search.value,
                    min_diameter_mm=self.minDiam.value, junction_margin_mm=self.junction.value,
                    end_trim_mm=self.trim.value, mismatch_mm=self.mismatch.value,
                    insert_hu=self.insertHU.value if self.insertCheck.checked else None)

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
        all_ = [r for r in lg.summary if r["vessel"] == "ALL"][0]
        self.status.text = (f"Registration: {reg['rotation_deg']:.1f}°, median wall distance {reg['median_wall_distance_mm']:.2f} mm, "
                            f"{100 * reg['wall_within_1mm']:.0f}% of wall within 1 mm.<br>"
                            f"All vessels: diameter error {all_['dD_mean_mm']:+.2f} ± {all_['dD_sd_mm']:.2f} mm "
                            f"({all_['dD_mean_pct']:+.1f}%), centerline error {all_['cl_module_mean_mm']:.2f} mm (module fit), "
                            f"{all_['cl_vessel_mean_mm']:.2f} mm (per-vessel fit), n = {all_['n_stations']} stations.")
        lg.showLayout()
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
        st = self.logic.inspect(name, int(value))
        f = lambda x, fmt: "–" if x is None or not np.isfinite(x) else format(x, fmt)
        self.info.text = (f"<b>{name}</b>, s = {st['s_mm']:.1f} mm &nbsp; ({st['excluded'] or 'included'})<br>"
                          f"Diameter: printed {f(st['D_print_mm'], '.2f')} mm, design {f(st['D_design_mm'], '.2f')} mm, "
                          f"error {f(st['dD_mm'], '+.2f')} mm ({f(st['dD_pct'], '+.1f')}%)<br>"
                          f"Centerline error: {f(st['cl_dev_module_mm'], '.2f')} mm (module fit), "
                          f"{f(st['cl_dev_vessel_mm'], '.2f')} mm (per-vessel fit)")

    def onTableChoice(self, idx):
        tv = slicer.mrmlScene.GetSingletonNode("VPETable", "vtkMRMLTableViewNode")
        node = self.logic.summaryTable if idx == 0 else self.logic.stationsTable
        if tv and node: tv.SetTableNodeID(node.GetID())

    def onSave(self):
        folder = qt.QFileDialog.getExistingDirectory(None, "Save tables to folder")
        if folder:
            for node in (self.logic.stationsTable, self.logic.summaryTable):
                slicer.util.saveNode(node, os.path.join(folder, node.GetName().replace(" ", "_") + ".csv"))
            self.status.text = f"Saved tables to {folder}"

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
        ct = self._ct(volNode, params.get("insert_hu"))
        curves = {}
        for c in curveNodes:
            P = slicer.util.arrayFromMarkupsCurvePoints(c, world=True)
            if len(P) < 2: continue
            curves[c.GetName()] = P
        self.volumeNode = volNode
        self.result = core.run(ct, curves, self._design(modelNode), params, progress)
        self.summary = core.summarise(self.result["stations"], params.get("station_spacing_mm", 0.5))
        self.stationsTable = self._table("VPE stations", self.result["stations"])
        self.summaryTable = self._table("VPE summary", self.summary)
        self._plots()

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
        chartD = self._node("vtkMRMLPlotChartNode", "VPE diameter error")
        chartC = self._node("vtkMRMLPlotChartNode", "VPE centerline error")
        for ch in (chartD, chartC): ch.RemoveAllPlotSeriesNodeIDs()
        chartD.SetTitle("Diameter error, print − design (excluded stations not drawn)")
        chartD.SetXAxisTitle("distance along centerline (mm)"); chartD.SetYAxisTitle("ΔD (mm)")
        chartC.SetTitle("Centerline error: solid = per-vessel fit, dashed = whole-module fit")
        chartC.SetXAxisTitle("distance along centerline (mm)"); chartC.SetYAxisTitle("deviation (mm)")
        self.seriesToPath = {}
        for k, (name, p) in enumerate(self.result["paths"].items()):
            rows = [s for s in self.result["stations"] if s["path"] == name]
            ok = p["ok"]; col = COLORS[k % len(COLORS)]
            masked = lambda key: [r[key] if ok[i] else np.nan for i, r in enumerate(rows)]
            t = self._table(f"VPE plot {name}", [dict(s_mm=r["s_mm"], dD_mm=a, cl_module_mm=b, cl_vessel_mm=c)
                                                  for r, a, b, c in zip(rows, masked("dD_mm"), masked("cl_dev_module_mm"), masked("cl_dev_vessel_mm"))])
            t.SetHideFromEditors(True)
            for chart, ycol, label, dashed in [(chartD, "dD_mm", name, False), (chartC, "cl_vessel_mm", f"{name} (vessel)", False),
                                               (chartC, "cl_module_mm", f"{name} (module)", True)]:
                s = self._node("vtkMRMLPlotSeriesNode", f"VPE {chart.GetName()[4:]} {label}", hidden=True)
                s.SetName(label); s.SetAndObserveTableNodeID(t.GetID())
                s.SetXColumnName("s_mm"); s.SetYColumnName(ycol)
                s.SetPlotType(slicer.vtkMRMLPlotSeriesNode.PlotTypeScatter)
                s.SetMarkerStyle(slicer.vtkMRMLPlotSeriesNode.MarkerStyleNone if dashed else slicer.vtkMRMLPlotSeriesNode.MarkerStyleCircle)
                s.SetMarkerSize(4)
                s.SetLineStyle(slicer.vtkMRMLPlotSeriesNode.LineStyleDash if dashed else slicer.vtkMRMLPlotSeriesNode.LineStyleSolid)
                s.SetLineWidth(2); s.SetColor(*col)
                chart.AddAndObservePlotSeriesNodeID(s.GetID()); self.seriesToPath[s.GetID()] = name
        self.chartD, self.chartC = chartD, chartC

    def showLayout(self):
        lm = slicer.app.layoutManager(); ln = lm.layoutLogic().GetLayoutNode()
        if not ln.IsLayoutDescription(LAYOUT_ID): ln.AddLayoutDescription(LAYOUT_ID, LAYOUT_XML)
        lm.setLayout(LAYOUT_ID); slicer.app.processEvents()
        for tag, chart in [("VPEPlotDiameter", self.chartD), ("VPEPlotCenterline", self.chartC)]:
            v = slicer.mrmlScene.GetSingletonNode(tag, "vtkMRMLPlotViewNode")
            if v: v.SetPlotChartNodeID(chart.GetID())
        tv = slicer.mrmlScene.GetSingletonNode("VPETable", "vtkMRMLTableViewNode")
        if tv: tv.SetTableNodeID(self.summaryTable.GetID())
        slicer.util.setSliceViewerLayers(background=self.volumeNode, fit=False)
        vd = self.volumeNode.GetDisplayNode()
        if vd and not self.volumeNode.GetAttribute("VPE.WindowSet"):          # once per volume: air -> resin contrast
            vd.AutoWindowLevelOff(); vd.SetWindowLevel(1800, -100); self.volumeNode.SetAttribute("VPE.WindowSet", "1")

    def inspect(self, name, i):
        """Reorient the Red slice view onto station i of vessel `name` and draw both contours."""
        p = self.result["paths"][name]
        row = [s for s in self.result["stations"] if s["path"] == name][i]
        c, n = p["C"][i], p["T"][i]
        U, _ = core.frames(p["T"][max(0, i - 1):i + 1]); u = U[-1]
        sn = slicer.mrmlScene.GetNodeByID("vtkMRMLSliceNodeRed")
        if sn is not None:
            sn.SetSliceToRASByNTP(*n, *u, *c, 0)
            ring = p["rings"][i]
            r = np.nanmax(np.linalg.norm(ring - c, axis=1)) if np.isfinite(ring).any() else 5.0
            fov = 2 * (r + 4.0); sn.SetFieldOfView(fov, fov, sn.GetFieldOfView()[2])
        self._contour("VPE printed contour", p["rings"][i] if np.isfinite(p["rings"][i]).all() else None, (0.13, 0.83, 0.93))
        self._contour("VPE design contour", p["design_loops_in_ct"][i], (0.96, 0.62, 0.04))
        return row

    def _contour(self, name, pts, color):
        node = self._node("vtkMRMLMarkupsClosedCurveNode", name, hidden=True)
        if node.GetDisplayNode() is None: node.CreateDefaultDisplayNodes()
        node.SetCurveTypeToLinear(); node.SetLocked(True)
        dn = node.GetDisplayNode()
        dn.SetPointLabelsVisibility(False); dn.SetPropertiesLabelVisibility(False)
        dn.SetUseGlyphScale(False); dn.SetGlyphSize(0.05)                     # hide control points, keep the line
        dn.SetCurveLineSizeMode(dn.UseLineDiameter); dn.SetLineDiameter(0.12)
        dn.SetColor(*color); dn.SetSelectedColor(*color)
        dn.SetSliceProjection(True); dn.SetSliceProjectionUseFiducialColor(True); dn.SetSliceProjectionOpacity(1.0)
        if pts is None or len(pts) < 3:
            node.RemoveAllControlPoints(); return
        pts = np.asarray(pts)
        if len(pts) > 200: pts = pts[np.linspace(0, len(pts) - 1, 200).astype(int)]
        slicer.util.updateMarkupsControlPointsFromArray(node, pts)
