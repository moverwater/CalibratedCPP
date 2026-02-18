import prior.lphystudio.viewer.CalibrationArrayViewer;

module calibratedcpp.lphy.studio {
    requires transitive lphystudio;
    requires calibratedcpp.lphy;

    uses lphystudio.app.graphicalmodelpanel.viewer.Viewer;
    provides lphystudio.app.graphicalmodelpanel.viewer.Viewer with CalibrationArrayViewer;

}