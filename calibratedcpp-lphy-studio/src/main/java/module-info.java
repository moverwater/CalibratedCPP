import prior.lphystudio.viewer.CalibrationArrayViewer;

module calibratedcpp.lphy.studio {
    requires transitive lphystudio;
    requires calibratedcpp.lphy;

    // Viewer SPI
    uses lphystudio.app.graphicalmodelpanel.viewer.Viewer;
    // declare what service interface the provider intends to use
    provides lphystudio.app.graphicalmodelpanel.viewer.Viewer with CalibrationArrayViewer;

    // Note: to adapt with the system not using Java module but using class path,
    // they need to be declared inside META-INF/services/lphystudio.app.graphicalmodelpanel.viewer.Viewer as well.
}