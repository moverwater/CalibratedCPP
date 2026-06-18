open module calibratedcpp.beast {
    requires beast.pkgmgmt;
    requires beast.base;
    requires beast.fx;
    requires hipparchus.core;

    // JavaFX modules used by CalibratedCPPInputEditor (beast.fx does not re-export these)
    requires javafx.graphics;
    requires javafx.controls;
    requires javafx.web;
    requires jdk.jsobject;
    requires org.apache.commons.statistics.distribution;

    exports calibratedcpp;
    exports calibration;
    exports calibrationprior;
    exports calibratedcpp.beauti;

    // CalibratedCoalescentPointProcess (abstract), CalibrationNode/Forest (no default ctor),
    // and CalibratedCPPInputEditor (not a BEASTInterface subtype) are registered via version.xml only.
    provides beast.base.core.BEASTInterface with
        calibratedcpp.CalibratedBirthDeathModel,
        calibratedcpp.CalibratedBirthDeathSkylineModel,
        calibratedcpp.CalibratedAgeDependentBirthDeathModel,
        calibratedcpp.SkylineParameter,
        calibrationprior.CalibrationPrior,
        calibrationprior.CalibrationCladePrior,
        calibration.CalibrationClade;
}
