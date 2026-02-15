module calibratedcpp.lphy {
    requires transitive lphy.base;
    requires jdk.jfr;

    exports calibratedcpp.lphy.tree;
    exports calibratedcpp.lphy.prior;

    uses lphy.core.spi.Extension;
    provides lphy.core.spi.Extension with calibratedcpp.lphy.spi.CalibratedcppImpl;

}