module calibratedcpp.lphy {
    requires transitive lphy.base;
    requires jdk.jfr;

    exports calibratedcpp.lphy.tree;
    exports calibratedcpp.lphy.prior;
    exports calibratedcpp.lphy.util;

    uses lphy.core.spi.Extension;
    provides lphy.core.spi.Extension with calibratedcpp.lphy.spi.CalibratedcppImpl;

}