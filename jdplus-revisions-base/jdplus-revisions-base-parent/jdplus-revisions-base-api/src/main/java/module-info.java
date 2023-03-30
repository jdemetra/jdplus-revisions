module jdplus.revisions.base.api {

    requires static lombok;
    requires static nbbrd.design;
    requires static nbbrd.service;
    requires static org.checkerframework.checker.qual;
    requires jdplus.toolkit.base.api;

    exports jdplus.revisions.base.api.parametric;
    exports jdplus.revisions.base.api.timeseries;
}