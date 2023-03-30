module jdplus.revisions.base.core {

    requires static lombok;
    requires static nbbrd.design;
    requires static nbbrd.service;
    requires static org.checkerframework.checker.qual;

    requires transitive jdplus.revisions.base.api;
    requires jdplus.toolkit.base.api;
    requires jdplus.toolkit.base.core;

    exports jdplus.revisions.base.core.parametric;
}