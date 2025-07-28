module jdplus.revisions.base.r {

    requires static lombok;
    requires static nbbrd.design;
    requires static nbbrd.service;
    requires static org.jspecify;

    requires transitive jdplus.revisions.base.api;
    requires jdplus.toolkit.base.api;
    requires jdplus.revisions.base.core;
    requires jdplus.toolkit.base.core;

    exports jdplus.revisions.base.r;
}