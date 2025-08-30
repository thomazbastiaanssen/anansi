#' @rawNamespace if (getRversion() < "4.3.0") importFrom("S7", "@")
NULL


rlang::on_load({

    S7::S4_register(anansi:::AnansiWeb);
    S7::S4_register(anansi:::MultiFactor);

    S7::methods_register();

    S7::method(dimnames, AnansiWeb) <- function(x, ...) `dimnames.anansi::AnansiWeb`(x);
    S7::method(dim, AnansiWeb)  <- function(x, ...) `dim.anansi::AnansiWeb`(x);
    S7::method(names, AnansiWeb) <- function(x, ...) `names.anansi::AnansiWeb`(x);

    S7::method(dimnames, MultiFactor) <- function(x) `dimnames.anansi::MultiFactor`(x);
    S7::method(dim, MultiFactor)  <- function(x) `dim.anansi::MultiFactor`(x);
    S7::method(names, MultiFactor) <- function(x) `names.anansi::MultiFactor`(x);
    S7::method(levels, MultiFactor) <- function(x) `levels.anansi::MultiFactor`(x);
    S7::method(`levels<-`, MultiFactor) <- function(x, value) `levels<-.anansi::MultiFactor`(x, value);
    S7::method(droplevels, MultiFactor) <- function(x, exclude = NULL, select = NULL) `droplevels.anansi::MultiFactor`(x, exclude = exclude, select = exclude);


    S7::method(`[`, MultiFactor) <- function(x, ...) `[.anansi::MultiFactor`(x, ...);

    S7::method(`[[`, MultiFactor) <- function(x, ...) `[[.anansi::MultiFactor`(x, ...);

    S7::method(`[<-`, MultiFactor) <- function(x, ..., value) `[<-.anansi::MultiFactor`(x, ..., value);

    S7::method(`[[<-`, MultiFactor) <- function(x, ..., value) `[[<-.anansi::MultiFactor`(x, ..., value)

})

.onLoad <- function(...) {
    rlang::run_on_load()
}
