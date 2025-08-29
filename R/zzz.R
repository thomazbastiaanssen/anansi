#' @rawNamespace if (getRversion() < "4.3.0") importFrom("S7", "@")
NULL


rlang::on_load({

    S7::S4_register(anansi:::AnansiWeb);
    S7::S4_register(anansi:::MultiFactor);

    S7::methods_register();

    S7::method(`[`, MultiFactor) <- function(x, ...) `[.anansi::MultiFactor`(x, ...);

    S7::method(`[[`, MultiFactor) <- function(x, ...) `[[.anansi::MultiFactor`(x, ...);

    S7::method(`[<-`, MultiFactor) <- function(x, ..., value) `[<-.anansi::MultiFactor`(x, ..., value);

    S7::method(`[[<-`, MultiFactor) <- function(x, ..., value) `[[<-.anansi::MultiFactor`(x, ..., value)

})

.onLoad <- function(...) {
    rlang::run_on_load()
}
