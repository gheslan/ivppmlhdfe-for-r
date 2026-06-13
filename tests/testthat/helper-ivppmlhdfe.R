# Shared test helpers ------------------------------------------------------

# fixest is an Imports dependency (used for HDFE demeaning), so fepois is
# always available; we use it as the independent numerical benchmark for
# plain PPML point estimates and (dof-aligned) sandwich matrices.
have_fixest <- requireNamespace("fixest", quietly = TRUE)

skip_no_fixest <- function() {
  testthat::skip_if_not(have_fixest, "fixest not available")
}

get_class <- function(cl) {
  e <- new.env()
  utils::data(list = paste0("ivppmlhdfe_Class", cl), package = "ivppmlhdfe",
              envir = e)
  get(paste0("ivppmlhdfe_Class", cl), envir = e)
}
