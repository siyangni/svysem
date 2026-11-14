# Record the test methods on this package. UseMethod finds them in the S3
# table of stats::coef. The stats namespace registration itself is untouched.
registerS3method(
  "coef", "svysem_fake",
  function(object, ...) object$coefficients,
  envir = asNamespace("svysem")
)
registerS3method(
  "vcov", "svysem_fake",
  function(object, ...) object$variance,
  envir = asNamespace("svysem")
)

fits_from_qu <- function(Qhat, Uhat) {
  nms <- rownames(Qhat)
  lapply(seq_len(ncol(Qhat)), function(i) {
    b <- Qhat[, i]
    names(b) <- nms
    V <- Uhat[, , i, drop = TRUE]
    V <- as.matrix(V)
    dimnames(V) <- list(nms, nms)
    structure(list(coefficients = b, variance = V), class = "svysem_fake")
  })
}

# UseMethod finds generic.class in the package namespace before the S3
# methods table. Replacing only the table does not hide the real method
# once survey.design2 methods exist. Put the stub in both places.
# assignInNamespace() re-reads the NAMESPACE S3 table and errors here, so
# the binding is replaced directly.
register_design_method <- function(generic, fn) {
  ns <- asNamespace("svysem")
  name <- paste(generic, "survey.design2", sep = ".")
  had_fn <- exists(name, envir = ns, inherits = FALSE)
  old_fn <- if (had_fn) get(name, envir = ns) else NULL
  locked <- had_fn && bindingIsLocked(name, ns)
  if (locked) {
    unlockBinding(name, ns)
  }
  if (had_fn) {
    assign(name, fn, envir = ns)
  }
  if (locked) {
    lockBinding(name, ns)
  }
  registerS3method(generic, "survey.design2", fn, envir = ns)
  function() {
    if (!had_fn) {
      return(invisible(NULL))
    }
    if (bindingIsLocked(name, ns)) {
      unlockBinding(name, ns)
    }
    assign(name, old_fn, envir = ns)
    if (locked) {
      lockBinding(name, ns)
    }
    registerS3method(generic, "survey.design2", old_fn, envir = ns)
  }
}
