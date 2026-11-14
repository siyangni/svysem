# Registered here so UseMethod finds them from the package namespace.
registerS3method(
  "coef", "svysem_fake",
  function(object, ...) object$coefficients,
  envir = asNamespace("stats")
)
registerS3method(
  "vcov", "svysem_fake",
  function(object, ...) object$variance,
  envir = asNamespace("stats")
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

# Restore whatever svysem method for survey.design2 was registered, if any.
register_design_method <- function(generic, fn) {
  ns <- asNamespace("svysem")
  key <- paste(generic, "survey.design2", sep = ".")
  old <- NULL
  if (exists(".__S3MethodsTable__.", envir = ns, inherits = FALSE)) {
    table <- get(".__S3MethodsTable__.", envir = ns)
    if (exists(key, envir = table, inherits = FALSE)) {
      old <- get(key, envir = table)
    }
  }
  registerS3method(generic, "survey.design2", fn, envir = ns)
  function() {
    table <- get(".__S3MethodsTable__.", envir = ns)
    if (is.null(old)) {
      if (exists(key, envir = table, inherits = FALSE)) {
        rm(list = key, envir = table)
      }
    } else {
      assign(key, old, envir = table)
    }
  }
}
