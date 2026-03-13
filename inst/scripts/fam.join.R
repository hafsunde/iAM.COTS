#' Attach Nuclear Families to Related Pairs
#'
#' Splits `extended` by `fam_link` and joins each subset to nuclear family
#' identifiers using `fam.join.inner()`.
#'
#' @param extended Data.frame/data.table with related pairs and a `fam_link`
#'   column.
#' @param nuclear Data.frame/data.table of nuclear families.
#'
#' @return A data.table with family IDs and group labels.
#'
#' @keywords internal
#' @noRd
fam.join <- function(extended, nuclear) {
  ext <- data.table::as.data.table(data.table::copy(extended))
  if (!"fam_link" %in% names(ext)) {
    warning("fam_link column does not exist. If one family structure, use fam.join.inner().")
    return(NULL)
  }

  parts <- lapply(split(ext, ext$fam_link), fam.join.inner, nuclear = nuclear)
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) {
    return(data.table::as.data.table(data.frame()))
  }
  data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
}
