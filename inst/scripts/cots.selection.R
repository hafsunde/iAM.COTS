#' Create a canonical pair identifier
#'
#' Builds a stable pair id using full identifiers (no substring truncation).
#'
#' @param x Vector of ids for person 1.
#' @param y Vector of ids for person 2.
#' @param sep Separator used in the pair id.
#'
#' @return Character vector with canonicalized pair ids.
#'
#' @keywords internal
#' @noRd
cots_pair_id <- function(x, y, sep = "_") {
  x_chr <- as.character(x)
  y_chr <- as.character(y)
  first <- ifelse(x_chr <= y_chr, x_chr, y_chr)
  second <- ifelse(x_chr <= y_chr, y_chr, x_chr)
  paste(first, second, sep = sep)
}


#' Deterministic greedy selector for non-overlapping rows
#'
#' Keeps rows so that identifiers listed in `id_cols` are not reused across
#' selected rows.
#'
#' @param dt Data.frame/data.table with candidate rows.
#' @param id_cols Character vector of columns containing ids that must be unique
#'   across selected rows.
#' @param order_cols Optional character vector used for deterministic ordering
#'   before selection.
#' @param deterministic Logical. If `FALSE`, randomizes row order before greedy
#'   selection (useful for diagnostics).
#' @param prioritize_n_ids Logical. If `TRUE`, rows with more non-missing IDs
#'   in `id_cols` are prioritized before greedy selection.
#'
#' @return A data.table containing selected non-overlapping rows.
#'
#' @keywords internal
#' @noRd
cots_select_nonoverlap <- function(dt,
                                   id_cols,
                                   order_cols = NULL,
                                   deterministic = TRUE,
                                   prioritize_n_ids = TRUE) {
  x <- data.table::as.data.table(data.table::copy(dt))
  if (!nrow(x)) {
    return(x)
  }

  missing_cols <- setdiff(id_cols, names(x))
  if (length(missing_cols)) {
    stop("Missing id columns: ", paste(missing_cols, collapse = ", "))
  }

  if (prioritize_n_ids) {
    x[, ..n_present_ids := rowSums(!is.na(.SD) & .SD != "", na.rm = TRUE), .SDcols = id_cols]
  }

  order_cols <- order_cols[order_cols %in% names(x)]
  if (prioritize_n_ids && length(order_cols)) {
    data.table::setorderv(
      x,
      cols = c("..n_present_ids", order_cols),
      order = c(-1L, rep(1L, length(order_cols))),
      na.last = TRUE
    )
  } else if (prioritize_n_ids) {
    data.table::setorderv(x, cols = "..n_present_ids", order = -1L, na.last = TRUE)
  } else if (length(order_cols)) {
    data.table::setorderv(x, cols = order_cols, na.last = TRUE)
  } else if (!deterministic) {
    x <- x[sample.int(nrow(x))]
  }

  # Build per-row ID bundles using only parent columns (id_cols), then map the
  # IDs to dense integer codes for fast O(1) membership checks during the greedy
  # pass.
  id_mat <- as.matrix(x[, ..id_cols])
  id_vec <- as.character(as.vector(t(id_mat)))
  row_id <- rep(seq_len(nrow(x)), each = ncol(id_mat))
  valid <- !is.na(id_vec) & nzchar(id_vec)

  if (!any(valid)) {
    return(x)
  }

  id_vec <- id_vec[valid]
  row_id <- row_id[valid]
  unique_ids <- unique(id_vec)
  key <- match(id_vec, unique_ids)
  idx_split <- split(key, row_id)

  idx_list <- vector("list", nrow(x))
  idx_names <- as.integer(names(idx_split))
  idx_list[idx_names] <- idx_split
  idx_list <- lapply(idx_list, function(v) {
    if (is.null(v)) integer() else as.integer(unique(v))
  })

  keep <- logical(nrow(x))
  seen <- rep(FALSE, length(unique_ids))

  for (i in seq_along(idx_list)) {
    v <- idx_list[[i]]
    if (!length(v)) {
      keep[i] <- TRUE
      next
    }
    if (any(seen[v])) {
      next
    }
    keep[i] <- TRUE
    seen[v] <- TRUE
  }

  out <- x[keep]
  if ("..n_present_ids" %in% names(out)) {
    out[, ..n_present_ids := NULL]
  }
  out
}


#' Collect non-missing ids from selected columns
#'
#' @param dt Data.frame/data.table.
#' @param id_cols Character vector of columns containing ids.
#'
#' @return Character vector of unique non-missing ids.
#'
#' @keywords internal
#' @noRd
cots_collect_ids <- function(dt, id_cols) {
  x <- data.table::as.data.table(dt)
  present <- intersect(id_cols, names(x))
  if (!length(present) || !nrow(x)) {
    return(character())
  }
  vals <- as.character(unlist(x[, ..present], use.names = FALSE))
  unique(vals[!is.na(vals) & nzchar(vals)])
}
