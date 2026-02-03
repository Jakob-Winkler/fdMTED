
#' @export
BLP_solver_Cpp <- function(cost_Vec, constraint_Matrix, off_set, idx_tree1, idx_tree2, W) {
    model            <- list()
    model$A          <- constraint_Matrix
    model$obj        <- cost_Vec
    model$modelsense <- 'min'
    model$rhs        <- rep(1, nrow(constraint_Matrix))
    model$sense      <- c('<')
    model$objcon     <- off_set
    model$lb <- 0
    model$ub <- 1
    model$vtype <- "B"


    result <- gurobi(model, params=list(LogToConsole = 0, Threads = 1)) # Is multithreading worth it?

    result$objval
}




#' Split a string
#'
#' @param x A character vector with one element.
#' @param split What to split on.
#'
#' @return A character vector.
#' @export
#'
#' @examples
#' x <- "alfa,bravo,charlie,delta"
#' strsplit1(x, split = ",")
TED <- function(tree1, tree2) {

  tree_edit_distance(tree1$parent_vec, tree1$child_1_vec, tree1$child_2_vec, tree1$weight_vec,
                     tree2$parent_vec, tree2$child_1_vec, tree2$child_2_vec, tree2$weight_vec)

}

#' @export
dist_matrix <- function(x, dist_fun, diag_val = 0, verbose = TRUE, ...) {
  stopifnot(is.list(x), length(x) >= 1L, is.function(dist_fun))

  n <- length(x)
  nm <- names(x)
  if (is.null(nm) || any(nm == "")) nm <- as.character(seq_len(n))

  # allocate result
  D <- matrix(NA_real_, n, n, dimnames = list(nm, nm))
  diag(D) <- diag_val

  if (n < 2L) return(D)

  # number of unique pairs
  n_pairs <- n * (n - 1) / 2

  pb <- NULL
  if (isTRUE(verbose)) {
    pb <- utils::txtProgressBar(min = 0, max = n_pairs, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }

  k <- 0L
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      d <- dist_fun(x[[i]], x[[j]], ...)
      D[i, j] <- d
      D[j, i] <- d

      if (isTRUE(verbose)) {
        k <- k + 1L
        utils::setTxtProgressBar(pb, k)
      }
    }
  }

  D
}

#' @export
dist_matrix_parallel <- function(x, dist_fun,
                                 diag_val = 0,
                                 ncores = max(1L, parallel::detectCores() - 1L),
                                 verbose = TRUE, ...) {
  stopifnot(is.list(x), length(x) >= 1L, is.function(dist_fun))

  n <- length(x)
  nm <- names(x)
  if (is.null(nm) || any(nm == "")) nm <- as.character(seq_len(n))

  D <- matrix(NA_real_, n, n, dimnames = list(nm, nm))
  diag(D) <- diag_val
  if (n < 2L) return(D)

  pairs <- utils::combn(n, 2, simplify = FALSE)
  n_pairs <- length(pairs)

  cl <- parallel::makeCluster(min(ncores, n))
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterEvalQ(cl, {
    NULL
  })

  # Export needed objects
  parallel::clusterExport(cl, varlist = c("x", "dist_fun"), envir = environment())

  pb <- NULL
  if (isTRUE(verbose)) {
    pb <- utils::txtProgressBar(min = 0, max = n_pairs, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }

  # Use load-balanced apply; update progress as chunks return
  res <- vector("list", n_pairs)
  idx <- 0L

  # parLapplyLB returns results in order of input, but finishes chunks out of order internally.
  # So we do a simple loop updating the bar after it returns (still gives feedback).
  res <- parallel::parLapplyLB(cl, pairs, function(ij) {
    i <- ij[1]; j <- ij[2]
    d <- dist_fun(x[[i]], x[[j]], ...)
    list(i=i, j=j, d=d)
  })

  if (isTRUE(verbose)) utils::setTxtProgressBar(pb, n_pairs)

  for (r in res) {
    D[r$i, r$j] <- r$d
    D[r$j, r$i] <- r$d
  }

  D
}
