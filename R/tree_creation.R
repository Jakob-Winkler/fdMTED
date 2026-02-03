
tf_where_all_value <- function (f, cond, arg)
{
  if (missing(arg)) {
    arg <- tf_arg(f)
  }
  # assert_arg(arg, f) TODO should not be removed
  cond_call <- substitute(cond)
  parent <- parent.frame()
  where_at <- map(f[, arg, matrix = FALSE], function(x) subset(x,eval(cond_call, envir = x, enclos = parent))[["value"]])
  where_at[is.na(f)] <- NA

  return(where_at)

}

prep_extrem <- function(data){
  # TODO checks
  if (data[1] > data[2]) data <- data[-1]
  last_idx <- length(data)
  if (data[last_idx - 1]  < data[last_idx]) data <- data[-last_idx]
  return(data)
}
#' @export
get_prep_extrems <- function(data) {
  extrems <- tf_where_all_value(
    data,
    (sign(c(diff(value)[1], diff(value))) !=
       sign(c(diff(value), tail(diff(value), 1))))|(arg == min(arg))|(arg == max(arg))
  )
  return(lapply(extrems,prep_extrem))
}


# I think it works, not validated, no checks
#' @export
tf_mergetree <- function(data){
  #browser()
  number_nodes <- length(data)

  # Pre-allocat output list
  out <- replicate(number_nodes, list(parent = NULL, childs = NULL, weight = 0), simplify = FALSE)
  names(out) <- 1:number_nodes
  # Handel first minima
  out[[1]]["parent"] <- 2
  out[[1]]["weight"] <- data[2]-data[1]

  out[[2]][["childs"]] <- 1

  # Init stack
  stack <- vector("list", number_nodes)
  top_stack <- 1
  # c(index, height)
  stack[[top_stack]] <- c(idx = 2, height = data[2])

  repeat {
    new_idx <- stack[[top_stack]][["idx"]] + 2
    if (new_idx > number_nodes) {
      #browser()
      # output edit

      # Edit on minima
      out[[stack[[top_stack]]["idx"] + 1]]["parent"] <- stack[[top_stack]]["idx"]
      out[[stack[[top_stack]]["idx"] + 1]]["weight"] <- stack[[top_stack]]["height"]-data[stack[[top_stack]]["idx"] + 1]

      # Edit on maxima
      out[[stack[[top_stack]]["idx"]]][["childs"]] <- c(out[[stack[[top_stack]]["idx"]]][["childs"]], stack[[top_stack]][["idx"]] + 1)
      # Edits based on stack
      while(top_stack != 1) {

        # Pop
        merge_idx <- stack[[top_stack]][["idx"]]
        stack[[top_stack]] <- list(NULL)
        top_stack <- top_stack-1

        # Merge merge_idx with canidet
        out[[merge_idx]]["parent"] <- stack[[top_stack]]["idx"]
        out[[merge_idx]]["weight"] <- stack[[top_stack]]["height"]-data[merge_idx]

        out[[stack[[top_stack]]["idx"]]][["childs"]] <- c(out[[stack[[top_stack]]["idx"]]][["childs"]], merge_idx)

      }
      break

    }
    if (stack[[top_stack]]["height"] < data[new_idx]) {

      # output edit

      # Edit on minima
      out[[stack[[top_stack]]["idx"] + 1]]["parent"] <- stack[[top_stack]]["idx"]
      out[[stack[[top_stack]]["idx"] + 1]]["weight"] <- stack[[top_stack]]["height"]-data[stack[[top_stack]]["idx"] + 1]

      # Edit on maxima
      out[[stack[[top_stack]]["idx"]]][["childs"]] <- c(out[[stack[[top_stack]]["idx"]]][["childs"]], stack[[top_stack]][["idx"]] + 1)

      # Edits based on stack
      loop_end <- F
      while(top_stack != 1) {

        # Pop
        merge_idx <- stack[[top_stack]][["idx"]]
        stack[[top_stack]] <- list(NULL)
        top_stack <- top_stack-1

        if (stack[[top_stack]]["height"] < data[new_idx]) {
          # Merge merge_idx with canidet
          out[[merge_idx]]["parent"] <- stack[[top_stack]]["idx"]
          out[[merge_idx]]["weight"] <- stack[[top_stack]]["height"]-data[merge_idx]

          out[[stack[[top_stack]]["idx"]]][["childs"]] <- c(out[[stack[[top_stack]]["idx"]]][["childs"]], merge_idx)

        } else {
          # Merge merge_idx with new_idx
          out[[merge_idx]]["parent"] <- new_idx
          out[[merge_idx]]["weight"] <- data[new_idx] - data[merge_idx]

          out[[new_idx]][["childs"]] <- c(out[[new_idx]][["childs"]], merge_idx)

          # new_idx add to stack
          # Push
          top_stack <- top_stack + 1
          stack[[top_stack]] <- c(idx = new_idx, height = data[new_idx])
          loop_end <- T
          break # find new node
        }
      }
      if (loop_end) next
      # If top is root, new is then root
      #+ Merge top stack with new
      out[[stack[[top_stack]]["idx"]]]["parent"] <- new_idx
      out[[stack[[top_stack]]["idx"]]]["weight"] <- data[new_idx]-stack[[top_stack]]["height"]

      out[[new_idx]][["childs"]] <- c(out[[new_idx]][["childs"]], stack[[top_stack]][["idx"]])

      #+ New is top stack
      stack[[top_stack]] <- c(idx = new_idx, height = data[new_idx])
    }
    else {
      # output edit

      # Edit on minima
      out[[stack[[top_stack]]["idx"] + 1]]["parent"] <- new_idx
      out[[stack[[top_stack]]["idx"] + 1]]["weight"] <- data[new_idx] - data[stack[[top_stack]]["idx"] + 1]

      # Edit on maxima
      out[[new_idx]][["childs"]] <- c(out[[new_idx]][["childs"]], stack[[top_stack]][["idx"]] + 1)

      # Add new node to stack
      # Push
      top_stack <- top_stack + 1
      stack[[top_stack]] <- c(idx = new_idx, height = data[new_idx])
    }

  }

  return(out)
}
#' @export
tf_mergetree_list <- function(data) {
  lapply(data,tf_mergetree)
}


node_as_vector <- function(node) {
  c(ifelse(is.null(node[["parent"]]), NA, node[["parent"]]),
    ifelse(is.null(node[["childs"]]), NA, node[["childs"]][1]),
    ifelse(is.null(node[["childs"]]), NA, node[["childs"]][2]),
    node[["weight"]])
}
#' @export
mergetree_as_vector <- function(tree) { # This returns index for cpp
  number_nodes <- length(tree)
  parent_vec <- vector(mode = "integer", number_nodes)
  child_1_vec <- vector(mode = "integer", number_nodes)
  child_2_vec <- vector(mode = "integer", number_nodes)
  weight_vec <- vector(mode = "numeric", number_nodes)
  # The first index in reserved for the root, so all are shifted until the root is found
  for (x in 1:number_nodes) {
    node_vec <- node_as_vector(tree[[x]])
    idx <- x
    parent_vec[idx] <- node_vec[1] - 1
    child_1_vec[idx] <- node_vec[2] - 1
    child_2_vec[idx] <- node_vec[3] - 1
    weight_vec[idx] <- node_vec[4]
  }

  return(list(parent_vec = parent_vec,
              child_1_vec = child_1_vec,
              child_2_vec = child_2_vec,
              weight_vec = weight_vec))
}
