#' Estimate parameters from Real Datasets
#'
#' This is the core function of estimating parameters from real single-cell RNA
#' sequencing datasets. Users can get many estimation results by pointing out
#' more than two methods and inputting many real datasets.
#'
#' @param ref_data A matrix for one dataset or a list of datasets with their own
#' names.
#' @param method A character or a string of methods.
#' @param other_prior A list with names of certain parameters. Some methods need
#' extra parameters to execute the estimation step, so you must input them. You
#' will get the instructions if the error occurs.
#' @param seed A random see.
#' @param verbose Logical. Whether to return messages or not.
#' @param use_docker Logical. Default is FALSE. Whether to execute the estimation
#' step by establishing a Docker container or not. Suggest that you would better
#' use Docker to execute the step in order to avoid the issues of coding environment.
#'
#' @return A list of estimated parameters and detection results
#' @importFrom tidyr crossing
#' @importFrom simmethods get_method
#' @importFrom simutils write_h5files random_seed
#' @importFrom assertthat assert_that
#' @importFrom purrr map
#' @importFrom stats setNames
#' @export
#'
estimate_parameters <- function(
  ref_data,
  method,
  other_prior = list(),
  seed = simutils::random_seed(),
  verbose = TRUE,
  use_docker = FALSE
){
  # Check-----------------------------------------------------------------------
  if(is.matrix(ref_data)){
    ref_data <- list(ref_data = ref_data)
  }
  assertthat::assert_that(is.list(ref_data))

  # Prepare methods-------------------------------------------------------------
  all_methods <- simmethods::get_method()
  if(method == "all"){
    method <- names(all_methods)
  }
  assertthat::assert_that(all(method %in% names(all_methods)))

  # Estimation design-----------------------------------------------------------
  design <- tidyr::crossing(dataset_id = names(ref_data),
                            method_id = method)
  # Run methods with each dataset-----------------------------------------------
  result <- purrr::map(
    .x = seq_len(nrow(design)),
    .f = function(er) {
      seed <- ifelse(is.null(seed), random_seed(), seed)
      if(use_docker){
        method_execute_container_estimate(
          ref_data = ref_data[[design$dataset_id[er]]],
          method = design$method_id[er],
          other_prior = other_prior,
          seed = seed,
          verbose = verbose)
      }else{
        method_execute_function_estimate(
          ref_data = ref_data[[design$dataset_id[er]]],
          method = design$method_id[er],
          other_prior = other_prior,
          seed = seed,
          verbose = verbose,
          env = env
        )
      }
    }
  )
  result_names <- paste0(design$dataset_id, "_", design$method_id)
  result <- stats::setNames(result, result_names)
  return(result)
}



a <- matrix(rpois(n = 10^6, lambda = 0.5), nrow = 1000)
colnames(a) <- paste0("cell_", 1:ncol(a))
rownames(a) <- paste0("gene_", 1:nrow(a))

b <- matrix(rpois(n = 10^6, lambda = 0.1), nrow = 1000)
colnames(b) <- paste0("cell_", 1:ncol(b))
rownames(b) <- paste0("gene_", 1:nrow(b))

ref_data <- list(a = a, b = b)
