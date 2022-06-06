simulate_datasets <- function(
  method,
  parameters,
  other_prior = list(),
  seed = simutils::random_seed(),
  return_format = "SingleCellExperiment",
  verbose = TRUE,
  use_docker = FALSE
){

  # Check-----------------------------------------------------------------------
  assertthat::assert_that(is.list(parameters))

  # Prepare methods-------------------------------------------------------------
  all_methods <- simmethods::get_method()
  if(method == "all"){
    method <- names(all_methods)
  }
  assertthat::assert_that(all(method %in% names(all_methods)))

  # Run methods with each estimation and each method----------------------------
  id_name <- names(estimate_output)

  result <- purrr::map(
    .x = id_name,
    .f = function(id) {
      seed <- ifelse(is.null(seed), random_seed(), seed)
      method_name <- stringr::str_split(id, pattern = "_", simplify = T)[2]

      if(use_docker){
        # method_execute_container_simulate(
        #   parameters = parameters[[id]][["estimate_result"]],
        #   method = method_name,
        #   other_prior = other_prior,
        #   return_format = return_format,
        #   seed = seed,
        #   verbose = verbose)
      }else{
        method_execute_function_simulate(
          parameters = parameters[[id]][["estimate_result"]],
          method = method_name,
          other_prior = other_prior,
          return_format = return_format,
          seed = seed,
          verbose = verbose
        )
      }
    }
  )
  result <- stats::setNames(result, id_name)
  return(result)
}
