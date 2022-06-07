method_execute_function_simulate <- function(
  method,
  parameters,
  other_prior,
  return_format,
  seed,
  verbose
){
  # Select the right method
  env <- asNamespace("simmethods")
  right_method <- paste0(method, "_simulation")
  assign(right_method, get(right_method, envir = env))

  # When the parameter is NULL
  if(is.null(parameters)){
    cat("Parameters will be set by default")
    parameters <- simutils::default_parameters(method)
  }

  # Change parameters except for "n"
  if(!is.null(other_prior) & is.null(other_prior[["n"]])){
    parameters <- simutils::set_parameters(parameters = parameters,
                                           other_prior = other_prior,
                                           method = method)
  }
  # Match arguments existing in the method function
  arguments <- dplyr::lst(parameters,
                          other_prior,
                          return_format,
                          verbose,
                          seed)
  arguments <- arguments[intersect(names(arguments), names(formals(right_method)))]

  # Execute estimation step by do.call and pass the arguments to the estimation function
  simulate_output <- do.call(right_method, arguments)

  # Output
  simulate_output
}

