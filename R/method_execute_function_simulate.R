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

  # Change parameters
  parameters <- set_parameters(parameters = parameters,
                               other_prior = other_prior,
                               method)

  # Match arguments existing in the method function
  arguments <- dplyr::lst(parameters,
                          other_prior,
                          return_format,
                          seed,
                          verbose)
  arguments <- arguments[intersect(names(arguments), names(formals(right_method)))]

  # Execute estimation step by do.call and pass the arguments to the estimation function
  simulate_output <- do.call(right_method, arguments)

  # Output
  simulate_output
}

