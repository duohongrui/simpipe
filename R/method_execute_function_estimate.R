method_execute_function_estimate <- function(
  ref_data,
  method,
  other_prior,
  seed,
  verbose
){
  # Select the right method
  env <- asNamespace("simmethods")
  right_method <- paste0(method, "_estimation")
  assign(right_method, get(right_method, envir = env))

  # Match arguments existing in the method function
  arguments <- dplyr::lst(ref_data, other_prior, seed, verbose)
  arguments <- arguments[intersect(names(arguments), names(formals(right_method)))]

  # Execute estimation step by do.call and pass the arguments to the estimation function
  estimate_output <- do.call(right_method, arguments)

  # Output
  estimate_output

}

