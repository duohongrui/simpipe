#' @importFrom simutils fix_path
#' @importFrom dynwrap test_docker_installation
#' @importFrom babelwhale list_docker_images pull_container
#' @importFrom tidyr unite
#' @importFrom dplyr pull
method_execute_container_simulate <- function(
  parameters,
  method,
  other_prior,
  return_format,
  seed,
  verbose
){

  # Process file path-----------------------------------------------------------
  ## 1. Create a temp directory to store .rds file or .h5 file for Python
  temp_dir <- tempdir()
  ## 2. Set a path for the input data
  temp_input_path <- file.path(temp_dir, "input.rds") %>% simutils::fix_path()

  # Process datasets------------------------------------------------------------
  ## 1. Convert parameters into .rds and save to input path
  # simutils::write_h5files(data = ref_data, file_path = temp_input_path)
  saveRDS(parameters, file = temp_input_path)

  # Prepare the input parameters-----------------------------------------------
  ## 1. docker image working directory
  wd <- "/home/admin/"
  ## 2. local directory of the mount point
  local_path <- temp_dir %>% simutils::fix_path()
  ## 3. docker image directory of the mount point
  docker_path <- "/home/admin/docker_path"
  ## 4. verbose
  verbose <- verbose
  ## 5. args
  args <- NULL
  ## 6. command
  command <- NULL
  ## 7. container id
  ### (1. Check docker installation
  docker_install <- dynwrap::test_docker_installation()
  if(!docker_install){
    stop("Docker has not been installed or started! Please check it!")
  }
  ### (2. Check the installation of simpipe docker image
  images <- babelwhale::list_docker_images() %>%
    tidyr::unite("Repository", "Tag", sep = ":", col = "Image") %>%
    dplyr::pull("Image")

  if(!"duohongrui/simpipe:latest" %in% images){
    # If not, pull duohongrui/simpipe:latest
    babelwhale::pull_container(container_id = "duohongrui/simpipe:latest")
  }
  ### (3. docker container id
  container_id <- "duohongrui/simpipe"

  ## 8. step (estimation or simulation)
  step <- "simulation"

  ## 9. Return format (list, SingleCellExperiment, Seurat, h5ad)
  return_format <- return_format

  # Save command parameters
  input_meta <- list(container_id = container_id,
                     command = command,
                     args = args,
                     volums = paste0(local_path, ":", docker_path),
                     workspace = wd,
                     verbose = verbose,
                     step = step,
                     seed = seed,
                     method = method,
                     return_format = return_format,
                     other_prior = other_prior)
  saveRDS(input_meta, file.path(local_path, "input_meta.rds"))

  # Run container---------------------------------------------------------------
  output <- babelwhale::run(container_id = input_meta$container_id,
                            command = input_meta$command,
                            args = input_meta$args,
                            volumes = input_meta$volums,
                            workspace = input_meta$workspace,
                            verbose = input_meta$verbose,
                            debug = FALSE)

  # Get result------------------------------------------------------------------
  if(verbose){
    cat("Output is saved to ", local_path, "\n")
    cat("Attempting to read output into R\n")
  }

  ## return output
  readRDS(file = file.path(local_path, "output.rds"))

}
