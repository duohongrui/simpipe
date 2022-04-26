#' Get A random Seed
#'
#' @export
#'
#' @examples
#' seed <- random_seed()
random_seed <- function(){
  sample.int(.Machine$integer.max, 1)
}

simulate_datasets <- function(method,
                              parameters,
                              seed = random_seed()){

  NULL

}
#' @examples
#' # Generate a reference data
#' set.seed(1)
#' ref_data <- matrix(rpois(n = 2500, lambda = 2), nrow = 50)
#' rownames(ref_data) <- paste0("cell_", 1:ncol(ref_data))
#' colnames(ref_data) <- paste0("gene_", 1:nrow(ref_data))
sim_splat <- function(data, nCells, nGenes, ...){

  NULL

}


# set.seed(1)
# ref_data <- matrix(rpois(n = 10^7, lambda = 2), nrow = 10000)
# colnames(ref_data) <- paste0("cell_", 1:ncol(ref_data))
# rownames(ref_data) <- paste0("gene_", 1:nrow(ref_data))
