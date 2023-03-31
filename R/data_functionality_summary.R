#' Summary
#'
#' @param data The matrix or data frame of gene expression profile or the output by `simpipe::simulate_datasets`.
#' @param group The characters of group names which cells belong to. Needed if data is a matrix or a dataframe and you want to evaluate the groups.
#' @param batch The characters of batch names which cells belong to. Needed if data is a matrix or a dataframe and you want to evaluate the batches
#' @param k k-nearest neighborhoods of cells.
#' @param DEGs A list of DEGs with the names of `xxxvsxxx`. Note that the names of DEGs are in the rownames of the matrix or the dataframe and the names of `xxx` is in the `group` characters. Needed if data is a matrix or a dataframe and you want to evaluate the DEGs.
#' @param verbose Whether the process massages are returned.
#' @importFrom SingleCellExperiment counts colData rowData
#' @importFrom simutils calculate_DEGs_properties
#'
#'
#' @return A list of three aspects of the data
#' @export
data_functionality_summary <- function(
  data,
  group = NULL,
  batch = NULL,
  k = NULL,
  DEGs = NULL,
  verbose = NULL

){
  ### data check and other imformation of cells and genes
  if("simpipe_simulation" %in% class(data)){
    data_number <- length(data)

    ## Check the return format
    if(class(data[[1]][["simulate_result"]]) == "SingleCellExperiment"){
      count_matrices <- purrr::map(
        .x = 1:data_number,
        .f = function(id){
          if(!is.matrix(SingleCellExperiment::counts(data[[id]]$simulate_result))){
            as.matrix(SingleCellExperiment::counts(data[[id]]$simulate_result))
          }else{
            SingleCellExperiment::counts(data[[id]]$simulate_result)
          }
        }
      ) %>% stats::setNames(names(data))

      data_colnames <- purrr::map(
        .x = 1:data_number,
        .f = function(id){
          colnames(data[[id]]$simulate_result)
        }
      ) %>% stats::setNames(names(data))

      data_rownames <- purrr::map(
        .x = 1:data_number,
        .f = function(id){
          rownames(data[[id]]$simulate_result)
        }
      ) %>% stats::setNames(names(data))

      data_coldata <- purrr::map(
        .x = 1:data_number,
        .f = function(id){
          as.data.frame(SingleCellExperiment::colData(data[[id]][["simulate_result"]]))
        }
      ) %>% stats::setNames(names(data))

      data_rowdata <- purrr::map(
        .x = 1:data_number,
        .f = function(id){
          as.data.frame(SingleCellExperiment::rowData(data[[id]][["simulate_result"]]))
        }
      ) %>% stats::setNames(names(data))
    }

    ## Seurat Object
    if(class(data[[1]][["simulate_result"]]) == "Seurat"){
      count_matrices <- purrr::map(
        .x = 1:data_number,
        .f = function(id){
          if(!is.matrix(data[[id]]$simulate_result@assays$originalexp@counts)){
            as.matrix(data[[id]]$simulate_result@assays$originalexp@counts)
          }else{
            data[[id]]$simulate_result@assays$originalexp@counts
          }
        }
      ) %>% stats::setNames(names(data))

      data_colnames <- purrr::map(
        .x = 1:data_number,
        .f = function(id){
          colnames(data[[id]]$simulate_result)
        }
      ) %>% stats::setNames(names(data))

      data_rownames <- purrr::map(
        .x = 1:data_number,
        .f = function(id){
          rownames(data[[id]]$simulate_result)
        }
      ) %>% stats::setNames(names(data))

      data_coldata <- purrr::map(
        .x = 1:data_number,
        .f = function(id){
          as.data.frame(data[[id]][["simulate_result"]]@meta.data)
        }
      ) %>% stats::setNames(names(data))

      data_rowdata <- purrr::map(
        .x = 1:data_number,
        .f = function(id){
          as.data.frame(data[[id]][["simulate_result"]]@assays$originalexp@meta.features)
        }
      ) %>% stats::setNames(names(data))
    }

    ## list
    if(class(data[[1]][["simulate_result"]]) == "list" &
       length(data[[1]][["simulate_result"]]) == 3){
      count_matrices <- purrr::map(
        .x = 1:data_number,
        .f = function(id){
          if(!is.matrix(data[[id]]$simulate_result$count_data)){
            as.matrix(data[[id]]$simulate_result$count_data)
          }else{
            data[[id]]$simulate_result$count_data
          }
        }
      ) %>% stats::setNames(names(data))

      data_colnames <- purrr::map(
        .x = 1:data_number,
        .f = function(id){
          colnames(data[[id]]$simulate_result$count_data)
        }
      ) %>% stats::setNames(names(data))

      data_rownames <- purrr::map(
        .x = 1:data_number,
        .f = function(id){
          rownames(data[[id]]$simulate_result$count_data)
        }
      ) %>% stats::setNames(names(data))

      data_coldata <- purrr::map(
        .x = 1:data_number,
        .f = function(id){
          as.data.frame(data[[id]][["simulate_result"]]$col_meta)
        }
      ) %>% stats::setNames(names(data))

      data_rowdata <- purrr::map(
        .x = 1:data_number,
        .f = function(id){
          as.data.frame(data[[id]][["simulate_result"]]$row_meta)
        }
      ) %>% stats::setNames(names(data))
    }

    ### Group information
    group <- purrr::map(
      .x = 1:data_number,
      .f = function(id){
        col_data <- data_coldata[[id]]
        if("plate" %in% colnames(col_data)){
          as.character(col_data$"plate")
        }else if("group" %in% colnames(col_data)){
          as.character(col_data$"group")
        }else{
          NULL
        }
      }
    ) %>% stats::setNames(names(data))
    ### Batch information
    batch <- purrr::map(
      .x = 1:data_number,
      .f = function(id){
        col_data <- data_coldata[[id]]
        if("batch" %in% colnames(col_data)){
          col_data$"batch"
        }else{
          NULL
        }
      }
    ) %>% stats::setNames(names(data))
    ### Group information
    DEGs <- purrr::map(
      .x = 1:data_number,
      .f = function(id){
        col_data <- data_coldata[[id]]
        row_data <- data_rowdata[[id]]
        if("plate" %in% colnames(col_data)){
          group <- col_data$"plate"
        }else if("group" %in% colnames(col_data)){
          group <- col_data$"group"
        }else{
          NULL
        }
        group <- as.character(group)
        group_combn <- utils::combn(unique(group), 2)
        de_genes <- data_rownames[[id]][which(row_data$"de_gene" == "yes")]
        ## iterate to extract the DEGs by every pair of groups
        sim_DEGs <- list()
        combn_number <- dim(group_combn)[2]
        for(w in 1:combn_number){
          conb1 <- group_combn[1, w]
          conb2 <- group_combn[2, w]
          conb_name <- paste0(group_combn[, w], collapse = "vs")
          fac1 <- row_data[, stringr::str_ends(colnames(row_data), pattern = conb1)]
          fac2 <- row_data[, stringr::str_ends(colnames(row_data), pattern = conb2)]
          index <- fac1 != fac2
          DEGs <- data_rownames[[id]][index]
          sim_DEGs[[conb_name]] <- DEGs
        }
        sim_DEGs
      }
    ) %>% stats::setNames(names(data))
  }

  ## dataframe or matrix
  if(class(data) == "matrix" | class(data) == "data.frame"){
    if(class(data) == "data.frame"){
      data <- list("data" = as.matrix(data))
    }else{
      data <- list("data" = data)
    }
    ## DEGs
    if(!is.null(DEGs)){
      message("The DEGs information is input...")
      if(is.null(group)){
        stop("The group information is necessary when DEGs have already input")
      }
      ## check names of DEGs
      DEGs_names <- names(DEGs)
      group_names <- unique(group)
      for(group_compare in DEGs_names){
        group1_name <- stringr::str_split(group_compare[1],
                                          pattern = "vs",
                                          simplify = TRUE)[1]
        group2_name <- stringr::str_split(group_compare[1],
                                          pattern = "vs",
                                          simplify = TRUE)[2]
        if(!group1_name %in% group_names | !group2_name %in% group_names){
          stop("The names of group in DEGs are not in the group information of cells")
        }
      }
      DEGs <- list("data" = DEGs)
    }
    ## group
    if(!is.null(group)){
      message("The group information is input...")
      group <- list("data" = group)
    }
    ## batch
    if(!is.null(batch)){
      message("The batch information is input...")
      batch <- list("data" = batch)
    }
  }

  ### DEGs evaluation
  DEGs_evaluation <- simutils::calculate_DEGs_properties(
    count_matrix = count_matrices,
    group = group,
    DEGs = DEGs
  )

  return(DEGs_evaluation)
}
