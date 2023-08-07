#' Summarize the Functionality in Group, Batch and DEGs of A Simulation Method
#'
#' @param data The matrix or data frame of gene expression profile or the output by `simpipe::simulate_datasets`.
#' @param group The characters of group names which cells belong to. Needed if data is a matrix or a dataframe and you want to evaluate the groups, DEGs and trajectories.
#' @param batch The characters of batch names which cells belong to. Needed if data is a matrix or a dataframe and you want to evaluate the batches
#' @param k k-nearest neighborhoods of cells.
#' @param DEGs A list of DEGs with the names of `xxxvsxxx`. Note that the names of DEGs are in the rownames of the matrix or the dataframe and the names of `xxx` is in the `group` characters. Needed if data is a matrix or a dataframe and you want to evaluate the DEGs.
#' @param DEA_method The DEA method to get the DEGs. Choices: edgeRQLF, edgeRQLFDetRate, MASTcpmDetRate, MASTtpmDetRate, MASTcpm, MASTtpm, limmatrend, limmavoom, ttest and wilcox. Default is edgeRQLFDetRate.
#' @param model_method The method to establish the model. SVM, Decision tree or RF (Random Forest). Default is SVM.
#' @param ref_data A matrix, data frame or list of real gene expression profiles.
#' @param ref_data_grouping The vector or list of group names which real cells belong to.
#' @param algorithm The algorithm for matching the real cells and the simulated cells. Hungarian (default) or Improved_Hungarian.
#' @param seed Random seed for trajectory inference.
#' @param verbose Whether the process massages are returned.
#' @param threads How many cores used for parallel computation.
#' @importFrom SingleCellExperiment counts colData rowData
#' @importFrom simutils calculate_DEGs_properties calculate_batch_properties calculate_cluster_properties
#' @importFrom methods is
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
  DEA_method = "edgeRQLFDetRate",
  model_method = "SVM",
  ref_data = NULL,
  ref_data_grouping = NULL,
  algorithm = "Hungarian",
  seed = 1,
  verbose = TRUE,
  threads = 1
){
  ### data check and other imformation of cells and genes
  if("simpipe_simulation" %in% class(data)){
    data_number <- length(data)

    ## Check the return format
    if(methods::is(data[[1]][["simulate_result"]], "SingleCellExperiment")){
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
    if(methods::is(data[[1]][["simulate_result"]], "Seurat")){
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
    if(is.list(data[[1]][["simulate_result"]]) &
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
          group_evl <- as.character(col_data$"plate")
        }else if("group" %in% colnames(col_data)){
          group_evl <- as.character(col_data$"group")
        }else{
          group_evl <- 0
        }
        if(length(unique(group_evl)) == 1){
          group_evl <- NULL
        }
        group_evl
      }
    ) %>% stats::setNames(names(data))
    ### Batch information
    batch <- purrr::map(
      .x = 1:data_number,
      .f = function(id){
        col_data <- data_coldata[[id]]
        if("batch" %in% colnames(col_data)){
          batch_evl <- col_data$"batch"
        }else{
          batch_evl <- 0
        }
        if(length(unique(batch_evl)) == 1){
          batch_evl <- NULL
        }
        batch_evl
      }
    ) %>% stats::setNames(names(data))
    ### DEGs information
    DEGs <- purrr::map(
      .x = 1:data_number,
      .f = function(id){
        col_data <- data_coldata[[id]]
        row_data <- data_rowdata[[id]]
        if("plate" %in% colnames(col_data)){
          group_eval <- col_data$"plate"
        }else if("group" %in% colnames(col_data)){
          group_eval <- col_data$"group"
        }else{
          group_eval <- 0
        }
        if(length(unique(group_eval)) > 1){
          group <- as.character(group_eval)
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
        }else{
          sim_DEGs <- NULL
        }
        sim_DEGs
      }
    ) %>% stats::setNames(names(data))
  }

  ## dataframe or matrix
  if(is.matrix(data) | is.data.frame(data)){
    if(is.data.frame(data)){
      count_matrices <- list("data" = as.matrix(data))
    }else{
      count_matrices <- list("data" = data)
    }
    data_number <- length(count_matrices)
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
  if(!is.null(DEGs)){
    cat("-------------------------------------------------\n")
    cat("                  Evaluating DEGs\n")
    cat("-------------------------------------------------\n")
    DEGs_evaluation <- simutils::calculate_DEGs_properties(
      count_matrix = count_matrices,
      group = group,
      DEGs = DEGs,
      DEA_method = DEA_method,
      model_method = model_method,
      verbose = verbose
    )
  }else{
    DEGs_evaluation <- NULL
  }

  ### batch evaluation
  if(!is.null(batch)){
    cat("-------------------------------------------------\n")
    cat("             Evaluating cell batches\n")
    cat("-------------------------------------------------\n")
    batch_evaluation <- purrr::map(
      .x = 1:data_number,
      .f = function(id){
        batch_data <- count_matrices[[id]]
        batch_info <- batch[[id]]
        if(is.null(batch_info)){
          NULL
        }else{
          simutils::calculate_batch_properties(
            data = batch_data,
            batch_info = batch_info,
            k = k,
            verbose = verbose
          )
        }
      }
    ) %>% stats::setNames(names(count_matrices))
  }else{
    batch_evaluation <- NULL
  }

  ### group evaluation
  if(!is.null(group)){
    cat("-------------------------------------------------\n")
    cat("             Evaluating cell groups\n")
    cat("-------------------------------------------------\n")
    group_evaluation <- purrr::map(
      .x = 1:data_number,
      .f = function(id){
        group_data <- count_matrices[[id]]
        group_info <- group[[id]]
        if(is.null(group_info)){
          NULL
        }else{
          simutils::calculate_cluster_properties(
            data = group_data,
            dist = NULL,
            cluster_info = group_info,
            threads = threads,
            verbose = verbose
          )
        }
      }
    ) %>% stats::setNames(names(count_matrices))
  }else{
    group_evaluation <- NULL
  }

  ### trajectory evaluation
  if(!is.null(ref_data)){
    cat("-------------------------------------------------\n")
    cat("             Evaluating trajectories\n")
    cat("-------------------------------------------------\n")

    ## list of reference dataset
    if(is.list(ref_data)){
      ref_count_matrices <- ref_data
      ref_data_number <- length(ref_count_matrices)
      if(is.list(ref_data_grouping)){
        if(length(ref_data) != length(ref_data_grouping)){
          stop("The length of ref_data must equal to that of information of cell groups.")
        }else{
          ref_group <- ref_data_grouping
        }
      }else{
        if(is.vector(ref_data_grouping)){
          ref_group <- purrr::map(ref_data_number, .f = function(x){ref_data_grouping})
        }else{
          ref_group <- ref_data_grouping
        }
      }
    }

    ## dataframe or matrix
    if(is.matrix(ref_data) | is.data.frame(ref_data)){
      if(is.data.frame(ref_data)){
        ref_count_matrices <- list("data" = as.matrix(ref_data))
      }else{
        ref_count_matrices <- list("data" = ref_data)
      }
      ref_data_number <- length(ref_count_matrices)
      ## group
      if(!is.null(ref_data_grouping)){
        message("The group information of reference data is input...")
        ref_group <- list("data" = ref_data_grouping)
      }else{
        ref_group <- ref_data_grouping
      }
    }

    ## length between reference and simulated datasets
    if(length(count_matrices) != length(ref_count_matrices)){
      stop("The number of reference datasets in a list must equal to that of simulated datasets in the list.")
    }

    ## evaluating trajectories
    trajectory_evaluation <- purrr::map(
      .x = 1:ref_data_number,
      .f = function(id){
        ref_traj_data <- ref_count_matrices[[id]]
        ref_traj_info <- ref_group[[id]]
        simutils::calculate_trajectory_properties(
          ref_data = ref_traj_data,
          ref_data_grouping = ref_traj_info,
          sim_data = count_matrices[[id]],
          sim_data_grouping = group[[id]],
          algorithm = algorithm,
          seed = seed,
          verbose = verbose
        )
      }
    ) %>% stats::setNames(names(count_matrices))
  }else{
    trajectory_evaluation <- NULL
  }

  dplyr::lst(group_evaluation,
             batch_evaluation,
             DEGs_evaluation,
             trajectory_evaluation)

}


# a <- powsimR::CELseq2_Gene_UMI_Counts
#
# estimate_output <- simpipe::estimate_parameters(ref_data = as.matrix(a),
#                                        method = "Splat",
#                                        seed = 10,
#                                        verbose = TRUE,
#                                        use_docker = FALSE)
# data <- simpipe::simulate_datasets(method = "Splat",
#                                      parameters = estimate_output,
#                                      seed = 110,
#                                      other_prior = list(batchCells = c(500, 500),
#                                                         group.prob = c(0.3, 0.3, 0.4),
#                                                         nGenes = 3000),
#                                      return_format = "list",
#                                      verbose = TRUE,
#                                      use_docker = FALSE)
# data2 <- data$refdata_Splat_1$simulate_result$count_data
# group <- data$refdata_Splat_1$simulate_result$col_meta$group
# batch <- data$refdata_Splat_1$simulate_result$col_meta$batch
# DEGs <- NULL

