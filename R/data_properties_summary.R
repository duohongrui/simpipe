.transformation <- function(x){
  if(is.matrix(x)){
    sort(x[upper.tri(x)])
  }else if(length(x) == 1){
    x
  }else{
    sort(x)
  }
}


data_properties_summary <- function(
  ref_data_cell_properties,
  sim_data_cell_properties,
  ref_data_gene_properties,
  sim_data_gene_properties
){
  ref_data_cell_properties <- lapply(ref_data_cell_properties, .transformation)
  sim_data_cell_properties <- lapply(sim_data_cell_properties, .transformation)
  ref_data_gene_properties <- lapply(ref_data_gene_properties, .transformation)
  sim_data_gene_properties <- lapply(sim_data_gene_properties, .transformation)

  ###------------------------------------------------------------------------###
  ###                           Cell Level
  ###------------------------------------------------------------------------###

  ## MAD
  MAD_library <- median(abs(ref_data_cell_properties[[1]] - sim_data_cell_properties[[1]]))
  MAD_cellzero <- median(abs(ref_data_cell_properties[[2]] - sim_data_cell_properties[[2]]))
  MAD_cellcor <- median(abs(ref_data_cell_properties[[3]] - sim_data_cell_properties[[3]]))
  MAD_TMM <- median(abs(ref_data_cell_properties[[4]] - sim_data_cell_properties[[4]]))
  MAD_elibrary <- median(abs(ref_data_cell_properties[[5]] - sim_data_cell_properties[[5]]))
  MAD_outcell <- median(abs(ref_data_cell_properties[[6]] - sim_data_cell_properties[[6]]))

  ## KS Distance
  if(!requireNamespace("provenance")){
    message("Installing provenance package...")
    install.packages("provenance")
  }
  KS_library <- provenance::KS.diss(ref_data_cell_properties[[1]], sim_data_cell_properties[[1]])
  KS_cellzero <- provenance::KS.diss(ref_data_cell_properties[[2]], sim_data_cell_properties[[2]])
  KS_cellcor <- provenance::KS.diss(ref_data_cell_properties[[3]], sim_data_cell_properties[[3]])
  KS_TMM <- provenance::KS.diss(ref_data_cell_properties[[4]], sim_data_cell_properties[[4]])
  KS_elibrary <- provenance::KS.diss(ref_data_cell_properties[[5]], sim_data_cell_properties[[5]])

  ## Correlation Distance
  cor_library <- cor(ref_data_cell_properties[[1]], sim_data_cell_properties[[1]], method = "pearson")
  cor_cellzero <- cor(ref_data_cell_properties[[2]], sim_data_cell_properties[[2]], method = "pearson")
  cor_cellcor <- cor(ref_data_cell_properties[[3]], sim_data_cell_properties[[3]], method = "pearson")
  cor_TMM <- cor(ref_data_cell_properties[[4]], sim_data_cell_properties[[4]], method = "pearson")
  cor_elibrary <- cor(ref_data_cell_properties[[5]], sim_data_cell_properties[[5]], method = "pearson")

  ## MAE
  if(!requireNamespace("MLmetrics")){
    message("Installing MLmetrics package...")
    install.packages("MLmetrics")
  }
  MAE_library <- MLmetrics::MAE(sim_data_cell_properties[[1]], ref_data_cell_properties[[1]])
  MAE_cellzero <- MLmetrics::MAE(sim_data_cell_properties[[2]], ref_data_cell_properties[[2]])
  MAE_cellcor <- MLmetrics::MAE(sim_data_cell_properties[[3]], ref_data_cell_properties[[3]])
  MAE_TMM <- MLmetrics::MAE(sim_data_cell_properties[[4]], ref_data_cell_properties[[4]])
  MAE_elibrary <- MLmetrics::MAE(sim_data_cell_properties[[5]], ref_data_cell_properties[[5]])
  MAE_outcell <- MLmetrics::MAE(sim_data_cell_properties[[6]], ref_data_cell_properties[[6]])

  ## RMSE
  RMSE_library <- MLmetrics::RMSE(sim_data_cell_properties[[1]], ref_data_cell_properties[[1]])
  RMSE_cellzero <- MLmetrics::RMSE(sim_data_cell_properties[[2]], ref_data_cell_properties[[2]])
  RMSE_cellcor <- MLmetrics::RMSE(sim_data_cell_properties[[3]], ref_data_cell_properties[[3]])
  RMSE_TMM <- MLmetrics::RMSE(sim_data_cell_properties[[4]], ref_data_cell_properties[[4]])
  RMSE_elibrary <- MLmetrics::RMSE(sim_data_cell_properties[[5]], ref_data_cell_properties[[5]])
  RMSE_outcell <- MLmetrics::RMSE(sim_data_cell_properties[[6]], ref_data_cell_properties[[6]])

  ## Bivariate
  if(!requireNamespace("RcppParallel")){
    message("Installing RcppParallel package...")
    install.packages("RcppParallel")
  }
  if(!requireNamespace("fasano.franceschini.test")){
    message("Installing fasano.franceschini.test package...")
    install.packages("fasano.franceschini.test")
  }
  libraryvscellzero_ref <- matrix(c(log10(ref_data_cell_properties$library_size)+1,
                                    ref_data_cell_properties$zero_fraction_cell), ncol = 2)
  libraryvscellzero_sim <- matrix(c(log10(sim_data_cell_properties$library_size)+1,
                                    sim_data_cell_properties$zero_fraction_cell), ncol = 2)
  libraryvscellzero <- fasano.franceschini.test::fasano.franceschini.test(
    libraryvscellzero_ref,
    libraryvscellzero_sim,
    threads = RcppParallel::defaultNumThreads(),
    method = "o")
  libraryvscellzero <- mean(libraryvscellzero$estimate)

  ###------------------------------------------------------------------------###
  ###                           Gene Level
  ###------------------------------------------------------------------------###

  # cv length
  cv_length <- min(length(ref_data_gene_properties$cv), length(sim_data_gene_properties$cv))
  ref_data_gene_properties$cv <- ref_data_gene_properties$cv[1:cv_length]
  sim_data_gene_properties$cv <- sim_data_gene_properties$cv[1:cv_length]
  # gene cor length
  genecor_length <- min(length(ref_data_gene_properties$gene_cor), length(sim_data_gene_properties$gene_cor))
  ref_data_gene_properties$gene_cor <- ref_data_gene_properties$gene_cor[1:genecor_length]
  sim_data_gene_properties$gene_cor <- sim_data_gene_properties$gene_cor[1:genecor_length]

  ## MAD
  MAD_mean <- median(abs(ref_data_gene_properties[[1]] - sim_data_gene_properties[[1]]))
  MAD_sd <- median(abs(ref_data_gene_properties[[2]] - sim_data_gene_properties[[2]]))
  MAD_cv <- median(abs(ref_data_gene_properties[[3]] - sim_data_gene_properties[[3]]))
  MAD_genecor <- median(abs(ref_data_gene_properties[[4]] - sim_data_gene_properties[[4]]))
  MAD_genezero <- median(abs(ref_data_gene_properties[[5]] - sim_data_gene_properties[[5]]))
  MAD_dispersion <- median(abs(ref_data_gene_properties[[6]] - sim_data_gene_properties[[6]]))
  MAD_outgene <- median(abs(ref_data_gene_properties[[7]] - sim_data_gene_properties[[7]]))

  ## KS Distance
  if(!requireNamespace("provenance")){
    message("Installing provenance package...")
    install.packages("provenance")
  }
  KS_mean <- provenance::KS.diss(ref_data_gene_properties[[1]], sim_data_gene_properties[[1]])
  KS_sd <- provenance::KS.diss(ref_data_gene_properties[[2]], sim_data_gene_properties[[2]])
  KS_cv <- provenance::KS.diss(ref_data_gene_properties[[3]], sim_data_gene_properties[[3]])
  KS_genecor <- provenance::KS.diss(ref_data_gene_properties[[4]], sim_data_gene_properties[[4]])
  KS_genezero <- provenance::KS.diss(ref_data_gene_properties[[5]], sim_data_gene_properties[[5]])
  KS_dispersion <- provenance::KS.diss(ref_data_gene_properties[[6]], sim_data_gene_properties[[6]])

  ## Correlation Distance
  cor_mean <- cor(ref_data_gene_properties[[1]], sim_data_gene_properties[[1]], method = "pearson")
  cor_sd <- cor(ref_data_gene_properties[[2]], sim_data_gene_properties[[2]], method = "pearson")
  cor_cv <- cor(ref_data_gene_properties[[3]], sim_data_gene_properties[[3]], method = "pearson")
  cor_genecor <- cor(ref_data_gene_properties[[4]], sim_data_gene_properties[[4]], method = "pearson")
  cor_genezero <- cor(ref_data_gene_properties[[5]], sim_data_gene_properties[[5]], method = "pearson")
  cor_dispersion <- cor(ref_data_gene_properties[[6]], sim_data_gene_properties[[6]], method = "pearson")

  ## MAE
  if(!requireNamespace("MLmetrics")){
    message("Installing MLmetrics package...")
    install.packages("MLmetrics")
  }
  MAE_mean <- MLmetrics::MAE(sim_data_gene_properties[[1]], ref_data_gene_properties[[1]])
  MAE_sd <- MLmetrics::MAE(sim_data_gene_properties[[2]], ref_data_gene_properties[[2]])
  MAE_cv <- MLmetrics::MAE(sim_data_gene_properties[[3]], ref_data_gene_properties[[3]])
  MAE_genecor <- MLmetrics::MAE(sim_data_gene_properties[[4]], ref_data_gene_properties[[4]])
  MAE_genezero <- MLmetrics::MAE(sim_data_gene_properties[[5]], ref_data_gene_properties[[5]])
  MAE_dispersion <- MLmetrics::MAE(sim_data_gene_properties[[6]], ref_data_gene_properties[[6]])
  MAE_outgene <- MLmetrics::MAE(sim_data_gene_properties[[7]], ref_data_gene_properties[[7]])

  ## RMSE
  RMSE_mean <- MLmetrics::RMSE(sim_data_gene_properties[[1]], ref_data_gene_properties[[1]])
  RMSE_sd <- MLmetrics::RMSE(sim_data_gene_properties[[2]], ref_data_gene_properties[[2]])
  RMSE_cv <- MLmetrics::RMSE(sim_data_gene_properties[[3]], ref_data_gene_properties[[3]])
  RMSE_genecor <- MLmetrics::RMSE(sim_data_gene_properties[[4]], ref_data_gene_properties[[4]])
  RMSE_genezero <- MLmetrics::RMSE(sim_data_gene_properties[[5]], ref_data_gene_properties[[5]])
  RMSE_dispersion <- MLmetrics::RMSE(sim_data_gene_properties[[6]], ref_data_gene_properties[[6]])
  RMSE_outgene <- MLmetrics::RMSE(sim_data_gene_properties[[7]], ref_data_gene_properties[[7]])

  ## Bivariate
  ### mean vs sd
  meanvssd_ref <- matrix(c(ref_data_gene_properties$mean_expression,
                           ref_data_gene_properties$sd), ncol = 2)
  meanvssd_sim <- matrix(c(sim_data_gene_properties$mean_expression,
                           sim_data_gene_properties$sd), ncol = 2)

  meanvssd <- fasano.franceschini.test::fasano.franceschini.test(
    meanvssd_ref,
    meanvssd_sim,
    threads = RcppParallel::defaultNumThreads(),
    method = "o")
  meanvssd <- mean(meanvssd$estimate)


  ### mean vs zero
  meanvszero_ref <- matrix(c(ref_data_gene_properties$mean_expression,
                             ref_data_gene_properties$zero_fraction_gene), ncol = 2)
  meanvszero_sim <- matrix(c(sim_data_gene_properties$mean_expression,
                             sim_data_gene_properties$zero_fraction_gene), ncol = 2)

  meanvszero <- fasano.franceschini.test::fasano.franceschini.test(
    meanvszero_ref,
    meanvszero_sim,
    threads = RcppParallel::defaultNumThreads(),
    method = "o")
  meanvszero <- mean(meanvszero$estimate)


  ### mean vs dispersion
  meanvsdispersion_ref <- matrix(c(ref_data_gene_properties$mean_expression,
                                   ref_data_gene_properties$dispersion), ncol = 2)
  meanvsdispersion_sim <- matrix(c(sim_data_gene_properties$mean_expression,
                                   sim_data_gene_properties$dispersion), ncol = 2)

  meanvsdispersion <- fasano.franceschini.test::fasano.franceschini.test(
    meanvsdispersion_ref,
    meanvsdispersion_sim,
    threads = RcppParallel::defaultNumThreads(),
    method = "o")
  meanvsdispersion <- mean(meanvsdispersion$estimate)


  return(list(MAD_library = MAD_library,
              MAD_zero-fraction-cell = MAD_cellzero,
              MAD_cell-correlation = MAD_cellcor,
              MAD_TMM = MAD_TMM,
              MAD_effective-library = MAD_elibrary,
              MAD_cell-outlier = MAD_outcell,
              KS_library = KS_library,
              KS_zero-fraction-cell = KS_cellzero,
              KS_cell-correlation = KS_cellcor,
              KS_TMM = KS_TMM,
              KS_effective-library = KS_elibrary))

}


data <- readRDS("/Users/duohongrui/Desktop/preprocessed_data/data1_GSE54006.rds")
ref_data_cell_properties <- simutils::cell_properties(data$data, verbose = TRUE)
ref_data_gene_properties <- simutils::gene_properties(data$data, verbose = TRUE)

sim <- readRDS("/Volumes/Elements/06-personal information/Evaluation of Sim Tools/Evaluation/splatter_package_methods/simulation_data/dropsim_data1_GSE54006.rds")
sim_data <- sim$sim_data$count_data
sim_data_cell_properties <- simutils::cell_properties(sim_data, verbose = TRUE)
sim_data_gene_properties <- simutils::gene_properties(sim_data, verbose = TRUE)
