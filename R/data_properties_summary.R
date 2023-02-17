.transformation <- function(x){
  if(is.matrix(x)){
    sort(x[upper.tri(x)])
  }else if(length(x) == 1){
    x
  }else{
    sort(x)
  }
}


#' Summary Data Properties
#'
#' @param ref_data_cell_properties Cell properties of reference data
#' @param sim_data_cell_properties Cell properties of simulated data
#' @param ref_data_gene_properties Gene properties of reference data
#' @param sim_data_gene_properties Gene properties of simulated data
#' @param ncore The number of CPU cores to use.
#'
#'
#' @return A list
#' @export
data_properties_summary <- function(
  ref_data_cell_properties,
  sim_data_cell_properties,
  ref_data_gene_properties,
  sim_data_gene_properties,
  ncore = 1
){
  ###------------------------------------------------------------------------###
  ###                           Check
  ###------------------------------------------------------------------------###
  cell_properties <- names(ref_data_cell_properties)
  for(i in cell_properties){
    ref_length <- length(ref_data_cell_properties[[i]])
    sim_length <- length(sim_data_cell_properties[[i]])
    if(ref_length != sim_length){
      stop(glue::glue("The length of property {i} is not equal in reference data properties and simulation data properties."))
    }
  }
  gene_properties <- names(ref_data_gene_properties)
  for(i in gene_properties){
    ref_length <- length(ref_data_gene_properties[[i]])
    sim_length <- length(sim_data_gene_properties[[i]])
    if(ref_length != sim_length){
      stop(glue::glue("The length of property {i} is not equal in reference data properties and simulation data properties."))
    }
  }

  ref_data_cell_properties <- lapply(ref_data_cell_properties, .transformation)
  sim_data_cell_properties <- lapply(sim_data_cell_properties, .transformation)
  ref_data_gene_properties <- lapply(ref_data_gene_properties, .transformation)
  sim_data_gene_properties <- lapply(sim_data_gene_properties, .transformation)


  ###------------------------------------------------------------------------###
  ###                           Cell Level
  ###------------------------------------------------------------------------###

  ## MAD
  message("Cell properties...")
  message("1-MAD")
  MAD_library <- median(abs(ref_data_cell_properties[[1]] - sim_data_cell_properties[[1]]))
  MAD_cellzero <- median(abs(ref_data_cell_properties[[2]] - sim_data_cell_properties[[2]]))
  MAD_cellcor <- median(abs(ref_data_cell_properties[[3]] - sim_data_cell_properties[[3]]))
  if(is.na(sim_data_cell_properties[[4]][1])){
    MAD_TMM <- NA
  }else{
    MAD_TMM <- median(abs(ref_data_cell_properties[[4]] - sim_data_cell_properties[[4]]))
  }
  MAD_elibrary <- median(abs(ref_data_cell_properties[[5]] - sim_data_cell_properties[[5]]))
  MAD_outcell <- median(abs(ref_data_cell_properties[[6]] - sim_data_cell_properties[[6]]))

  ## KS Distance
  if(!requireNamespace("provenance")){
    message("Installing provenance package...")
    install.packages("provenance")
  }
  message("2-KS")
  KS_library <- provenance::KS.diss(ref_data_cell_properties[[1]], sim_data_cell_properties[[1]])
  KS_cellzero <- provenance::KS.diss(ref_data_cell_properties[[2]], sim_data_cell_properties[[2]])
  KS_cellcor <- provenance::KS.diss(ref_data_cell_properties[[3]], sim_data_cell_properties[[3]])
  if(is.na(sim_data_cell_properties[[4]][1])){
    KS_TMM <- NA
  }else{
    KS_TMM <- provenance::KS.diss(ref_data_cell_properties[[4]], sim_data_cell_properties[[4]])
  }
  KS_elibrary <- provenance::KS.diss(ref_data_cell_properties[[5]], sim_data_cell_properties[[5]])

  ## MAE
  if(!requireNamespace("MLmetrics")){
    message("Installing MLmetrics package...")
    install.packages("MLmetrics")
  }
  message("3-MAE")
  MAE_library <- MLmetrics::MAE(sim_data_cell_properties[[1]], ref_data_cell_properties[[1]])
  MAE_cellzero <- MLmetrics::MAE(sim_data_cell_properties[[2]], ref_data_cell_properties[[2]])
  MAE_cellcor <- MLmetrics::MAE(sim_data_cell_properties[[3]], ref_data_cell_properties[[3]])
  if(is.na(sim_data_cell_properties[[4]][1])){
    MAE_TMM <- NA
  }else{
    MAE_TMM <- MLmetrics::MAE(sim_data_cell_properties[[4]], ref_data_cell_properties[[4]])
  }
  MAE_elibrary <- MLmetrics::MAE(sim_data_cell_properties[[5]], ref_data_cell_properties[[5]])
  MAE_outcell <- MLmetrics::MAE(sim_data_cell_properties[[6]], ref_data_cell_properties[[6]])

  ## RMSE
  message("4-RMSE")
  RMSE_library <- MLmetrics::RMSE(sim_data_cell_properties[[1]], ref_data_cell_properties[[1]])
  RMSE_cellzero <- MLmetrics::RMSE(sim_data_cell_properties[[2]], ref_data_cell_properties[[2]])
  RMSE_cellcor <- MLmetrics::RMSE(sim_data_cell_properties[[3]], ref_data_cell_properties[[3]])
  if(is.na(sim_data_cell_properties[[4]][1])){
    RMSE_TMM <- NA
  }else{
    RMSE_TMM <- MLmetrics::RMSE(sim_data_cell_properties[[4]], ref_data_cell_properties[[4]])
  }
  RMSE_elibrary <- MLmetrics::RMSE(sim_data_cell_properties[[5]], ref_data_cell_properties[[5]])
  RMSE_outcell <- MLmetrics::RMSE(sim_data_cell_properties[[6]], ref_data_cell_properties[[6]])

  ## OV
  if(!requireNamespace("overlapping")){
    message("Installing overlapping package...")
    install.packages("overlapping")
  }
  message("5-OV")
  OV_library <- as.numeric(overlapping::overlap(list(x = ref_data_cell_properties[[1]],
                                                     y = sim_data_cell_properties[[1]]))[["OV"]])
  OV_cellzero <- as.numeric(overlapping::overlap(list(x = ref_data_cell_properties[[2]],
                                                      y = sim_data_cell_properties[[2]]))[["OV"]])
  OV_cellcor <- as.numeric(overlapping::overlap(list(x = ref_data_cell_properties[[3]],
                                                     y = sim_data_cell_properties[[3]]))[["OV"]])
  if(is.na(sim_data_cell_properties[[4]][1])){
    OV_TMM <- NA
  }else{
    OV_TMM <- as.numeric(overlapping::overlap(list(x = ref_data_cell_properties[[4]],
                                                   y = sim_data_cell_properties[[4]]))[["OV"]])
  }
  OV_elibrary <- as.numeric(overlapping::overlap(list(x = ref_data_cell_properties[[5]],
                                                      y = sim_data_cell_properties[[5]]))[["OV"]])
  ## bhattacharyya distance
  if(!requireNamespace("philentropy")){
    message("Installing philentropy package...")
    install.packages("philentropy")
  }
  message("6-bhattacharyya distance")
  BH_library <- philentropy::distance(rbind(ref_data_cell_properties$library_size/sum(ref_data_cell_properties$library_size),
                                            sim_data_cell_properties$library_size/sum(sim_data_cell_properties$library_size)),
                                      method = "bhattacharyya")
  BH_cellzero <- philentropy::distance(rbind(ref_data_cell_properties$zero_fraction_cell/sum(ref_data_cell_properties$zero_fraction_cell),
                                             sim_data_cell_properties$zero_fraction_cell/sum(sim_data_cell_properties$zero_fraction_cell)),
                                       method = "bhattacharyya")
  BH_cellcor <- philentropy::distance(rbind(ref_data_cell_properties$cell_cor/sum(ref_data_cell_properties$cell_cor),
                                            sim_data_cell_properties$cell_cor/sum(sim_data_cell_properties$cell_cor)),
                                      method = "bhattacharyya")
  if(is.na(sim_data_cell_properties$TMM_factor)){
    BH_TMM <- NA
  }else{
    BH_TMM <- philentropy::distance(rbind(ref_data_cell_properties$TMM_factor/sum(ref_data_cell_properties$TMM_factor),
                                          sim_data_cell_properties$TMM_factor/sum(sim_data_cell_properties$TMM_factor)),
                                    method = "bhattacharyya")
  }
  BH_elibrary <- philentropy::distance(rbind(ref_data_cell_properties$effective_library_size/sum(ref_data_cell_properties$effective_library_size),
                                             sim_data_cell_properties$effective_library_size/sum(sim_data_cell_properties$effective_library_size)),
                                       method = "bhattacharyya")


  ## Bivariate
  if(!requireNamespace("RcppParallel")){
    message("Installing RcppParallel package...")
    install.packages("RcppParallel")
  }
  if(!requireNamespace("fasano.franceschini.test")){
    message("Installing fasano.franceschini.test package...")
    install.packages("fasano.franceschini.test")
  }
  if(!requireNamespace("ks")){
    message("Installing ks package...")
    install.packages("ks")
  }
  message("7-library size vs zero fraction of cells")
  libraryvscellzero_ref <- matrix(c(log10(ref_data_cell_properties$library_size)+1,
                                    ref_data_cell_properties$zero_fraction_cell), ncol = 2)
  libraryvscellzero_sim <- matrix(c(log10(sim_data_cell_properties$library_size)+1,
                                    sim_data_cell_properties$zero_fraction_cell), ncol = 2)
  libraryvscellzero <- fasano.franceschini.test::fasano.franceschini.test(
    libraryvscellzero_ref,
    libraryvscellzero_sim,
    threads = ncore)
  libraryvscellzero <- mean(libraryvscellzero$estimate)

  KDE_libraryvscellzero <- ks::kde.test(libraryvscellzero_ref, libraryvscellzero_sim)$zstat

  ###------------------------------------------------------------------------###
  ###                           Gene Level
  ###------------------------------------------------------------------------###

  # cv length
  add_number <- abs(length(ref_data_gene_properties$cv) - length(sim_data_gene_properties$cv))
  if(length(ref_data_gene_properties$cv) >= length(sim_data_gene_properties$cv)){
    sample_number <- sample(sim_data_gene_properties$cv, add_number, replace = TRUE)
    sim_data_gene_properties$cv <- c(sim_data_gene_properties$cv, sample_number)
    sim_data_gene_properties$cv <- sort(sim_data_gene_properties$cv)
  }else{
    sample_number <- sample(ref_data_gene_properties$cv, add_number, replace = TRUE)
    ref_data_gene_properties$cv <- c(ref_data_gene_properties$cv, sample_number)
    ref_data_gene_properties$cv <- sort(ref_data_gene_properties$cv)
  }


  ## MAD
  message("Gene properties...")
  message("1-MAD")
  MAD_mean <- median(abs(ref_data_gene_properties[[1]] - sim_data_gene_properties[[1]]))
  MAD_sd <- median(abs(ref_data_gene_properties[[2]] - sim_data_gene_properties[[2]]))
  MAD_cv <- median(abs(ref_data_gene_properties[[3]] - sim_data_gene_properties[[3]]))
  MAD_genezero <- median(abs(ref_data_gene_properties[[4]] - sim_data_gene_properties[[4]]))
  MAD_dispersion <- median(abs(ref_data_gene_properties[[5]] - sim_data_gene_properties[[5]]))
  MAD_outgene <- median(abs(ref_data_gene_properties[[6]] - sim_data_gene_properties[[6]]))

  ## KS Distance
  if(!requireNamespace("provenance")){
    message("Installing provenance package...")
    install.packages("provenance")
  }
  message("2-KS")
  KS_mean <- provenance::KS.diss(ref_data_gene_properties[[1]], sim_data_gene_properties[[1]])
  KS_sd <- provenance::KS.diss(ref_data_gene_properties[[2]], sim_data_gene_properties[[2]])
  KS_cv <- provenance::KS.diss(ref_data_gene_properties[[3]], sim_data_gene_properties[[3]])
  KS_genezero <- provenance::KS.diss(ref_data_gene_properties[[4]], sim_data_gene_properties[[4]])
  KS_dispersion <- provenance::KS.diss(ref_data_gene_properties[[5]], sim_data_gene_properties[[5]])

  ## MAE
  if(!requireNamespace("MLmetrics")){
    message("Installing MLmetrics package...")
    install.packages("MLmetrics")
  }
  message("3-MAE")
  MAE_mean <- MLmetrics::MAE(sim_data_gene_properties[[1]], ref_data_gene_properties[[1]])
  MAE_sd <- MLmetrics::MAE(sim_data_gene_properties[[2]], ref_data_gene_properties[[2]])
  MAE_cv <- MLmetrics::MAE(sim_data_gene_properties[[3]], ref_data_gene_properties[[3]])
  MAE_genezero <- MLmetrics::MAE(sim_data_gene_properties[[4]], ref_data_gene_properties[[4]])
  MAE_dispersion <- MLmetrics::MAE(sim_data_gene_properties[[5]], ref_data_gene_properties[[5]])
  MAE_outgene <- MLmetrics::MAE(sim_data_gene_properties[[6]], ref_data_gene_properties[[6]])

  ## RMSE
  message("4-RMSE")
  RMSE_mean <- MLmetrics::RMSE(sim_data_gene_properties[[1]], ref_data_gene_properties[[1]])
  RMSE_sd <- MLmetrics::RMSE(sim_data_gene_properties[[2]], ref_data_gene_properties[[2]])
  RMSE_cv <- MLmetrics::RMSE(sim_data_gene_properties[[3]], ref_data_gene_properties[[3]])
  RMSE_genezero <- MLmetrics::RMSE(sim_data_gene_properties[[4]], ref_data_gene_properties[[4]])
  RMSE_dispersion <- MLmetrics::RMSE(sim_data_gene_properties[[5]], ref_data_gene_properties[[5]])
  RMSE_outgene <- MLmetrics::RMSE(sim_data_gene_properties[[6]], ref_data_gene_properties[[6]])

  ## OV
  message("5-OV")
  OV_mean <- as.numeric(overlapping::overlap(list(x = sim_data_gene_properties[[1]],
                                                  y = ref_data_gene_properties[[1]]))[["OV"]])
  OV_sd <- as.numeric(overlapping::overlap(list(x = sim_data_gene_properties[[2]],
                                                y = ref_data_gene_properties[[2]]))[["OV"]])
  OV_cv <- as.numeric(overlapping::overlap(list(x = sim_data_gene_properties[[3]],
                                                y = ref_data_gene_properties[[3]]))[["OV"]])
  OV_genezero <- as.numeric(overlapping::overlap(list(x = sim_data_gene_properties[[4]],
                                                      y = ref_data_gene_properties[[4]]))[["OV"]])
  OV_dispersion <- as.numeric(overlapping::overlap(list(x = sim_data_gene_properties[[5]],
                                                        y = ref_data_gene_properties[[5]]))[["OV"]])

  ## bhattacharyya distance
  message("6-bhattacharyya distance")
  BH_mean <- philentropy::distance(rbind(ref_data_gene_properties$mean_expression/sum(ref_data_gene_properties$mean_expression),
                                         sim_data_gene_properties$mean_expression/sum(sim_data_gene_properties$mean_expression)),
                                   method = "bhattacharyya")
  BH_sd <- philentropy::distance(rbind(ref_data_gene_properties$sd/sum(ref_data_gene_properties$sd),
                                       sim_data_gene_properties$sd/sum(sim_data_gene_properties$sd)),
                                 method = "bhattacharyya")
  BH_cv <- philentropy::distance(rbind(ref_data_gene_properties$cv/sum(ref_data_gene_properties$cv),
                                       sim_data_gene_properties$cv/sum(sim_data_gene_properties$cv)),
                                 method = "bhattacharyya")
  BH_genezero <- philentropy::distance(rbind(ref_data_gene_properties$zero_fraction_gene/sum(ref_data_gene_properties$zero_fraction_gene),
                                             sim_data_gene_properties$zero_fraction_gene/sum(sim_data_gene_properties$zero_fraction_gene)),
                                       method = "bhattacharyya")
  BH_dispersion <- philentropy::distance(rbind(ref_data_gene_properties$dispersion/sum(ref_data_gene_properties$dispersion),
                                               sim_data_gene_properties$dispersion/sum(sim_data_gene_properties$dispersion)),
                                         method = "bhattacharyya")


  ## Bivariate
  ### mean vs sd
  message("7-mean expression vs sd")
  meanvssd_ref <- matrix(c(ref_data_gene_properties$mean_expression,
                           ref_data_gene_properties$sd), ncol = 2)
  meanvssd_sim <- matrix(c(sim_data_gene_properties$mean_expression,
                           sim_data_gene_properties$sd), ncol = 2)

  meanvssd <- fasano.franceschini.test::fasano.franceschini.test(
    meanvssd_ref,
    meanvssd_sim,
    threads = ncore)
  meanvssd <- mean(meanvssd$estimate)

  KDE_meanvssd <- ks::kde.test(meanvssd_ref, meanvssd_sim)$zstat


  ### mean vs zero
  message("8-mean expression vs zero fraction of genes")
  meanvszero_ref <- matrix(c(ref_data_gene_properties$mean_expression,
                             ref_data_gene_properties$zero_fraction_gene), ncol = 2)
  meanvszero_sim <- matrix(c(sim_data_gene_properties$mean_expression,
                             sim_data_gene_properties$zero_fraction_gene), ncol = 2)

  meanvszero <- fasano.franceschini.test::fasano.franceschini.test(
    meanvszero_ref,
    meanvszero_sim,
    threads = ncore)
  meanvszero <- mean(meanvszero$estimate)

  KDE_meanvszero <- ks::kde.test(meanvszero_ref, meanvszero_sim)$zstat


  ### mean vs dispersion
  message("9-mean expression vs dispersion")
  meanvsdispersion_ref <- matrix(c(ref_data_gene_properties$mean_expression,
                                   ref_data_gene_properties$dispersion), ncol = 2)
  meanvsdispersion_sim <- matrix(c(sim_data_gene_properties$mean_expression,
                                   sim_data_gene_properties$dispersion), ncol = 2)

  meanvsdispersion <- fasano.franceschini.test::fasano.franceschini.test(
    meanvsdispersion_ref,
    meanvsdispersion_sim,
    threads = ncore)
  meanvsdispersion <- mean(meanvsdispersion$estimate)

  KDE_meanvsdispersion <- ks::kde.test(meanvsdispersion_ref, meanvsdispersion_sim)$zstat


  return(list(MAD_library = MAD_library,
              `MAD_zero fraction of cells` = MAD_cellzero,
              `MAD_cell correlation` = MAD_cellcor,
              MAD_TMM = MAD_TMM,
              `MAD_effective library` = MAD_elibrary,
              `MAD_cell outlier` = MAD_outcell,
              KS_library = KS_library,
              `KS_zero fraction of cells` = KS_cellzero,
              `KS_cell correlation` = KS_cellcor,
              KS_TMM = KS_TMM,
              `KS_effective library` = KS_elibrary,
              MAE_library = MAE_library,
              `MAE_zero fraction of cells` = MAE_cellzero,
              `MAE_cell correlation` = MAE_cellcor,
              MAE_TMM = MAE_TMM,
              `MAE_effective library` = MAE_elibrary,
              `MAE_cell outlier` = MAE_outcell,
              RMSE_library = RMSE_library,
              `RMSE_zero fraction of cells` = RMSE_cellzero,
              `RMSE_cell correlation` = RMSE_cellcor,
              RMSE_TMM = RMSE_TMM,
              `RMSE_effective library` = RMSE_elibrary,
              `RMSE_cell outlier` = RMSE_outcell,
              OV_library = OV_library,
              `OV_zero fraction of cells` = OV_cellzero,
              `OV_cell correlation` = OV_cellcor,
              OV_TMM = OV_TMM,
              `OV_effective library` = OV_elibrary,
              `bhattacharyya distance_library` = unname(BH_library),
              `bhattacharyya distance_zero fraction of cells` = unname(BH_cellzero),
              `bhattacharyya distance_cell correlation` = unname(BH_cellcor),
              `bhattacharyya distance_TMM` = unname(BH_TMM),
              `bhattacharyya distance_effective library` = unname(BH_elibrary),
              `KS library size vs zero fraction of cells` = libraryvscellzero,
               `KDE library size vs zero fraction of cells` = KDE_libraryvscellzero,
              `MAD_mean expression` = MAD_mean,
              MAD_sd = MAD_sd,
              MAD_cv = MAD_cv,
              `MAD_zero fraction of genes` = MAD_genezero,
              MAD_dispersion = MAD_dispersion,
              `MAD_gene outlier` = MAD_dispersion,
              `KS_mean expression` = KS_mean,
              KS_sd = KS_sd,
              KS_cv = KS_cv,
              `KS_zero fraction of genes` = KS_genezero,
              KS_dispersion = KS_dispersion,
              `MAE_mean expression` = MAE_mean,
              MAE_sd = MAE_sd,
              MAE_cv = MAE_cv,
              `MAE_zero fraction of genes` = MAE_genezero,
              MAE_dispersion = MAE_dispersion,
              `MAE_gene outlier` = MAE_dispersion,
              `RMSE_mean expression` = RMSE_mean,
              RMSE_sd = RMSE_sd,
              RMSE_cv = RMSE_cv,
              `RMSE_zero fraction of genes` = RMSE_genezero,
              RMSE_dispersion = RMSE_dispersion,
              `RMSE_gene outlier` = RMSE_dispersion,
              `OV_mean expression` = OV_mean,
              OV_sd = OV_sd,
              OV_cv = OV_cv,
              `OV_zero fraction of genes` = OV_genezero,
              OV_dispersion = OV_dispersion,
              `bhattacharyya distance_mean expression` = unname(BH_mean),
              `bhattacharyya distance_sd` = unname(BH_sd),
              `bhattacharyya distance_cv` = unname(BH_cv),
              `bhattacharyya distance_zero fraction of genes` = unname(BH_genezero),
              `bhattacharyya distance_dispersion` = unname(BH_dispersion),
              `KS mean expression vs sd` = meanvssd,
              `KS mean expression vs zero fraction of genes` = meanvszero,
              `KS mean expression vs dispersion` = meanvsdispersion,
              `KDE mean expression vs sd` = KDE_meanvssd,
              `KDE mean expression vs zero fraction of genes` = KDE_meanvszero,
              `KDE mean expression vs dispersion` = KDE_meanvsdispersion))

}


# data <- readRDS("/Users/duohongrui/Desktop/preprocessed_data/data1_GSE54006.rds")
# ref_data_cell_properties <- simutils::cell_properties(data$data, verbose = TRUE)
# ref_data_gene_properties <- simutils::gene_properties(data$data, verbose = TRUE)
#
# sim <- readRDS("/Volumes/Elements/06-personal information/Evaluation of Sim Tools/Evaluation/splatter_package_methods/simulation_data/dropsim_data1_GSE54006.rds")
# sim_data <- sim$sim_data$count_data
# sim_data_cell_properties <- simutils::cell_properties(sim_data, verbose = TRUE)
# sim_data_gene_properties <- simutils::gene_properties(sim_data, verbose = TRUE)
