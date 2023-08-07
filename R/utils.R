#' Prepare the list of DEGs of each group pair
#' @param col_meta cell information
#' @param row_meta gene information
#' @export
prepare_DEGs_list <- function(col_meta, row_meta){
  ### Check groups and DEGs
  if(length(unique(col_meta$group)) >= 2 & "de_gene" %in% colnames(row_meta) |
     length(unique(col_meta$plate)) >= 2 & "de_gene" %in% colnames(row_meta) |
     length(unique(col_meta$group)) >= 2 & "DEstatus" %in% colnames(row_meta)){

    ### group information
    group <- as.character(col_meta$group)

    ### group combination
    group_combn <- utils::combn(unique(group), 2)

    ### special situation in scDD
    if("DEstatus" %in% colnames(row_meta)){
      de_genes_index <- stringr::str_starts(row_meta$"DEstatus", pattern = "^D")
      de_genes <- row_meta$gene_name[de_genes_index]
    }else{
      de_genes <- row_meta$gene_name[which(row_meta$"de_gene" == "yes")]
    }

    ### the list of DEGs
    sim_DEGs <- list()
    for(conb in 1:ncol(group_combn)){
      conb1 <- group_combn[1, conb]
      conb2 <- group_combn[2, conb]
      conb_name <- paste0(group_combn[, conb], collapse = "vs")
      message(conb_name)
      ### Splat, SCRIP, Lun, ESCO (every pair of groups has its own DEGs)
      if(any(stringr::str_detect(colnames(row_meta), "DEF"))){
        fac1 <- row_meta[, stringr::str_ends(colnames(row_meta), pattern = conb1)]
        fac2 <- row_meta[, stringr::str_ends(colnames(row_meta), pattern = conb2)]
        index <- fac1 != fac2
        DEGs <- row_meta$gene_name[index]
      }else if("DEstatus" %in% colnames(row_meta)){
        ### scDD
        DEGs <- de_genes
      }else{
        ### powsimR, muscat, scDesign, SPARSim, SPsimSeq, Lun2 (every pair of groups dose not have its own DEGs)
        if(ncol(group_combn) == 1){
          DEGs <- de_genes
        }else{
          DEGs <- de_genes
        }
      }
      sim_DEGs[[conb_name]] <- DEGs
    }
    return(sim_DEGs)
  }else{
    stop("There are no DEGs in the metadata of simulated genes or no cell groups in the metadata of simulated cells")
  }
}
