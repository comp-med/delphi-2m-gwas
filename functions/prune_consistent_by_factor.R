prune_consistent_by_factor <- function(dt, factor_col = "embedding", overlap_threshold = 0.2) {
  
  # Step A: Global Priority remains the same
  global_priority <- dt[, .(best_fdr = min(fdr)), by = trait][order(best_fdr)]
  priority_traits <- global_priority$trait
  
  # Step B: Modified Subset Pruning
  prune_subset <- function(sub_df) {
    if (nrow(sub_df) == 0) return(NULL)
    
    # 1. Sort by Global Rank
    sub_df[["global_rank"]] <- match(sub_df[["trait"]], priority_traits)
    sub_df <- sub_df[order(sub_df[["global_rank"]]), ]
    
    sub_df[["collapsed_traits"]] <- ""
    keep_indices <- c(1)
    
    if (nrow(sub_df) > 1) {
      for (i in 2:nrow(sub_df)) {
        current_snps <- sub_df$snp_list[[i]]
        current_name <- sub_df$trait[[i]]
        is_redundant <- FALSE
        
        for (j in keep_indices) {
          previous_snps <- sub_df$snp_list[[j]]
          
          # Intersection count
          inter_count <- length(intersect(current_snps, previous_snps))
          
          if (inter_count == 0) next
          
          # LOGIC FIX: 
          # Calculate how much of the CURRENT trait is covered by the PREVIOUS (higher rank) trait
          overlap_of_current <- inter_count / length(current_snps)
          
          # Also check if the PREVIOUS trait is actually just a tiny subset of the CURRENT one
          # (We don't want a 5-SNP global winner to prune a 100-SNP local winner)
          overlap_of_previous <- inter_count / length(previous_snps)
          
          # Only prune if the higher-ranked trait (j) covers the lower-ranked trait (i) 
          # AND the higher-ranked trait isn't significantly smaller than the current one locally.
          if (overlap_of_current > overlap_threshold && length(previous_snps) >= (length(current_snps) * 0.5)) {
            existing <- sub_df$collapsed_traits[j]
            sub_df$collapsed_traits[j] <- if(existing == "") current_name else paste(existing, current_name, sep = "|")
            is_redundant <- TRUE
            break
          }
        }
        if (!is_redundant) keep_indices <- c(keep_indices, i)
      }
    }
    return(as.data.table(sub_df[keep_indices, ]))
  }
  
  final_pruned <- dt[, prune_subset(as.data.frame(.SD)), by = factor_col]
  if("global_rank" %in% names(final_pruned)) final_pruned[, global_rank := NULL]
  
  return(final_pruned)
}