###########################################################
## function to prune enrichment results for one trait
## generated with Gemini 3

# 3. Clustered Pruning Function
prune_with_tracking <- function(dt, overlap_threshold = 0.2) {
  if (nrow(dt) == 0) return(dt)
  
  # Initialize a list to store collapsed trait names for each row
  dt[, collapsed_traits := ""]
  keep_indices <- c(1)
  
  for (i in 2:nrow(dt)) {
    current_snps <- dt$snp_list[[i]]
    current_name <- dt$trait[i]
    is_redundant <- FALSE
    
    for (j in keep_indices) {
      previous_snps <- dt$snp_list[[j]]
      
      # Jaccard Calculation
      inter_count <- length(intersect(current_snps, previous_snps))
      union_count <- length(union(current_snps, previous_snps))
      jaccard <- if(union_count == 0) 0 else inter_count / union_count
      
      if (!is.na(jaccard) && jaccard > overlap_threshold) {
        # Mark as redundant and add this trait name to the 'Lead' trait's list
        existing <- dt$collapsed_traits[j]
        dt$collapsed_traits[j] <- ifelse(existing == "", current_name, paste(existing, current_name, sep = "|"))
        is_redundant <- TRUE
        break
      }
    }
    
    if (!is_redundant) {
      keep_indices <- c(keep_indices, i)
    }
  }
  
  return(dt[keep_indices])
}
