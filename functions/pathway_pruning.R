library(data.table)

#' Prune redundant pathways using a greedy set cover approach.
#'
#' Within each embedding group (`emb`), iteratively selects the pathway that
#' covers the most yet-uncovered genes. A candidate pathway is skipped if its
#' gene set overlaps too heavily (> `jaccard_threshold`) with any already
#' selected pathway.  Selection priority is governed by `rank_by`:
#'   - "new_genes"  : most newly covered genes first (pure set-cover)
#'   - "p_value"    : most significant pathway first (significance-driven)
#'   - "intersection_size" : largest absolute intersection first
#'
#' @param dt              data.table as loaded from the GSEA result file.
#' @param jaccard_threshold Numeric [0, 1]. Maximum allowed Jaccard index
#'                          between a candidate pathway and any already
#'                          selected pathway.  Default 0.5.
#' @param min_new_genes   Integer. Minimum number of *new* (not yet covered)
#'                        genes a pathway must contribute to be selected.
#'                        Default 1.
#' @param rank_by         Character. One of "new_genes", "p_value", or
#'                        "intersection_size". Default "new_genes".
#' @param gene_sep        Character used to split `intersection_genes`.
#'                        Default "|".
#'
#' @return A data.table of pruned pathways (one row per selected pathway per
#'         embedding), with two extra columns:
#'           - `n_new_genes`    : number of newly covered genes this pathway adds
#'           - `cumulative_genes`: all unique genes covered up to and including
#'                                 this pathway (within the embedding)
prune_pathways <- function(dt,
                           jaccard_threshold = 0.5,
                           min_new_genes     = 1,
                           rank_by           = "new_genes",
                           gene_sep          = "|") {
  
  ## ── input checks ──────────────────────────────────────────────────────────
  stopifnot(is.data.table(dt))
  stopifnot(all(c("emb", "native", "name", "p_value",
                  "intersection_genes") %in% names(dt)))
  stopifnot(rank_by %in% c("new_genes", "p_value", "intersection_size"))
  stopifnot(jaccard_threshold >= 0, jaccard_threshold <= 1)
  
  ## ── helper: Jaccard index between two character vectors ───────────────────
  jaccard <- function(a, b) {
    inter <- length(intersect(a, b))
    if (inter == 0L) return(0)
    inter / length(union(a, b))
  }
  
  ## ── parse gene lists once ─────────────────────────────────────────────────
  # strsplit is vectorised; fixed = TRUE is faster than regex
  dt <- copy(dt)
  dt[, gene_list := strsplit(intersection_genes, gene_sep, fixed = TRUE)]
  
  ## ── process each embedding independently ──────────────────────────────────
  results <- lapply(unique(dt$emb), function(emb_id) {
    
    sub <- dt[emb == emb_id]
    
    # Pre-compute intersection_size from gene_list (ground truth from file)
    sub[, n_genes := lengths(gene_list)]
    
    covered_genes  <- character(0)   # union of genes covered so far
    selected_rows  <- integer(0)     # row indices (within sub) chosen
    selected_lists <- list()         # gene vectors of selected pathways
    
    ## greedy loop: keep iterating until no pathway meets the criteria
    repeat {
      remaining <- setdiff(seq_len(nrow(sub)), selected_rows)
      if (length(remaining) == 0L) break
      
      candidates <- sub[remaining]
      
      ## count new genes each candidate would add
      candidates[, n_new_genes := vapply(
        gene_list,
        function(gl) length(setdiff(gl, covered_genes)),
        integer(1L)
      )]
      
      ## filter by minimum new-gene contribution
      candidates <- candidates[n_new_genes >= min_new_genes]
      if (nrow(candidates) == 0L) break
      
      ## filter by Jaccard overlap against every already-selected pathway
      if (length(selected_lists) > 0L) {
        keep <- vapply(candidates$gene_list, function(gl) {
          all(vapply(selected_lists,
                     function(sel) jaccard(gl, sel),
                     numeric(1L)) <= jaccard_threshold)
        }, logical(1L))
        candidates <- candidates[keep]
        if (nrow(candidates) == 0L) break
      }
      
      ## rank candidates
      if (rank_by == "new_genes") {
        setorder(candidates, -n_new_genes, p_value)
      } else if (rank_by == "p_value") {
        setorder(candidates, p_value, -n_new_genes)
      } else {                          # "intersection_size"
        setorder(candidates, -n_genes, p_value)
      }
      
      ## select top candidate
      winner <- candidates[1L]
      
      ## update state
      winner_genes   <- winner$gene_list[[1L]]
      covered_genes  <- union(covered_genes, winner_genes)
      selected_lists <- c(selected_lists, list(winner_genes))
      
      # find original row index in sub
      winner_idx <- which(sub$native == winner$native &
                            sub$name   == winner$name)[1L]
      selected_rows <- c(selected_rows, winner_idx)
    }
    
    ## collect selected rows and annotate
    if (length(selected_rows) == 0L) return(NULL)
    
    out <- sub[selected_rows]
    
    ## recompute n_new_genes and cumulative gene set in selection order
    cov <- character(0)
    out[, n_new_genes     := 0L]
    out[, cumulative_genes := ""]
    
    for (i in seq_len(nrow(out))) {
      gl  <- out$gene_list[[i]]
      new <- setdiff(gl, cov)
      cov <- union(cov, gl)
      set(out, i, "n_new_genes",     length(new))
      set(out, i, "cumulative_genes", paste(sort(cov), collapse = gene_sep))
    }
    out
  })
  
  ## ── combine and clean up ──────────────────────────────────────────────────
  out_dt <- rbindlist(results, use.names = TRUE)
  out_dt[, gene_list := NULL]   # drop list column (not serialisable to TSV)
  out_dt[]
}


# ## ══════════════════════════════════════════════════════════════════════════════
# ##  Example usage
# ## ══════════════════════════════════════════════════════════════════════════════
# 
# ## -- load data ----------------------------------------------------------------
# dt <- fread("Pathway_enrichment_effector_genee_by_embedding_20260422.txt",
#             sep = "\t", quote = '"')
# 
# ## -- run pruning --------------------------------------------------------------
# pruned <- prune_pathways(
#   dt,
#   jaccard_threshold = 0.5,   # drop pathways sharing >50 % of genes with a
#   #   previously selected one
#   min_new_genes     = 1,     # must add at least 1 novel gene
#   rank_by           = "new_genes"
# )
# 
# ## -- inspect ------------------------------------------------------------------
# cat(sprintf("Input rows : %d\n", nrow(dt)))
# cat(sprintf("Pruned rows: %d\n", nrow(pruned)))
# cat(sprintf("Reduction  : %.1f %%\n",
#             100 * (1 - nrow(pruned) / nrow(dt))))
# 
# print(pruned[, .(emb, source, name, p_value,
#                  intersection_size, n_new_genes)])
# 
# ## -- save ---------------------------------------------------------------------
# fwrite(pruned, "pruned_pathways.tsv", sep = "\t")