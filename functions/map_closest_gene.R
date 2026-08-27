library(data.table)

map_nearest_gene <- function(variants, genes, keep_ties = FALSE) {
  variants <- as.data.table(copy(variants))
  genes    <- as.data.table(copy(genes))
  
  # checks
  req_var  <- c("id", "CHROM", "GENPOS")
  req_gene <- c("CHROM", "start", "end", "gene_id", "gene_name")
  
  miss_var  <- setdiff(req_var, names(variants))
  miss_gene <- setdiff(req_gene, names(genes))
  
  if (length(miss_var)) {
    stop("variants is missing: ", paste(miss_var, collapse = ", "))
  }
  if (length(miss_gene)) {
    stop("genes is missing: ", paste(miss_gene, collapse = ", "))
  }
  
  # keep only usable genes
  genes <- genes[!is.na(CHROM) & !is.na(start) & !is.na(end),
                 .(CHROM, start, end, gene_id, gene_name)]
  
  # normalize types
  variants[, `:=`(
    CHROM = as.integer(as.character(CHROM)),
    GENPOS = as.integer(GENPOS)
  )]
  
  genes[, `:=`(
    CHROM = as.integer(as.character(CHROM)),
    start = as.integer(start),
    end   = as.integer(end)
  )]
  
  # preserve original row order
  variants[, row_id__ := .I]
  
  # process chromosome by chromosome
  res <- rbindlist(lapply(split(variants, by = "CHROM", keep.by = TRUE), function(vchr) {
    chr <- vchr$CHROM[1]
    gchr <- genes[CHROM == chr]
    
    # no genes on that chromosome
    if (nrow(gchr) == 0L) {
      return(vchr[, .(
        id, CHROM, GENPOS,
        nearest_gene_id = NA_character_,
        nearest_gene_name = NA_character_,
        distance = NA_integer_,
        row_id__
      )])
    }
    
    # sort genes for interval search
    setorder(gchr, start, end)
    
    out <- vector("list", nrow(vchr))
    
    for (k in seq_len(nrow(vchr))) {
      pos <- vchr$GENPOS[k]
      
      # genes containing the variant
      inside_idx <- which(gchr$start <= pos & gchr$end >= pos)
      
      if (length(inside_idx) > 0L) {
        hit <- gchr[inside_idx]
        
        if (keep_ties) {
          gene_id_val   <- paste(hit$gene_id, collapse = "|")
          gene_name_val <- paste(hit$gene_name, collapse = "|")
        } else {
          o <- order(hit$gene_name, hit$gene_id)
          gene_id_val   <- hit$gene_id[o][1]
          gene_name_val <- hit$gene_name[o][1]
        }
        
        out[[k]] <- data.table(
          id = vchr$id[k],
          CHROM = chr,
          GENPOS = pos,
          nearest_gene_id = gene_id_val,
          nearest_gene_name = gene_name_val,
          distance = 0L,
          row_id__ = vchr$row_id__[k]
        )
        
      } else {
        # nearest gene to the left: largest end < pos
        left_idx <- which(gchr$end < pos)
        # nearest gene to the right: smallest start > pos
        right_idx <- which(gchr$start > pos)
        
        left_dist <- if (length(left_idx)) pos - max(gchr$end[left_idx]) else Inf
        right_dist <- if (length(right_idx)) min(gchr$start[right_idx]) - pos else Inf
        
        if (left_dist < right_dist) {
          hit <- gchr[left_idx][end == max(end)]
        } else if (right_dist < left_dist) {
          hit <- gchr[right_idx][start == min(start)]
        } else {
          # tie
          if (is.infinite(left_dist) && is.infinite(right_dist)) {
            hit <- gchr[0]
          } else {
            hit_left  <- gchr[left_idx][end == max(end)]
            hit_right <- gchr[right_idx][start == min(start)]
            hit <- unique(rbind(hit_left, hit_right, fill = TRUE))
          }
        }
        
        if (nrow(hit) == 0L) {
          gene_id_val   <- NA_character_
          gene_name_val <- NA_character_
          dist_val      <- NA_integer_
        } else if (keep_ties) {
          gene_id_val   <- paste(hit$gene_id, collapse = "|")
          gene_name_val <- paste(hit$gene_name, collapse = "|")
          dist_val      <- as.integer(min(left_dist, right_dist))
        } else {
          o <- order(hit$gene_name, hit$gene_id)
          gene_id_val   <- hit$gene_id[o][1]
          gene_name_val <- hit$gene_name[o][1]
          dist_val      <- as.integer(min(left_dist, right_dist))
        }
        
        out[[k]] <- data.table(
          id = vchr$id[k],
          CHROM = chr,
          GENPOS = pos,
          nearest_gene_id = gene_id_val,
          nearest_gene_name = gene_name_val,
          distance = dist_val,
          row_id__ = vchr$row_id__[k]
        )
      }
    }
    
    rbindlist(out)
  }), use.names = TRUE, fill = TRUE)
  
  setorder(res, row_id__)
  res[, row_id__ := NULL]
  res[]
}