##################################################
## function to extract genomic region using polars
## provided by Carl Beuchel + augmented with ChatGPT

# Formatted file directory: <path_to_file>

extract_genomic_data <- function(
    parquet_location = fs::path(),
    chr = NULL,
    bp_min = NULL,
    bp_max = NULL,
    snp_list = NULL,
    chr_column = "CHR_ENSEMBL",
    bp_column = "BP",
    snp_column = "SNP",
    id_ucsc_column = "ID_UCSC"
) {
  
  # Helper: return an empty data.table with same schema as parquet
  empty_dt_like <- function(parquet_location) {
    pl_empty <- polars::pl$scan_parquet(parquet_location)$limit(0L)$collect()
    data.table::as.data.table(pl_empty)
  }
  
  ## -----------------------------
  ## Basic file and mode checks
  ## -----------------------------
  stopifnot(
    "Input .parquet file does not exist" = fs::file_exists(parquet_location)
  )
  
  region_mode <- !is.null(chr) && !is.null(bp_min) && !is.null(bp_max)
  snp_mode    <- !is.null(snp_list)
  
  if (!region_mode && !snp_mode) {
    stop("Must provide either (chr, bp_min, bp_max) for region query OR snp_list for SNP query")
  }
  
  if (region_mode && snp_mode) {
    stop("Cannot query by both region and SNP list simultaneously. Please choose one method.")
  }
  
  ## -----------------------------
  ## REGION-BASED QUERY
  ## -----------------------------
  if (region_mode) {
    
    # Normalise requested chromosome (strip "chr", lower-case)
    chr_norm <- gsub("^chr", "", trimws(tolower(as.character(chr))))
    checkmate::assert_character(
      as.character(chr_norm),
      min.chars = 1,
      max.chars = 2,
      len = 1,
      any.missing = FALSE
    )
    
    bp_min <- assertCount(as.integer(bp_min), positive = TRUE, coerce = TRUE)
    bp_max <- assertCount(as.integer(bp_max), positive = TRUE, coerce = TRUE)
    stopifnot("bp_min must be smaller than bp_max" = bp_min < bp_max)
    
    # Get all chromosomes present in the file
    all_chromosomes <- tryCatch(
      {
        polars::pl$scan_parquet(parquet_location)$
          select(chr_column)$
          unique()$
          collect()$
          to_series()$
          to_r_vector()
      },
      error = function(e) {
        warning("Failed to read chromosome column from parquet: ", conditionMessage(e),
                "\nReturning empty result.")
        return(character())  # handled below
      }
    )
    
    if (length(all_chromosomes) == 0L) {
      warning("No rows or no chromosome information found in parquet file – returning empty data.table with schema.")
      return(empty_dt_like(parquet_location))
    }
    
    # Normalise file chromosomes in the same way as the requested one
    all_chr_norm <- gsub("^chr", "", trimws(tolower(as.character(all_chromosomes))))
    
    # Find corresponding value in the file for the requested chromosome
    match_idx <- which(all_chr_norm == chr_norm)
    
    if (!length(match_idx)) {
      warning(paste0(
        "Requested chromosome '", chr, "' not found in file. ",
        "Available chromosomes (normalised): ",
        paste(sort(unique(all_chr_norm)), collapse = ", "),
        ". Returning empty data.table with schema."
      ))
      return(empty_dt_like(parquet_location))
    }
    
    # Use the actual value as stored in the file (e.g. "1" or "chr1")
    chr_value_in_file <- all_chromosomes[match_idx[1]]
    
    # Build and run Polars query
    filter_result <- tryCatch(
      {
        polars::pl$scan_parquet(parquet_location)$
          filter(
            polars::pl$col(chr_column) == chr_value_in_file
          )$
          # Cast BP to UInt64 but allow failures (non-numeric become NULL)
          with_columns(
            polars::pl$col(bp_column)$cast(polars::pl$UInt64, strict = FALSE)$alias(bp_column)
          )$
          # Drop non-numeric BP values and apply range filter
          filter(
            polars::pl$col(bp_column)$is_not_null() &
              polars::pl$col(bp_column)$ge(bp_min) &
              polars::pl$col(bp_column)$le(bp_max)
          )$
          collect()
      },
      error = function(e) {
        warning("Region query failed in Polars: ", conditionMessage(e),
                "\nReturning empty data.table with schema.")
        return(NULL)
      }
    )
    
    if (is.null(filter_result)) {
      return(empty_dt_like(parquet_location))
    }
    
    if (nrow(filter_result) == 0L) {
      warning(paste0(
        "Chromosome '", chr, "' is present, but no variants in BP range [",
        bp_min, ", ", bp_max, "]. Returning empty data.table with schema."
      ))
      return(empty_dt_like(parquet_location))
    }
    
    return(data.table::as.data.table(filter_result))
  }
  
  ## -----------------------------
  ## SNP-BASED QUERY
  ## -----------------------------
  if (snp_mode) {
    
    ## Validate SNP list
    if (!is.character(snp_list) && !is.vector(snp_list)) {
      stop("snp_list must be a character vector")
    }
    
    if (length(snp_list) == 0) {
      stop("snp_list is empty")
    }
    
    ## Remove duplicates and NAs
    snp_list <- unique(na.omit(snp_list))
    
    ## Detect SNP format: rsID or ID_UCSC
    sample_snps <- head(snp_list, min(10, length(snp_list)))
    
    is_rsid <- any(grepl("^rs\\d+", sample_snps))
    is_ucsc <- any(grepl("^chr.+:.+_.+_.+", sample_snps))
    
    ## Determine which column to use
    if (is_rsid && !is_ucsc) {
      query_column <- snp_column
      message(paste0("Querying by rsID using column: ", snp_column))
    } else if (is_ucsc && !is_rsid) {
      query_column <- id_ucsc_column
      message(paste0("Querying by ID_UCSC using column: ", id_ucsc_column))
    } else if (is_rsid && is_ucsc) {
      # Mixed format - try both columns
      warning("SNP list contains mixed formats (rsID and ID_UCSC). Will query both columns.")
      
      result_rsid <- tryCatch(
        {
          polars::pl$scan_parquet(parquet_location)$
            filter(
              polars::pl$col(snp_column)$is_in(
                polars::pl$lit(snp_list[grepl("^rs\\d+", snp_list)])
              )
            )$
            collect()
        },
        error = function(e) {
          warning("rsID query failed in Polars: ", conditionMessage(e),
                  "\nTreating as no matches for rsIDs.")
          return(NULL)
        }
      )
      
      result_ucsc <- tryCatch(
        {
          polars::pl$scan_parquet(parquet_location)$
            filter(
              polars::pl$col(id_ucsc_column)$is_in(
                polars::pl$lit(snp_list[grepl("^chr.+:.+_.+_.+", snp_list)])
              )
            )$
            collect()
        },
        error = function(e) {
          warning("ID_UCSC query failed in Polars: ", conditionMessage(e),
                  "\nTreating as no matches for ID_UCSC.")
          return(NULL)
        }
      )
      
      if (is.null(result_rsid) && is.null(result_ucsc)) {
        warning("Both rsID and ID_UCSC sub-queries failed. Returning empty data.table with schema.")
        return(empty_dt_like(parquet_location))
      }
      
      non_null <- Filter(Negate(is.null), list(result_rsid, result_ucsc))
      
      combined <- polars::pl$concat(non_null)$unique()
      
      if (nrow(combined) == 0L) {
        warning("None of the requested SNPs were found in the file. Returning empty data.table with schema.")
        return(empty_dt_like(parquet_location))
      }
      
      return(data.table::as.data.table(combined))
    } else {
      stop("Could not detect SNP format. Expected rsID (e.g., 'rs123456') or ID_UCSC (e.g., 'chr1:12345_A_G')")
    }
    
    ## Filter based on SNP list (single format)
    filter_result <- tryCatch(
      {
        polars::pl$scan_parquet(parquet_location)$
          filter(
            polars::pl$col(query_column)$is_in(polars::pl$lit(snp_list))
          )$
          collect()
      },
      error = function(e) {
        warning("SNP query failed in Polars: ", conditionMessage(e),
                "\nReturning empty data.table with schema.")
        return(NULL)
      }
    )
    
    if (is.null(filter_result)) {
      return(empty_dt_like(parquet_location))
    }
    
    n_found     <- nrow(filter_result)
    n_requested <- length(snp_list)
    
    if (n_found == 0L) {
      warning(paste0(
        "No matches found for the provided SNP list (", n_requested,
        " SNPs). Returning empty data.table with schema."
      ))
      return(empty_dt_like(parquet_location))
    } else if (n_found < n_requested) {
      warning(paste0(
        "Found ", n_found, " out of ", n_requested,
        " requested SNPs. Returning only the matched SNPs."
      ))
    } else {
      message(paste0("Found all ", n_found, " requested SNPs"))
    }
    
    return(data.table::as.data.table(filter_result))
  }
}

extract_genomic_region <- function(
    parquet_location = fs::path(),
    chr = character(),
    bp_min = integer(),
    bp_max = integer(),
    chr_column = "CHR_ENSEMBL",
    bp_column = "BP"
) {
  extract_genomic_data(
    parquet_location = parquet_location,
    chr = chr,
    bp_min = bp_min,
    bp_max = bp_max,
    chr_column = chr_column,
    bp_column = bp_column
  )
}
