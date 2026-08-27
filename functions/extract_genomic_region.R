##################################################
## function to extract genomic region using polars
## provided by Carl Beuchel

# Formatted file directory: <path_to_file>

extract_genomic_region <- function(
    parquet_location = fs::path(),
    chr = character(),
    bp_min = integer(),
    bp_max = integer(),
    chr_column = "CHR_ENSEMBL",
    bp_column = "BP"
) {
  chr <- gsub("^chr", "", trimws(tolower(as.character(chr))))
  checkmate::assert_character(
    as.character(chr),
    min.chars = 1,
    max.chars = 2,
    len = 1,
    any.missing = FALSE
  )
  bp_min <- assertCount(as.integer(bp_min), positive = TRUE, coerce = TRUE)
  bp_max <- assertCount(as.integer(bp_max), positive = TRUE, coerce = TRUE)
  stopifnot(
    "Input .parquet file does not exist" = fs::file_exists(parquet_location)
  )
  stopifnot(
    "`bp_min` must be smaller than `bp_max`" = bp_min < bp_max
  )
  all_chromosomes <- polars::pl$scan_parquet(parquet_location)$select(
    (chr_column)
  )$unique()$collect()$to_series()$to_r_vector()
  if (!any(chr %in% all_chromosomes)) {
    stop(paste0(
      "Chromosome not found. Available chromosomes:",
      paste(sort(as.integer(all_chromosomes)), collapse = ", ")
    ))
  }
  
  ## filter based on genomic region
  filter_result <- polars::pl$scan_parquet(parquet_location)$cast(
    BP = pl$UInt64
  )$filter(pl$col(chr_column) == chr)$filter(pl$col(bp_column)$ge(
    bp_min
  ))$filter(pl$col(bp_column)$le(bp_max))$collect()
  return(data.table::as.data.table(filter_result))
}
