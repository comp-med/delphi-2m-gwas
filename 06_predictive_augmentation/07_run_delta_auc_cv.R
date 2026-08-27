#!/usr/bin/env Rscript

## script to compute delta-AUC (exposure gain over the Delphi LP) in UKB
## now also returns a cross-fitted (out-of-sample) delta-AUC: delta.auc.cv
## CONTINUOUS exposures only -- categorical/questionnaire items are not considered
## Maik Pietzner 18/06/2026
rm(list = ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly = T)

## little options
options(stringsAsFactors = F)

## correct directory
setwd("<path_to_file>")

## packages needed
require(data.table)
require(arrow)

## import variables needed
target.pop <- args[1] ## 'all', 'men' or 'women'
target.ehr <- args[2] ## 'all' for the delta-AUC (full eval set)
## outcome of interest
outc       <- args[3]
## Delphi token (= column of the LP parquet) for this outcome
token      <- args[4]

## --> baseline-age column in ukb.dat used for the age x sex strata (match delphi.delta.auc) <-- ##
age.col    <- "age"

# ## for testing purpose
# target.pop <- "all"
# target.ehr <- "all"
# outc       <- "f.130008.0.0"
# token      <- "A04 Other bacterial intestinal infections"

cat("Computing delta-AUC for", outc, "(token", token, ") in", target.pop, "\n")
cat("--------------------------------------------------\n")

##################################
####   import relevant data   ####
##################################

## labels and exposure data (same inputs as the Cox driver)
lab.set <- fread("input/UKB.labels.prep.ICD10.code.Cox.models.20260626.txt")
ukb.dat <- fread("input/UKB.data.prep.ICD10.code.Cox.models.20260626.txt")
## the shortlist built by 04 (which exposures to run for this outcome)
short   <- fread("input/delphi.shortlist.txt")

#---------------------------------#
##-- transform to ease analysis --#
#---------------------------------#

## columns to deal with
cols    <- c(lab.set[type %in% c("character", "factor"), short_name], "centre")
## most frequent level as reference
ukb.dat[, (cols) := lapply(.SD, function(x) {
  levs <- names(sort(table(x), decreasing = TRUE))
  factor(x, levels = levs)
}), .SDcols = cols]

#---------------------------------#
##--    reduce to target pop.   --#
#---------------------------------#

## subset accordingly (primary care x sex) -- identical logic to 03
if(target.pop == "all" & target.ehr == "ehr"){
  ukb.dat <- ukb.dat[ primary.care == T]
}else if(target.pop == "all" & target.ehr == "no.ehr"){
  ukb.dat <- ukb.dat[ primary.care == F]
}else if(target.pop == "men" & target.ehr == "no.ehr"){
  ukb.dat <- ukb.dat[ primary.care == F & sex == "Male"]
}else if(target.pop == "men" & target.ehr == "ehr"){
  ukb.dat <- ukb.dat[ primary.care == T & sex == "Male"]
}else if(target.pop == "men" & target.ehr == "all"){
  ukb.dat <- ukb.dat[ sex == "Male"]
}else if(target.pop == "women" & target.ehr == "no.ehr"){
  ukb.dat <- ukb.dat[ primary.care == F & sex == "Female"]
}else if(target.pop == "women" & target.ehr == "ehr"){
  ukb.dat <- ukb.dat[ primary.care == T & sex == "Female"]
}else if(target.pop == "women" & target.ehr == "all"){
  ukb.dat <- ukb.dat[ sex == "Female"]
}

#---------------------------------#
##--    import outcome data     --#
#---------------------------------#

## survival data for this outcome (same prep as the Cox driver)
dt.cox  <- fread(paste0("input/", outc, ".UKB.first.occurrence.parsed.20260520.txt"))
## delete early and pre-existing cases
dt.cox  <- dt.cox[ pre.baseline.event == 0 & !(event.occurred == 1 & t.follow <= 365/2)]

#---------------------------------#
##--   import Delphi LP        --##
#---------------------------------#

## long LP table from the extraction step: columns  eid | token | delphi.lp
delphi.lp     <- fread("<path_to_file>",
                       select = c("eid", "delphi_lp_baseline", "name"))

## subset to token needed
lp.dat        <- delphi.lp[name == token, .(eid, delphi_lp_baseline)]
## edit names
names(lp.dat) <- c("f.eid", "delphi.lp")

## if something went wrong
if(nrow(lp.dat) == 0){
  cat("No Delphi LP matched token '", token, "' -- name mismatch, not sparsity.\n", sep = "")
  quit(save = "no", status = 0)
}

## delete what is no longer needed
rm(delphi.lp); gc(reset=T)

#---------------------------------#
##--         combine            --#
#---------------------------------#

## merge explicitly on f.eid (do NOT rely on shared-column auto-merge)
ukb.dat <- merge(ukb.dat, dt.cox,  by = "f.eid")
ukb.dat <- merge(ukb.dat, lp.dat,  by = "f.eid")

#---------------------------------#
##--   exposures for this outc --##
#---------------------------------#

## the shortlisted exposures for this outcome and population
expo.set <- short[ event == outc & sex == target.pop]

## CONTINUOUS exposures only -- categorical / questionnaire items are not considered
cont.set <- lab.set[ type %in% c("numeric", "integer"), short_name]
expo.set <- expo.set[ short_name %in% cont.set]
cat("Continuous exposures to run:", nrow(expo.set), "\n")

## guard: nothing to do
if(nrow(expo.set) == 0){
  cat("No continuous shortlisted exposures for this outcome x population -- exiting.\n")
  quit(save = "no", status = 0)
}

##################################
####        delta-AUC         ####
##################################

## one cross-fitted estimator for everything (INT, dl.cv, delphi.delta.auc.cv, biomarker.gain.cv)
source("../functions/delphi_delta_auc_cv.R")

## minimum incident cases per age x sex cell (full cohort -> 20; the global fallback rescues the rest)
min.cell <- 20L

#-------------------------------------#
##--      loop the exposures       --##
#-------------------------------------#

res.delta <- rbindlist(lapply(1:nrow(expo.set), function(x){
  
  ## cross-fitted delta-AUC over the Delphi LP, with component AUCs (auc.base.cv -> auc.comb.cv = delta)
  cv  <- delphi.delta.auc.cv(ukb.dat,
                             exposure      = expo.set$short_name[x],
                             dlp.col       = "delphi.lp",
                             time.col      = "t.follow",
                             event.col     = "event.occurred",
                             age.col       = age.col,
                             sex.col       = "sex",
                             inv.transform = expo.set$inv.transform[x],
                             horizon.years = 10,
                             age.bin       = 5,
                             n.folds       = 5L,
                             min.cell      = min.cell,
                             min.strata    = 2L)
  if(is.null(cv)) return(NULL)
  
  ## the exposure's OWN gain over demographics (age + sex), same machinery (auc.dem.cv -> auc.demB.cv = gain)
  gn  <- biomarker.gain.cv(ukb.dat,
                           exposure      = expo.set$short_name[x],
                           time.col      = "t.follow",
                           event.col     = "event.occurred",
                           age.col       = age.col,
                           sex.col       = "sex",
                           inv.transform = expo.set$inv.transform[x],
                           horizon.years = 10,
                           age.bin       = 5,
                           n.folds       = 5L,
                           min.cell      = min.cell,
                           min.strata    = 2L)
  if(is.null(gn)){
    gn <- data.table(gain.cv = NA_real_, gain.se.cv = NA_real_, gain.lci.cv = NA_real_, gain.uci.cv = NA_real_,
                     gain.pval.cv = NA_real_, auc.dem.cv = NA_real_, auc.demB.cv = NA_real_,
                     gain.n.strata.cv = NA_integer_, gain.mode = NA_character_)
  }else{
    gn <- gn[, .(gain.cv, gain.se.cv, gain.lci.cv, gain.uci.cv, gain.pval.cv,
                 auc.dem.cv, auc.demB.cv, gain.n.strata.cv, gain.mode)]
  }
  
  ## carry the cleaned name + outcome for collation
  cbind(event = outc, sex = target.pop,
        short_name = expo.set$short_name[x], short_name_new = expo.set$short_name_new[x],
        cv, gn)
}), fill = T)

#-------------------------------------#
##--          output results       --##
#-------------------------------------#

## nothing survived the per-stratum cell-count filter
if(is.null(res.delta) || nrow(res.delta) == 0){
  cat("No estimable delta-AUC (cell counts too small) -- exiting.\n")
  quit(save = "no", status = 0)
}

## write to file
write.table(
  res.delta,
  paste("output/delta_auc", outc, target.pop, target.ehr, "txt", sep = "."),
  row.names = F,
  sep = "\t",
  quote = F
)

cat("Finished file\n")
cat("--------------------------------------------------\n")