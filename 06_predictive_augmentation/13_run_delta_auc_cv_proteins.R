#!/usr/bin/env Rscript

## delta-AUC (protein gain over the Delphi LP) in the UKB Olink subcohort (~45k)
## cross-fitted; global age/sex-adjusted fallback for sparse outcomes
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

## import variables needed
target.pop <- args[1] ## 'all', 'men' or 'women'
target.ehr <- args[2] ## 'all' for the full eval set
outc       <- args[3] ## outcome of interest

# ## for testing purpose
# target.pop <- "all"; target.ehr <- "all"; outc <- "f.131286.0.0"

cat("Computing protein delta-AUC for", outc, "in", target.pop, "\n")
cat("--------------------------------------------------\n")

##################################
####   import relevant data   ####
##################################

## covariates only (proteins carry the exposures)
ukb.dat  <- fread("input/UKB.data.prep.ICD10.code.Cox.models.20260626.txt",
                  select = c("f.eid", "age", "sex", "centre", "primary.care"))

#---------------------------------#
##--    reduce to target pop.   --#
#---------------------------------#

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
##--    import protein data     --#
#---------------------------------#

ukb.prot <- fread("input/UKB.Olink.delphi.2M.LP.20260618.txt")
lab.prot <- fread("input/Olink.proteins.delphi.2M.LP.20260618.txt")
lab.prot[, short_name := id]
lab.prot[, inv.transform := F]                    ## already applied

## inner-merge -> retains the Olink subcohort only (~45k)
ukb.dat  <- merge(ukb.dat, ukb.prot, by.x = "f.eid", by.y = "eid")
rm(ukb.prot); gc(reset = T)

#---------------------------------#
##--    import outcome data     --#
#---------------------------------#

dt.cox   <- fread(paste0("input/", outc, ".UKB.first.occurrence.parsed.20260520.txt"))
dt.cox   <- dt.cox[ pre.baseline.event == 0 & !(event.occurred == 1 & t.follow <= 365/2)]
lab.outc <- fread("input/Event.count.incident.ICD.codes.UKB.20260616.txt")
ukb.dat  <- merge(ukb.dat, dt.cox)

#---------------------------------#
##--   attach Delphi LP         --#
#---------------------------------#

delphi.lp  <- fread("<path_to_file>",
                    select = c("eid", "delphi_lp_baseline", "name"))
tmp.delphi <- fread("<path_to_file>")
tmp.delphi[, icd10.code := {
  m <- regexpr("[A-Z][0-9]{2}", Name)
  fifelse(m == -1L, NA_character_, substr(Name, m, m + 2L))
}]
lab.outc   <- merge(lab.outc, tmp.delphi[ !is.na(icd10.code), .(icd10.code, Name)])
this.token <- lab.outc[ event == outc ]$Name

if(length(this.token) != 1){
  cat("No unique Delphi token for", outc, "-- exiting.\n"); quit(save = "no", status = 0)
}
lp.k     <- delphi.lp[ name == this.token, .(eid, delphi_lp_baseline)]
ukb.dat  <- merge(ukb.dat, lp.k, by.x = "f.eid", by.y = "eid")
rm(delphi.lp); gc(reset = T)

#---------------------------------#
##--   proteins for this outc   --#
#---------------------------------#

## shortlist (built from the protein Cox results) if present, else all proteins
f.short  <- "input/protein.shortlist.txt"     ## cols: event, short_name, sex
if(file.exists(f.short)){
  short    <- fread(f.short)
  prot.set <- short[ event == outc & sex == target.pop, short_name]
} else {
  prot.set <- lab.prot$short_name
}
prot.set <- intersect(prot.set, names(ukb.dat))
cat("Proteins to test:", length(prot.set), "\n")
if(length(prot.set) == 0){
  cat("No proteins for this outcome x population -- exiting.\n"); quit(save = "no", status = 0)
}

##################################
####        delta-AUC         ####
##################################

## one cross-fitted estimator for everything (INT, dl.cv, delphi.delta.auc.cv, biomarker.gain.cv)
source("../functions/delphi_delta_auc_cv.R")

## minimum incident cases per age x sex cell (45k subcohort -> 10; the global fallback rescues the rest)
min.cell <- 10L

#-------------------------------------#
##--      loop the proteins        --##
#-------------------------------------#

res.delta <- rbindlist(lapply(prot.set, function(p){
  
  ## cross-fitted delta-AUC over the Delphi LP, with component AUCs (auc.base.cv -> auc.comb.cv = delta)
  cv  <- delphi.delta.auc.cv(ukb.dat,
                             exposure      = p,
                             dlp.col       = "delphi_lp_baseline",
                             time.col      = "t.follow",
                             event.col     = "event.occurred",
                             age.col       = "age",
                             sex.col       = "sex",
                             inv.transform = F,
                             horizon.years = 10,
                             age.bin       = 5,
                             n.folds       = 5L,
                             min.cell      = min.cell,
                             min.strata    = 2L)
  if(is.null(cv)) return(NULL)
  
  ## the protein's OWN gain over demographics (age + sex), same machinery (auc.dem.cv -> auc.demB.cv = gain)
  gn  <- biomarker.gain.cv(ukb.dat,
                           exposure      = p,
                           time.col      = "t.follow",
                           event.col     = "event.occurred",
                           age.col       = "age",
                           sex.col       = "sex",
                           inv.transform = F,
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
  
  cbind(event = outc, sex = target.pop, cv, gn)
}), fill = T)


#-------------------------------------#
##--          output results       --##
#-------------------------------------#

if(is.null(res.delta) || nrow(res.delta) == 0){
  cat("No estimable delta-AUC -- exiting.\n"); quit(save = "no", status = 0)
}

write.table(
  res.delta,
  paste("output/delta_auc_proteins", outc, target.pop, target.ehr, "txt", sep = "."),
  row.names = F, sep = "\t", quote = F
)

cat("Finished file\n")
cat("--------------------------------------------------\n")