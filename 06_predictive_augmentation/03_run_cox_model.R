#!/usr/bin/env Rscript

## script to run Cox regression models in UKB
## Maik Pietzner 20/05/2026
rm(list = ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly = T)

## little options
options(stringsAsFactors = F)
## avoid conversion of numbers
# options(scipen = 1)
# options(rgl.useNULL = TRUE)
# print(R.Version())

## correct directory
setwd("<path_to_file>")

## packages needed
require(data.table)
require(doMC)
require(survival)

## import variables needed
target.pop <- args[1] ## can be 'all' or subset
target.ehr <- args[2]
## outcome of interest
outc       <- args[3]
## adjustment variables ('supplied as string')
adj        <- args[4]

# ## for testing purpose
# target.pop <- "all" ## can be 'all' or subset
# target.ehr <- "all"
# outc       <- "f.130016.0.0"
# adj        <- "+ age + sex + centre + smoking + alcohol"

cat("Performing Cox models for", outc, "in", target.pop, "using", adj, "\n")
cat("--------------------------------------------------\n")

##################################
####   import relevant data   ####
##################################

## import labels of variables of interest
lab.set <- fread("input/UKB.labels.prep.ICD10.code.Cox.models.20260626.txt")
## import data set
ukb.dat <- fread("input/UKB.data.prep.ICD10.code.Cox.models.20260626.txt")

# ## pull in phecode labels to account for sex-specific diseases
# lab.phe <- fread("<path_to_file>", sep="\t", header=T)
# ## add identifier to match with data set
# lab.phe[, id := paste0("bin_", phecode)]
# ## replace missing ones
# lab.phe[, sex := ifelse(sex == "", "Both", sex)]
# 
# ## create indicator for INV
# lab.set[, inv.transform := category %in% c("Biomarker", "Blood cell counts", "Body composition", "Bone", "Technical", "Cardiovascular", "Pollution")]
# ## add information on possible sex-specific outcomes
# lab.set <- merge(lab.set, lab.phe[, .(id, sex)], by.x = "short_name", by.y = "id", all.x = T)
# ## replace NA accordingly
# lab.set[, sex := ifelse(!is.na(sex), sex, "Both")]

#---------------------------------#
##-- transform to ease analysis --#
#---------------------------------#

## columns to deal with
cols    <- c(lab.set[type %in% c("character", "factor"), short_name], "centre")
## apply transformation
ukb.dat[, (cols) := lapply(.SD, function(x) {
  levs <- names(sort(table(x), decreasing = TRUE))
  factor(x, levels = levs)
}), .SDcols = cols]

#---------------------------------#
##--    reduce to target pop.   --#
#---------------------------------#

## subset accordingly (primary care x sex)
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

## import survival data
dt.cox   <- fread(paste0("input/", outc, ".UKB.first.occurrence.parsed.20260520.txt"))
## delete early and pre-existing cases
dt.cox   <- dt.cox[ pre.baseline.event == 0 & !(event.occurred == 1 & t.follow <= 365/2)]
## import outcome label
lab.outc <- fread("input/Event.count.incident.ICD.codes.UKB.20260616.txt")

#---------------------------------#
##--         combine            --#
#---------------------------------#

## create combined data set
ukb.dat <- merge(ukb.dat, dt.cox)

#---------------------------------#
##--   attach Delphi LP         --#
#---------------------------------#

## long LP table from the extraction step: columns  eid | token | delphi.lp
delphi.lp  <- fread("<path_to_file>",
                    select = c("eid", "delphi_lp_baseline", "name"))
## import information to ease mapping to token
tmp.delphi <- fread("<path_to_file>")
## add ICD-10 code
tmp.delphi[, icd10.code := {
  m <- regexpr("[A-Z][0-9]{2}", Name)
  fifelse(m == -1L, NA_character_, substr(Name, m, m + 2L))
}]
## add name to count data
lab.outc   <- merge(lab.outc, tmp.delphi[ !is.na(icd10.code), .(icd10.code, Name)])
## get the token of interest
this.token <- lab.outc[ event == outc ]$Name

## extract only what is needed for this iteration
if (length(this.token) == 1) {
  ## get what is of interest
  lp.k    <- delphi.lp[name == this.token, .(eid, delphi_lp_baseline)]
  ukb.dat <- merge(ukb.dat, lp.k, by.x = "f.eid", by.y = "eid", all.x = TRUE)   # use YOUR id column name
  d.var   <- "delphi_lp_baseline"
} else {
  d.var   <- NULL                                   # falls back to age/sex-only behaviour
}

##################################
####         Cox models       ####
##################################

## import function to do so
source("../functions/run_cox_models.R")

#-------------------------------------#
##--           run the model       --##
#-------------------------------------#

## run model depending on target population
res.cox <- run.cox( ukb.dat,
                    "t.follow",
                    "event.occurred",
                    hide.adj = T,
                    verbose = T,
                    adj,
                    lab.set,
                    delphi.var = d.var,
                    centre.var = "centre", 
                    centre.method = "strata")

# ## drop not needed; may well need additional cleaning...
# res.cox <- res.cox[
#   !(short_name_new %in%
#       c("age", "sexMale", 
#         grep("smoking", res.cox$short_name_new, value = T),
#         grep("alcohol", res.cox$short_name_new, value = T), 
#         grep("centre", res.cox$short_name_new, value = T)))
# ]

## edit names to match with labels; account for transformation
res.cox[, short_name_new := gsub("^INT\\((.*)\\)$", "\\1", short_name_new)]
res.cox[, short_name_new := gsub("[[:punct:]]| ", "_", short_name_new)]

## add the outcome for convience
res.cox[, event := outc]

#-------------------------------------#
##--          output results       --##
#-------------------------------------#

## write to file
write.table(
  res.cox,
  paste(
    "output/cox",
    outc,
    target.pop,
    target.ehr,
    gsub("\\s*\\+\\s*", ".", gsub("^\\+\\s*", "", adj)),
    "txt",
    sep = "."
  ),
  row.names = F,
  sep = "\t"
)

cat("Finished file\n")
cat("--------------------------------------------------\n")
