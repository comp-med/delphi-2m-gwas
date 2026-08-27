#######################################################
#### collate input for the delta-AUC job array      ####
#### Maik Pietzner                      16/06/2026 ####
#######################################################

rm(list=ls())
setwd("<path_to_file>")
options(stringsAsFactors = F)
load(".RData")

## --> packages needed <-- ##
require(data.table)

###############################################
####       collated Cox results (v3)       ####
###############################################

## res.cox.minimal must come from the *v3* Cox run, i.e. carry c.delta.delphi
## (exposure gain over the Delphi LP). If the .RData still holds the pre-v3
## collation, re-run 01's collation step first. One row per outcome x exposure x sex.
res.cox.delphi      <- copy(res.cox.minimal)

#---------------------------------#
##--      exposure labels      --##
#---------------------------------#

## bring in the data-column name (short_name) and the INT flag; the delta-AUC
## function needs the DATA column, not the cleaned short_name_new
lab.set             <- fread("UKB.labels.prep.ICD10.code.Cox.models.20260626.txt")
## clean key to match res.cox.minimal$short_name_new
lab.set[, short_name_new := gsub("[[:punct:]]| ", "_", short_name)]
## keep only what we need
lab.key             <- unique(lab.set[, .(short_name, short_name_new, category, inv.transform)])

#---------------------------------#
##--   select the shortlist    --##
#---------------------------------#

## thresholds are deliberately exposed -- edit to taste
p.thresh            <- 1.9e-7   ## study-wide significance used throughout
delta.thresh        <- 0.005    ## minimal gain over Delphi worth confirming

## keep significant, positive-gain rows
res.cox.delphi      <- res.cox.delphi[ !is.na(c.delta.delphi) &
                                       c.delta.delphi > delta.thresh &
                                       `Pr(>|z|)` < p.thresh]

## drop exposures Delphi could never sensibly ingest as a clean signal
## (definitional / reverse-causal / iatrogenic), reusing the AI annotation
res.annot           <- fread("Results_Cox_models_first_occurrence_ICD10_UKB_annotated_20260615.txt")
res.cox.delphi      <- merge(res.cox.delphi,
                             res.annot[, .(event, short_name_new, sex, novelty, mechanism_category, adr_flag)],
                             by = c("event", "short_name_new", "sex"), all.x = T)
res.cox.delphi      <- res.cox.delphi[ !(novelty %in% "Artefact") &
                                       (adr_flag == F | is.na(adr_flag)) &
                                       !grepl("Definitional|Reverse causation", mechanism_category)]

## the delta-AUC function takes a single continuous/binary exposure -> keep those
## (multi-level categoricals would need the function extended; drop them here)
res.cox.delphi      <- merge(res.cox.delphi, lab.key, by = c("short_name_new"), all.x = T)
res.cox.delphi      <- res.cox.delphi[ inv.transform == T |
                                       category %in% c("Biomarker", "Blood cell counts", "Drugs", "Disease",
                                                       "Body composition", "Bone", "Pulmonary", "Cardiovascular")]

#---------------------------------#
##--   map outcome -> token    --##
#---------------------------------#

## token = the column name under which the python LP is stored for this outcome.
## here we key the LP wide-parquet by ICD-10 code; if your postdoc keyed it by
## the Delphi token index/name instead, swap icd10.code for that identifier.
res.cox.delphi[, token := icd10.code]
## drop outcomes with no Delphi token (not in the vocabulary)
res.cox.delphi      <- res.cox.delphi[ !is.na(token)]

###############################################
####             write inputs              ####
###############################################

#---------------------------------#
##--     per-pair shortlist    --##
#---------------------------------#

## read by 06 to know which exposures (and INT flag) to run for each outcome
shortlist           <- unique(res.cox.delphi[, .(event, icd10.code, token,
                                                  short_name, short_name_new,
                                                  sex, inv.transform, c.delta.delphi)])
## one exposure can win in 'all' and a sex-stratum; keep the most informative row per (event, exposure, sex)
shortlist           <- shortlist[ order(event, short_name, sex, -c.delta.delphi)]
fwrite(shortlist, "delphi.shortlist.txt", sep = "\t", row.names = F, quote = F, na = NA)
## n pairs
nrow(shortlist)

#---------------------------------#
##--     job-array index       --##
#---------------------------------#

## one task = one outcome x population (06 loops the shortlisted exposures inside).
## target.ehr fixed to 'all' -- discrimination is judged on the full eval set,
## matching how the Delphi AUC itself was computed (not the EHR split).
arr                 <- unique(shortlist[, .(target_pop = sex, target_ehr = "all", outc = event, token)])
## NB: NO header row -- 05's awk uses 'NR == var', and run order must be stable
fwrite(arr, "Input.delta.auc.txt", sep = "\t", row.names = F, quote = F, col.names = F, na = NA)
## tell us the array size to put in 05_submit_delta_auc.sh
cat("Input.delta.auc.txt rows (set --array=1-N to):", nrow(arr), "\n")
## e.g. --array=1-N
