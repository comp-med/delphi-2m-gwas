#######################################################
#### phenome-wide association testing first occur. ####
#### Maik Pietzner                      13/05/2026 ####
#######################################################

rm(list=ls())
setwd("<path_to_file>")
options(stringsAsFactors = F)
load(".RData")

## --> packages needed <-- ##
require(data.table)
require(arrow)
require(readxl)
require(Hmisc)
require(tidyverse)
require(doMC)
require(susieR)

###############################################
####            import ukb data            ####
###############################################

## import labels of variables of interest
lab.set <- fread("../../01_phenotype_preparation//data/UKB.labels.Delphi.embeddings.20260129.txt")
## import data set
ukb.dat <- fread("../../01_phenotype_preparation//data/UKB.data.prep.Delphi.embeddings.20260129.txt")

#---------------------------------#
##-- transform to ease analysis --#
#---------------------------------#

## convert to data frame to ease coding
ukb.dat <- as.data.frame(ukb.dat)
## choose most frequent factor as a reference for all other traits
for(j in c(lab.set[ type %in% c("character", "factor")]$short_name, "centre")){
  ## get frequency
  jj            <- table(ukb.dat[, j])
  ## redefine
  ukb.dat[, j] <- factor(ukb.dat[, j], levels = names(jj[order(jj, decreasing = T)]), ordered = F)
}
## convert back
ukb.dat <- as.data.table(ukb.dat)

## create second set with extended labels
lab.ext <- lapply(1:nrow(lab.set), function(x){
  ## check whether factor
  if(lab.set$type[x] == "character" | lab.set$short_name[x] == "centre" | lab.set$type[x] == "factor"){
    ## get the relevant variable
    n <- lab.set$short_name[x]
    ## get all levels
    l <- levels(as.factor(unlist(ukb.dat[, ..n])))
    print(l)
    ## replace odd characters
    l <- gsub("[[:punct:]]| ", "_", l)
    ## extend labels accordingly
    return(data.table(lab.set[x,], short_name_new = paste0(n, l[-1]),
                      label_new = paste(lab.set$label[x], l[-1]), reference = l[1]))
  }else{
    ## extend labels accordingly
    return(data.table(lab.set[x,], short_name_new = gsub("[[:punct:]]| ", "_", lab.set$short_name[x]), label_new = lab.set$label[x], reference = NA))
  }
})
lab.ext <- rbindlist(lab.ext)

## create indicator for INV
lab.set[, inv.transform := category %in% c("Biomarker", "Blood cell counts", "Body composition", "Bone", "Technical", "Cardiovascular", "Pollution")]
lab.ext[, inv.transform := category %in% c("Biomarker", "Blood cell counts", "Body composition", "Bone", "Technical", "Cardiovascular", "Pollution")]

#-------------------------#
##-- drop phecode data --##
#-------------------------#

## columns to be deleted
jj <- intersect(names(ukb.dat), 
                c(lab.set[ category == "Diseases"]$short_name,
                  gsub("bin", "date",lab.set[ category == "Diseases"]$short_name)))

## now from the data
ukb.dat[, (jj) := NULL]

## drop from labels next
lab.set <- lab.set[ category != "Diseases"]
lab.ext <- lab.ext[ category != "Diseases"]

## do some cleaning
gc(reset=T)

###############################################
####          first occurrence data        ####
###############################################

#------------------------------#
##--       import data      --##
#------------------------------#

## import first occurrence: filtered by Wenhuan 
ukb.first.occurrence <- fread("<path_to_file>")

## import labels
lab.first.occurrence <- fread("<path_to_file>")
## subset to what is included in data
lab.first.occurrence <- lab.first.occurrence[ id.ukbb %in% names(ukb.first.occurrence)]
## whether date or source column
lab.first.occurrence[, source.column := description %in% grep("Source of", lab.first.occurrence$description, value = T)]
lab.first.occurrence[, date.column := description %in% grep("Date ", lab.first.occurrence$description, value = T)]

## drop 'source' columns
ukb.first.occurrence <- ukb.first.occurrence[, c("f.eid", lab.first.occurrence[ date.column == T]$id.ukbb), with = F]

#------------------------------#
##--    prep. for survival  --##
#------------------------------#

## convert to wide format to minimize data storage
ukb.first.occurrence <- melt(ukb.first.occurrence,
                             id.vars       = "f.eid",
                             measure.vars  = lab.first.occurrence[ date.column == T]$id.ukbb,
                             variable.name = "event",
                             value.name    = "date",
                             na.rm         = TRUE)

## add baseline date: this will also drop participants that have withdrawn consent and have no genetic data
ukb.first.occurrence <- merge(ukb.first.occurrence, ukb.dat[, .(f.eid, baseline_date)])
## generate an event code
ukb.first.occurrence[, event.code := ifelse(date > baseline_date, 1, -1)]

## import death as competing event
ukb.death            <- fread("<path_to_file>")
## make unique
ukb.death            <- unique(ukb.death[, .(eid, date_of_death)])
## convert date
ukb.death[, date_of_death := as.IDate(date_of_death, format = "%d/%m/%Y")]
## augment to all participants
ukb.death            <- merge(ukb.death, ukb.dat[, .(f.eid, baseline_date)], 
                              by.x = "eid", by.y = "f.eid", all.y = T)
## adjust columns
ukb.death[, end.date.observation := as.IDate(ifelse(!is.na(date_of_death), as.character(date_of_death), "2020-08-31"))]

## derive full cross-join of all cohort participants × all observed event types
dt.skeleton          <- CJ(f.eid = ukb.death$eid,
                           event = ukb.first.occurrence[, unique(event)],
                           unique = TRUE)

## Merge Follow-up and Baseline Dates onto Skeleton
dt.cox               <- merge(dt.skeleton,
                              ukb.death[, .(f.eid = eid, baseline_date, end.date.observation)],
                              by = "f.eid", all.x = TRUE)

## isolate post-baseline (incident) events with their date
dt.post              <- ukb.first.occurrence[event.code ==  1L, .(f.eid, event, event_date = date)]

## isolate pre-baseline (prevalent) events as binary covariate
dt.pre               <- ukb.first.occurrence[event.code == -1L, .(f.eid, event, pre.baseline.event = 1L)]

## merge incident event dates
dt.cox               <- merge(dt.cox, dt.post, by = c("f.eid", "event"), all.x = TRUE)

## merge prevalent event indicator; unmatched participants coded as event-naive
dt.cox               <- merge(dt.cox, dt.pre,  by = c("f.eid", "event"), all.x = TRUE)
dt.cox[is.na(pre.baseline.event), pre.baseline.event := 0L]

## binary incident event status: 1 = post-baseline event observed, 0 = censored
dt.cox[, event.occurred := as.integer(!is.na(event_date))]

## follow-up time in days from baseline to earliest of: incident event, death, or censoring
dt.cox[, t.follow := as.numeric(
  fifelse(!is.na(event_date), event_date, end.date.observation) - baseline_date
)]

## final table in long format
dt.cox               <- dt.cox[, .(f.eid,
                                   event,
                                   t.follow,             # days from baseline to event, death, or censoring
                                   event.occurred,       # incident event indicator (post-baseline)
                                   pre.baseline.event,   # prevalent event indicator (pre-baseline)
                                   baseline_date,
                                   end.date.observation)]

## remove what is no longer needed
rm(all.events); rm(all.ids); rm(dt.post); rm(dt.pre); rm(dt.skeleton); rm(ukb.first.occurrence); rm(ukb.death); gc(reset = T)

## recode event
dt.cox[, event := as.character(event)]

## add sex
dt.cox               <- merge(dt.cox, ukb.dat[, .(f.eid, sex)])

#------------------------------#
##--    prevalent events    --##
#------------------------------#

## count prevalent events
event.count.prev     <- dt.cox[ , .(num.events.female = sum(pre.baseline.event[ sex == "Female"] == 1),
                                    num.events.male = sum(pre.baseline.event[ sex == "Male"] == 1)), 
                                by = "event"]
## create combined number
event.count.prev[, num.events := num.events.female + num.events.male]
## create sex column
event.count.prev[, sex := ifelse(num.events > 0 & num.events.male == 0, "Female",
                                 ifelse(num.events.female == 0, "Male", "Both"))]

## subset to sufficient number of events
event.count.prev     <- event.count.prev[ num.events >= 50]
## add some explanation
event.count.prev     <- merge(event.count.prev, lab.first.occurrence[, .(id.ukbb, description)],
                              by.x = "event", by.y = "id.ukbb")

## create wide data set
icd10.prev.wide      <- dcast(dt.cox[ event %in% event.count.prev$event], f.eid ~ event, value.var = "pre.baseline.event")

## add to ukb
ukb.dat              <- merge(ukb.dat, icd10.prev.wide)

## augment label set accordingly
tmp                  <- event.count.prev[, .(event, description, sex)]
names(tmp)           <- c("short_name", "label", "sex")
tmp[, type := "integer"]
tmp[, inv.transform := F]
tmp[, released := T]
tmp[, category := "Diseases"]
tmp[, miss.per := 0]

## add sex column, including for medications
lab.set[, sex := "Both"]
tmp2                 <- melt(ukb.dat[, lapply(.SD, sum, na.rm = TRUE), by = sex, .SDcols = lab.set[ category == "Drugs"]$short_name],
                             id.vars       = "sex",
                             variable.name = "medication",
                             value.name    = "n.users")
tmp2                 <- dcast(tmp2, medication ~ sex, value.var = "n.users", sep = "_")
tmp2[, n.users := Female + Male]
## add sex variable based on those
tmp2[, sex := ifelse(n.users > 0 & Male == 0, "Female",
                     ifelse(Female == 0, "Male", "Both"))]
## combine
lab.set              <- rbind(lab.set, tmp, fill = T)
## parse medications
lab.set              <- merge(lab.set, tmp2[, .(medication, sex)], by.x = "short_name", by.y = "medication",
                              all.x = T)
lab.set[, sex := ifelse(!is.na(sex.y), sex.y, sex.x)]
## delete what is no longer needed
lab.set[, sex.x := NULL]
lab.set[, sex.y := NULL]

## write to file for imputation/downstream analysis
fwrite(ukb.dat, "UKB.data.prep.ICD10.code.Cox.models.20260626.txt", sep = "\t", row.names = F, na = NA)
## corresponding labels
fwrite(lab.set, "UKB.labels.prep.ICD10.code.Cox.models.20260626.txt", sep = "\t", row.names = F, na = NA)

#------------------------------#
##--          store         --##
#------------------------------#

## count of events
event.count          <- dt.cox[ , .(num.events = sum(event.occurred[pre.baseline.event == 0])), 
                                by = "event"]
## subset to events with at least 100 counts
event.count          <- event.count[ num.events >= 100]
## add explanation
event.count          <- merge(event.count, lab.first.occurrence[ date.column == T, .(id.ukbb, description)],
                              by.x = "event", by.y = "id.ukbb")
## add ICD10 count
event.count[, icd10.code := {
  m <- regexpr("[A-Z][0-9]{2}", description)
  fifelse(m == -1L, NA_character_, substr(description, m, m + 2L))
}]

## write to file
write.table(event.count, "Event.count.incident.ICD.codes.UKB.20260616.txt", sep = "\t", row.names = F)

## reduce complexity of the outcome data
dt.cox               <- dt.cox[ event %in% event.count$event]

## write to file: by event to ease import later on
for(j in event.count$event){
  fwrite(dt.cox[ event == j], paste0(j, ".UKB.first.occurrence.parsed.20260520.txt"), sep = "\t", row.names = F, na = NA)
}

## delete what is no longer needed for now
rm(dt.cox); rm(ukb.dat); rm(icd10.prev.wide); gc(reset=T)

###############################################
####  compute partial correlation network  ####
###############################################

## file to run
write.table(c("all", "Female", "Male"), "Input.partial.correlation.txt", sep = "\t", row.names = F, col.names = F)

#--------------------------------#
##--       import results     --##
#--------------------------------#

## import results across all layers
res.pcor.all <- rbindlist(lapply(c("all", "Female", "Male"), function(x){
  ## import
  tmp <- fread(paste0("../output/partial.correlation.", x, "..txt"))
  ## add cohort
  tmp[, cohort := x]
  ## return
  return(tmp)
}))
## have a look at the results
View(res.pcor.all)

## write to file for review
res.pcor.all[ cohort == "all" & estimate > .2]

## simple sample
res.pcor.all[ sample(1:nrow(res.pcor.all), 20), ]

#----------------------------#
##-- flag significance     --#
#----------------------------#

## computed ONCE per edge, per cohort (before the undirected stacking, so the
## multiple-testing correction sees each pair exactly once)
res.pcor.all[, sig.bonf := pval < (.05 / .N),      by = cohort]
res.pcor.all[, sig.mag  := abs(estimate) >= .05]
res.pcor.all[, sig      := sig.bonf & sig.mag]                     # headline definition

#----------------------------#
##-- undirected node table --#
#----------------------------#

## each edge contributes to BOTH of its endpoints
long <- rbind(
  res.pcor.all[, .(cohort, icd10 = icd10.Var1, partner = icd10.Var2,
                   estimate, pearson, sig, sig.bonf, sig.mag)],
  res.pcor.all[, .(cohort, icd10 = icd10.Var2, partner = icd10.Var1,
                   estimate, pearson, sig, sig.bonf, sig.mag)]
)

## empty-set-safe summaries (a node may have no significant partners)
smean <- function(x) if (length(x)) mean(x)   else NA_real_
smed  <- function(x) if (length(x)) median(x) else NA_real_
smax  <- function(x) if (length(x)) max(x)    else NA_real_

#----------------------------#
##-- per-code metrics      --#
#----------------------------#

## compute per code metrics
parcor.metrics <- long[, .(
  ## --- across ALL partners (complete graph minus pruned pairs) -----------
  n.partners        = .N,                                   # degree
  mean.estimate     = mean(estimate),                       # signed
  median.estimate   = median(estimate),
  mean.abs.estimate = mean(abs(estimate)),                  # magnitude
  median.abs.est    = median(abs(estimate)),
  sd.estimate       = sd(estimate),
  min.estimate      = min(estimate),                        # strongest negative
  max.estimate      = max(estimate),                        # strongest positive
  strength          = sum(abs(estimate)),                   # weighted degree
  mean.pearson      = mean(pearson),                        # marginal reference
  ## --- significance counts ----------------------------------------------
  n.sig             = sum(sig),                             # headline (FDR & |est|>=cutoff)
  prop.sig          = mean(sig),
  n.sig.bonf        = sum(sig.bonf),
  n.mag             = sum(sig.mag),                         # |est| >= cutoff regardless of p
  n.pos.sig         = sum(sig & estimate > 0),              # significant & positive
  n.neg.sig         = sum(sig & estimate < 0),              # significant & negative (rare)
  ## --- restricted to SIGNIFICANT partners -------------------------------
  mean.estimate.sig   = smean(estimate[sig]),
  median.estimate.sig = smed(estimate[sig]),
  mean.abs.est.sig    = smean(abs(estimate[sig])),
  max.abs.est.sig     = smax(abs(estimate[sig])),
  strength.sig        = sum(abs(estimate[sig]))
), by = .(cohort, icd10)]

## replace missing values
parcor.metrics[ is.na(parcor.metrics)] <- 0

## delete what is no longer needed
rm(long); gc(reset=T)

## add some labels
lab.first.occurrence[, icd10.code := {
  m <- regexpr("[A-Z][0-9]{2}", description)
  fifelse(m == -1L, NA_character_, substr(description, m, m + 2L))
}]

## add
parcor.metrics <- merge(parcor.metrics, lab.first.occurrence[ date.column == T, .(icd10.code, description)], 
                        by.x = "icd10", by.y = "icd10.code")

## write to file
write.table(parcor.metrics, "Partial.correlation.metrics.ICD10.codes.UKB.20260617.txt", sep = "\t", row.names = F)

###############################################
####            import results             ####
###############################################

#------------------------------#
##--  file for associations --##
#------------------------------#

## define array of jobs to be run
tmp <- expand.grid(target.pop = c("all", "men", "women"), 
                   target.ehr = c("all"),
                   outc = event.count$event,
                   adj = c("+ age + sex + centre"), stringsAsFactors = F)
tmp <- as.data.table(tmp)
tmp[, adj := ifelse(target.pop != "all", gsub(" \\+ sex", "", adj), adj)]
## write to file
write.table(tmp, "Input.Cox.models.txt",quote = F, row.names = F, col.names = F, sep = "\t")

#-------------------------------#
##--      minimal model      --##
#-------------------------------#

## do in parallele
registerDoMC(10)

## import results: --> lot's missing--
res.cox.minimal <- grep("cox", dir("../output/"), value = T)

## some are missing; find out which
tmp[, file_name := paste("cox", outc, target.pop, target.ehr, gsub("\\s*\\+\\s*", ".", gsub("^\\+\\s*", "", adj)), "txt", sep = ".")]
tmp[, completed := file_name %in% res.cox.minimal]
tmp[ completed == F]
tmp             <- merge(tmp, lab.first.occurrence[, .(id.ukbb, description)],
                         by.x = "outc", by.y = "id.ukbb",
                         all.x = T)
## only curated, none ICD-10 coded, outcomes are missing; ignore for now

## import results
res.cox.minimal <- rbindlist(mclapply(res.cox.minimal, function(x){
  ## import
  tmp <- fread(paste0("../output/", x))
  ## add population
  tmp[, sex := strsplit(x, "\\.")[[1]][6]]
} , mc.cores = 10))

## adopt labels
lab.ext         <- rbind(lab.ext,
                         lab.set[ category == "Diseases"],
                         fill = T)
## adopt some columns
lab.ext[, short_name_new := ifelse(is.na(short_name_new), gsub("\\.", "_", short_name), short_name_new)]
lab.ext[, label_new := ifelse(is.na(label_new), label, label_new)]

## add explanations: drops fractured heel
res.cox.minimal <- merge(res.cox.minimal, lab.ext[, .(short_name_new, short_name, label_new, category)], 
                         by = "short_name_new")
## the outcome
res.cox.minimal <- merge(res.cox.minimal, lab.first.occurrence[, .(id.ukbb, description)],
                         by.x = "event", by.y = "id.ukbb")

## add information on sex-specificty to split, if needed
# res.cox.minimal <- merge(res.cox.minimal, event.count)

# ## drop rare events
# res.cox.minimal <- res.cox.minimal[ nevent >= 50]

## view results
View(res.cox.minimal[ `Pr(>|z|)` < .05/nrow(res.cox.minimal) & sex == "all" & short_name != "sex"])
View(res.cox.minimal[ `Pr(>|z|)` < .05/nrow(res.cox.minimal) & category == "Diseases"])

## write to file: use lenient p-value for all traits, but the one later used for evaluating improved performance
fwrite(res.cox.minimal[ `Pr(>|z|)` <  3.287959e-06], "Results.Cox.models.first.occurrence.ICD10.UKB.20260618.txt",
       sep = "\t", row.names = F, na = NA)

## flag exposures Delphi could never sensibly ingest as a clean signal
## (definitional / reverse-causal / iatrogenic), reusing the AI annotation: does only affect continuous outcomes
res.annot       <- fread("Results_Cox_reannotated_continuous_20260618.txt")
## add to res.cox.minimal
res.cox.minimal <- merge(res.cox.minimal, res.annot[, .(event, short_name_new, sex, novelty, mechanism_category,
                                                        adr_flag, icd10.code, reverse_causation_flag, annotation.source, direction)],
                         by = c("event", "short_name_new", "sex"),
                         all.x = T)

###############################################
####   add delphi predictions and gen cor. ####
###############################################

#---------------------------------#
##--       import results      --##
#---------------------------------#

## import information from Wenhuan
delphi.stats             <- fread("<path_to_file>")
## import parsed genetic correlation stats
tmp.delphi               <- fread("<path_to_file>")
## flag correlation with broad range of codes
tmp.delphi[, icd10.specific := as.integer(!grepl("^ICD10\\s+[A-Z][0-9]{2}-[A-Z]?[0-9]{2}", `Reported trait`))]

## import updated stats from Wenhuan
auc.stats                <- fread("<path_to_file>")
## get ICD-10 code
auc.stats[, icd10.code := fifelse(grepl("^[A-Z][0-9]{2}", Name), substr(Name, 1L, 3L), NA_character_)]

#---------------------------------#
##-- convert association table --##
#---------------------------------#

## create matrix accordingly
tmp.cox                  <- dcast(res.cox.minimal[ sex == "all"], event + description ~ short_name_new, value.var = "c.delta", fill = 0L)
## fill missing values
tmp.cox[ is.na(tmp.cox)] <- 0

## extract ICD-10 code from description strings
tmp.cox[, icd10.code := {
  m <- regexpr("[A-Z][0-9]{2}", description)
  fifelse(m == -1L, NA_character_, substr(description, m, m + 2L))
}]

## add to cox results as well
res.cox.minimal[, icd10.code := {
  m <- regexpr("[A-Z][0-9]{2}", description)
  fifelse(m == -1L, NA_character_, substr(description, m, m + 2L))
}]

#---------------------------------#
##--        combine data       --##
#---------------------------------#

## add to delphi.stats
auc.pred.dat             <- merge(delphi.stats, tmp.cox,
                                  by.x = "ICD-10 code", 
                                  by.y = "icd10.code") 
## looses some ICD-10 codes, either being vast, or no incident events after baseline
## n = 636 overlapping

## test for what is not overlapping
tmp.cox[ !(icd10.code %in% auc.pred.dat$`ICD-10 code`), 1:5]
## includes only manually curated outcome

## add sample size
auc.pred.dat             <- merge(auc.pred.dat, res.cox.minimal[, .(number_cases = max(nevent, na.rm = T)), by = "icd10.code"],
                                  by.x = "ICD-10 code", 
                                  by.y = "icd10.code")

## add parsed delphi stats: omit duplications due to category levels
auc.pred.dat             <- merge(auc.pred.dat, tmp.delphi[ icd10.specific == 1, .(`ICD-10 code`, mean_rg, mean_abs_rg, absolute_mean_rg)],
                                  by = "ICD-10 code", all.x = T,
                                  suffixes = c(".all", ".parsed"))

## add evaluation data set to get baseline for AUPRC
auc.pred.dat             <- merge(auc.pred.dat, auc.stats[, .(icd10.code, n_diseased, n_healthy)],
                                  by.x = "ICD-10 code", 
                                  by.y = "icd10.code")

## partial correlation metrics as well
auc.pred.dat             <- merge(auc.pred.dat, parcor.metrics[ cohort == "all", !("description")],
                                  by.x = "ICD-10 code", 
                                  by.y = "icd10", all.x = T)

## write to file
write.table(auc.pred.dat, "Data.AUC.prediction.Delphi2M.c.delta.Cox.models.20260617.txt", sep = "\t", row.names = F)

#---------------------------------#
##--    association testing    --##
#---------------------------------#

## import function to do so
source("../functions/regression_adaptive.R")

## run association testing: N.B.: adjustment for multiple testing still needs to be implemented properly
## omit brier score, since not meaningful, if not scaled
res.auc.pred             <- regression_analysis(
  dt           = auc.pred.dat,
  outcome      = c("auc", "auprc"),
  features     = c(names(tmp.cox)[-c(1:2, ncol(tmp.cox))], 
                   "mean_rg.all", "mean_abs_rg.all", "absolute_mean_rg.all",
                   "mean_rg.parsed", "mean_abs_rg.parsed", "absolute_mean_rg.parsed",
                   names(parcor.metrics)[-c(1, ncol(parcor.metrics))]),
  covariates   = ~ log10(n_diseased),
  outcome_rint = T,
  p_adjust     = "BH",
  spline_type  = "rcs",
  knots        = 3L,
  verbose      = F
)

## add label
res.auc.pred             <- merge(lab.ext[, .(short_name_new, label_new, category)], res.auc.pred, all.y = T,
                                  by.x = "short_name_new", by.y = "feature")
## rename
names(res.auc.pred)      <- gsub("short_name_new", "feature", names(res.auc.pred))

#---------------------------------#
##--    fine-map associations  --##
#---------------------------------#

## rank-based inverse-normal transform (as used elsewhere)
INT                      <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))

## strip the log10(n_diseased) gradient from each metric
## na.exclude pads residuals back to full length so the := assignment aligns
auc.pred.dat[, auc.resid       := residuals(lm(auc   ~ log10(n_diseased), data = auc.pred.dat, na.action = na.exclude))]
auc.pred.dat[, auprc.resid     := residuals(lm(auprc ~ log10(n_diseased), data = auc.pred.dat, na.action = na.exclude))]

## inverse-normal the residuals -> Gaussian response for susie (matches outcome_rint = T)
auc.pred.dat[, auc.resid.int   := INT(auc.resid)]
auc.pred.dat[, auprc.resid.int := INT(auprc.resid)]

## --> restrict features <-- ##

## features to be considered (omit drugs and diseases entirely)
feat.susie               <- unique(res.auc.pred[ !(category %in% c("Drugs", "Diseases")) & 
                                                   feature != "cohort" & !is.na(beta)]$feature)
## omit genetic correlation feature as well
feat.susie               <- feat.susie[-which(feat.susie %in% c("mean_rg.all", "mean_abs_rg.all", "absolute_mean_rg.all",
                                                                "mean_rg.parsed", "mean_abs_rg.parsed", "absolute_mean_rg.parsed"))]
## n = 276

## run susie with specified parameters
set.seed(42)
res.susie.restricted     <- lapply(c("auc.resid.int", "auprc.resid.int"), function(x){
  ## run susie: drop one outcome with no partial correlation data
  susie(as.matrix(auc.pred.dat[ !is.na(n.sig), feat.susie, with=F]),
        unlist(auc.pred.dat[ !is.na(n.sig), ..x]),
        L=10, max_iter = 1000, coverage = .99, min_abs_corr = .1)
})

## get selected credible sets
res.credible.restricted  <- lapply(res.susie.restricted, function(x){
  ## obtain credible set information
  cred.sel  <- data.table(summary(x)$vars)
  cred.sel[, short_name_new := feat.susie[variable]]
  cred.sel  <- cred.sel[ order(cs)]
  ## add label
  cred.sel  <- merge(cred.sel, rbind(lab.ext[, .(short_name_new, label_new, category)],
                                     data.table(short_name_new = names(parcor.metrics)[-c(1:2, ncol(parcor.metrics))], 
                                                label_new = names(parcor.metrics)[-c(1:2, ncol(parcor.metrics))],
                                                category = "partial_correlation_network")),
                     by = "short_name_new")
  ## order again
  cred.sel  <- cred.sel[ order(cs, -variable_prob)]
  ## return relevant results
  return(cred.sel[ cs > 0])
})

## write results to file
res.credible.restricted[[1]]$outcome <-"auc"
res.credible.restricted[[2]]$outcome <-"auprc"
write.table(rbindlist(res.credible.restricted), "Results.AUC.prediction.fine.mapping.restricted.set.of.features.Delphi2M.20260617.txt", sep = "\t", row.names = F)

## combine with association statistics
jj  <- res.susie.restricted[[1]]$pip
## add to tmp results
tmp <- res.auc.pred[ outcome == "auc" & term == "linear" & feature %in% names(jj)]
## add PIP
tmp[, pip := jj[feature]]
## add credible set
tmp <- merge(tmp, res.credible.restricted[[1]][, .(short_name_new, cs)],
             by.x = "feature", by.y = "short_name_new", all.x = T)
## write to file
fwrite(tmp, "Results.fine.mapping.AUC.associations.20260630.txt", sep = "\t", row.names = F, na = NA)

#---------------------------------#
##--   association testing II  --##
#---------------------------------#

## run association testing with additional adjustment: genetic only
res.auc.pred.adjusted <- regression_analysis(
  dt           = auc.pred.dat,
  outcome      = c("auc", "auprc"),
  features     = c("mean_rg.all", "mean_abs_rg.all", "absolute_mean_rg.all",
                   "mean_rg.parsed", "mean_abs_rg.parsed", "absolute_mean_rg.parsed"),
  covariates   = ~ log10(n_diseased) + n.sig.bonf + Long_standing_illness__disability_or_infirmityYes + median.estimate + Falls_in_the_last_yearMore_than_one_fall,
  outcome_rint = T,
  p_adjust     = "BH",
  spline_type  = "rcs",
  knots        = 3L,
  verbose      = F
)

## look into results
View(res.auc.pred.adjusted[ term == "overall"])
## genetic associations persist

## boxplot on whether or not genetic associations did exist
boxplot(auc ~ is.na(mean_rg.parsed), auc.pred.dat)
t.test(auc ~ is.na(mean_rg.parsed), auc.pred.dat) ## p-value = 0.002939
wilcox.test(auc ~ is.na(mean_rg.parsed), auc.pred.dat) ## p-value = 0.003047

###############################################
####       embedding on risk factors       ####
###############################################

## import reconstructed exposures
recon.exposure <- fread("../output/Reconstruction_R2_exposures_by_sex_Delphi_embeddings.20260617.txt")
## filter to variables cleanly predictable from OLS
recon.exposure <- recon.exposure[ !(short_name %in% c("age", "centre")) & category != "Technical"]

## how many exposure
uniqueN(recon.exposure$short_name)
## N = 117

###############################################
####   compute delta-AUC to match Delphi   ####
###############################################

#---------------------------------#
##--   select the shortlist    --##
#---------------------------------#

## thresholds are deliberately exposed -- edit to taste
p.thresh            <- 3.287959e-06   ## study-wide significance used throughout
delta.thresh        <- 0.005    ## minimal gain over Delphi worth confirming

## keep significant, positive-gain rows: do not subset for delphi
res.cox.delphi      <- res.cox.minimal[ short_name %in% recon.exposure$short_name & 
                                          !is.na(c.delta.delphi) &
                                          # c.delta.delphi > delta.thresh &
                                          `Pr(>|z|)` < p.thresh & nevent >= 100]

## the delta-AUC function takes a single continuous/binary exposure -> keep those
res.cox.delphi      <- merge(res.cox.delphi, lab.ext[, .(short_name_new, type, inv.transform)], by = c("short_name_new"), all.x = T)

#---------------------------------#
##--   map outcome -> token    --##
#---------------------------------#

## import information from Wenhuan: enable mapping
delphi.stats         <- fread("<path_to_file>")
## add ICD-10 code
delphi.stats[, icd10.code := {
  m <- regexpr("[A-Z][0-9]{2}", Name)
  fifelse(m == -1L, NA_character_, substr(Name, m, m + 2L))
}]

## add relevant information
res.cox.delphi       <- merge(res.cox.delphi, delphi.stats[, .(icd10.code, Name)], by = "icd10.code")

## token = the column name under which the python LP is stored for this outcome
res.cox.delphi[, token := Name]

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
## n = 44,548

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
## n = 1260

#---------------------------------#
##--      import results       --##
#---------------------------------#

## do in parallel
registerDoMC(10)

## import results: --> lot's missing--
res.auc.delta   <- grep("delta_auc\\.", dir("../output/"), value = T)

## import results
res.auc.delta   <- rbindlist(mclapply(res.auc.delta, function(x){
  ## import
  tmp <- fread(paste0("../output/", x))
} , mc.cores = 10))

## add labels
res.auc.delta   <- merge(res.auc.delta, lab.first.occurrence[ date.column == T, .(id.ukbb, description, icd10.code)],
                         by.x = "event", by.y = "id.ukbb") 
## same for exposure
res.auc.delta   <- merge(res.auc.delta, lab.ext[, .(short_name_new, category, label_new)],
                         by = "short_name_new")

## write results to file
write.table(res.auc.delta, "Results.delta.AUC.Delphi.2M.Cox.model.informed.20260622.txt", sep = "\t", row.names = F)

###############################################
####           protein improvement         ####
###############################################

#-------------------------------#
##--  import and parse data  --##
#-------------------------------#

## import protein data
ukb.prot                <- fread("<path_to_file>")
## import label
lab.prot                <- fread("<path_to_file>")

#------------------------#
##--    do some QC    --##
#------------------------#

## count missing values by participant (start here, since one entire batch of proteins is missing for a subset of individuals)
eid.miss.prot           <- apply(ukb.prot[, lab.prot$id, with=F], 1, function(x) sum(is.na(x))/nrow(lab.prot))
eid.miss.prot           <- data.frame(eid=ukb.prot$eid, miss.per = eid.miss.prot)
## drop participants with >50% missing values
ukb.prot                <- ukb.prot[ eid %in% subset(eid.miss.prot, miss.per <= .5)$eid]

## same for proteins
lab.prot[, miss.per := sapply(id, function(x) nrow(ukb.prot[is.na(eval(as.name(x)))]))/nrow(ukb.prot)*100]
## drop proteins with >20% missing values
lab.prot                <- lab.prot[ miss.per <= 20 ]
## n = 2919 proteins

## apply normalization
ukb.prot                <- as.data.frame(ukb.prot)
ukb.prot[, lab.prot$id] <- apply(ukb.prot[, lab.prot$id], 2, function(x){
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
})

#----------------------------#
##--        export        --##
#----------------------------#

## export data to run association testing
fwrite(ukb.prot, "UKB.Olink.delphi.2M.LP.20260618.txt", sep="\t", row.names=F, na = NA)
write.table(lab.prot, "Olink.proteins.delphi.2M.LP.20260618.txt", sep="\t", row.names=F)
## delete what is no longer needed
rm(ukb.prot); gc(reset=T)

#----------------------------#
##--    import results    --##
#----------------------------#

## look at results:
res.cox.proteins <- grep("cox.prot", dir("../output/"), value = T)
res.cox.proteins <- rbindlist(mclapply(res.cox.proteins, function(x){
  ## import
  tmp <- fread(paste0("../output/", x))
  ## add population
  tmp[, sex := strsplit(x, "\\.")[[1]][7]]
} , mc.cores = 10))

## the outcome
res.cox.proteins <- merge(res.cox.proteins, lab.first.occurrence[, .(id.ukbb, description)],
                          by.x = "event", by.y = "id.ukbb")
## order by delta
res.cox.proteins <- res.cox.proteins[ order(-c.delta)]

## omit those with too few cases to sensibly estimates anything
res.cox.proteins <- res.cox.proteins[ nevent >= 100]

#----------------------------#
##--      prep for AUC    --##
#----------------------------#

## create new map
lab.prot[, short_name_new := id]
id.map     <- lab.prot[, .(short_name_new, short_name = id)]

## subset of interest
short.full <- res.cox.proteins[ !is.na(`Pr(>|z|)`) & `Pr(>|z|)` < p.thresh,
                                .(event, short_name_new, sex)]
short.full <- merge(short.full, id.map, by = "short_name_new", all.x = T)

## (a) shortlist consumed by the runner: event | short_name | sex  (header required)
shortlist  <- unique(short.full[, .(event, short_name, sex)])
fwrite(shortlist, "protein.shortlist.txt", sep = "\t", row.names = F, quote = F, na = NA)

## (b) array list for the submit script: target_pop  target_ehr  outc  (NO header, awk reads by position)
short.full[, target_ehr := "all"]
arr        <- unique(short.full[, .(target_pop = sex, target_ehr, outc = event)])
fwrite(arr, "Input.delta.auc.proteins.txt", sep = "\t", col.names = F, row.names = F, quote = F, na = NA)

cat("shortlisted protein-disease pairs:", nrow(shortlist),
    "| proteins:", uniqueN(shortlist$short_name),
    "| outcomes:", uniqueN(shortlist$event), "\n")
## shortlisted protein-disease pairs: 68548 | proteins: 2029 | outcomes: 323 
cat("array rows (set --array=1-", nrow(arr), " in 12_submit_delta_auc_proteins.sh)\n", sep = "")
## array rows (set --array=1-780 in 12_submit_delta_auc_proteins.sh)

#----------------------------#
##--     import results   --##
#----------------------------#

## do in parallele
registerDoMC(10)

## import results: --> lot's missing--
res.auc.proteins   <- grep("delta_auc_prot", dir("../output/"), value = T)

## import results
res.auc.proteins   <- rbindlist(mclapply(res.auc.proteins, function(x){
  ## import
  tmp <- fread(paste0("../output/", x))
} , mc.cores = 10))

## add labels
res.auc.proteins   <- merge(res.auc.proteins, lab.first.occurrence[ date.column == T, .(id.ukbb, description, icd10.code)],
                            by.x = "event", by.y = "id.ukbb") 

## drop very small ones
res.auc.proteins   <- res.auc.proteins[ n.case >= 25]

## order
res.auc.proteins   <- res.auc.proteins[ order(-delta.auc.cv)]

## write results to file
write.table(res.auc.proteins, "Results.delta.AUC.Delphi.2M.Cox.model.informed.Olink.proteins.20260622.txt", sep = "\t", row.names = F)

#----------------------------#
##--    reconstruction    --##
#----------------------------#

## import 
recon.proteins <- fread("../output/Reconstruction_R2_proteins_by_sex_Delphi_embeddings.20260619.txt")

###############################################
####      reporting for the manuscript     ####
###############################################

#--------------------------------#
##--    predicting the AUC    --##
#--------------------------------#

## how many diseases in the partial correlation network
nrow(parcor.metrics[ cohort == "all"])
## n = 780

## risk factors for delta c-index (restricted set)
uniqueN(lab.ext[ !(category %in% c("Drugs", "Diseases")) & 
                   short_name_new %in% res.auc.pred$feature &
                   !(short_name %in% c("age", "sex", "centre"))]$short_name)
table(lab.ext[ !(category %in% c("Drugs", "Diseases")) & 
                 short_name_new %in% res.auc.pred$feature &
                 !(short_name %in% c("age", "sex", "centre"))]$category)

## credible sets
res.credible.restricted

## specific estimates
res.auc.pred[ term == "linear" & feature %in% res.credible.restricted[[1]]$short_name_new & outcome == "auc",
              .(feature, beta, p_value)]
## correlation among predictors
feature <- res.auc.pred[ term == "linear" & !is.na(beta)]$feature
cor.tmp <- cor(as.matrix(auc.pred.dat[, ..feature]), use = "p")
# convert the matrix into a long pairwise table
cor.tmp <- as.data.table(as.table(cor.tmp))
setnames(cor.tmp, c("var1", "var2", "correlation"))
# keep only unique pairs (drop self-correlations and mirror duplicates)
cor.tmp <- cor.tmp[var1 < var2][order(-abs(correlation))]
## look into those of interest
View(cor.tmp[ var1 == "Long_standing_illness__disability_or_infirmityPrefer_not_to_answer" | var2 == "Long_standing_illness__disability_or_infirmityPrefer_not_to_answer"])

## selected ones
res.auc.pred[ term == "linear" & feature %in% c("log_crp_cleaned", "chol_cleaned", "whr") & outcome == "auc",
              .(feature, beta, p_value)]

## AUC results after adjustment
res.auc.pred.adjusted[ term == "linear" & outcome == "auc", .(feature, beta, p_value)]

#--------------------------------#
##--     improving the AUC    --##
#--------------------------------#

## things reconstructed
uniqueN(recon.exposure$short_name)
## n = 117

## build a common data frame for reporting
cox.auc.reporting <- merge(res.cox.minimal[ short_name %in% recon.exposure$short_name & `Pr(>|z|)` < p.thresh & nevent >= 100 &
                                              novelty != "Artefact" & reverse_causation_flag == "FALSE" &
                                              adr_flag == "FALSE" & !grepl("efinition", mechanism_category)],
                           res.auc.delta,
                           by = intersect(names(res.cox.minimal), names(res.auc.delta)), 
                           all.x = T)

## write
write.table(cox.auc.reporting, "../Summary.Cox.AUC.Delphi.2M.biomarker.comparison.20260622.txt", sep = "\t", row.names = F)

## --> ability to reconstruct <-- ##

## reconstruction summary by sex
recon.exposure[, .(median_r2 = median(recon.r2),
                   p25_r2 = quantile(recon.r2, .25),
                   p75_r2 = quantile(recon.r2, .75))
               , by = "sex"]
#         sex  median_r2     p25_r2    p75_r2
#      <char>      <num>      <num>     <num>
#   1:    all 0.09149446 0.03982120 0.3551401
#   2:    men 0.06143545 0.02849070 0.1221995
#   3:  women 0.06243574 0.02591665 0.1544749

## how many reached >20%
uniqueN(recon.exposure[ sex != "all" & recon.r2 > .2]$short_name)/uniqueN(recon.exposure$short_name)
# 0.2136752

## plasma proteins
recon.proteins[, .(median_r2 = median(recon.r2),
                   p25_r2 = quantile(recon.r2, .25),
                   p75_r2 = quantile(recon.r2, .75))
               , by = "sex"]
#         sex  median_r2     p25_r2     p75_r2
#      <char>      <num>      <num>      <num>
#   1:    all 0.02631116 0.01064608 0.07429466
#   2:    men 0.02155471 0.01264014 0.05554138
#   3:  women 0.02279142 0.01067687 0.06139833

## reconstruction summary by sex and group
recon.exposure[ category == "Body composition", .(median_r2 = median(recon.r2),
                                                  p25_r2 = quantile(recon.r2, .25),
                                                  p75_r2 = quantile(recon.r2, .75))
                , by = "sex"]
#         sex median_r2    p25_r2    p75_r2
#      <char>     <num>     <num>     <num>
#   1:    all 0.5295595 0.4728285 0.6436127
#   2:    men 0.3226034 0.2346098 0.3776787
#   3:  women 0.3738817 0.2369022 0.4564600

## reconstruction blood lipids
recon.exposure[ short_name %in% c("apoa_cleaned", "chol_cleaned", "hdl_chol_cleaned", "ldl_chol_cleaned", "log_lipo_a_cleaned",
                                  "log_total_tg_cleaned"), 
                .(median_r2 = median(recon.r2),
                  p25_r2    = quantile(recon.r2, .25),
                  p75_r2    = quantile(recon.r2, .75))
                , by = "sex"]
#         sex median_r2     p25_r2    p75_r2
#      <char>     <num>      <num>     <num>
#   1:    all 0.1588418 0.14463746 0.2104923
#   2:    men 0.1016931 0.07496379 0.1458779
#   3:  women 0.1335494 0.10916275 0.1455188

## --> reconstruction vs headroom (per pair, predictive, case-floored) <-- ##
source("../functions/reconstruction_vs_headroom_perpair.R")

## phenotypes: cox.auc.reporting carries short_name, delta.auc.cv, n.case (+ gain.cv after the 07 rerun)
recon.vs.headroom(dt.pairs = cox.auc.reporting[ gain.mode == "stratified" & n.strata.cv >= 2], dt.recon = recon.exposure,
                  id.col = "short_name", case.floor = 100L, group = "exposures")

## proteins: res.auc.proteins keys on 'exposure'; recon.proteins keys on 'short_name'
recon.vs.headroom(dt.pairs = res.auc.proteins, dt.recon = recon.proteins,
                  id.col = "exposure", case.floor = 100L, group = "proteins")

## --> ability to improve <-- ##

## number of diseases improving: biomarkers
median(cox.auc.reporting[ sex != "all" & !is.na(gain.cv), gain.cv], na.rm = F)
nrow(cox.auc.reporting[ sex != "all" & !is.na(gain.cv)])

## --> what is worth adding to Delphi-2M: points 1-3, biomarkers AND proteins <-- ##

## exposure / short_name_new / short_name.
source("../functions/improvement_report.R")
res.improve <- improvement.report(
  sets = list(
    biomarker = list(auc = res.auc.delta,    cox = res.cox.minimal,     recon = recon.exposure,
                     auc.id = "short_name",   cox.id = "short_name",     recon.id = "short_name",
                     highlight = c("log_tbil_cleaned", "ca_cleaned", "log_ggt_cleaned",
                                   "hba1c_cleaned", "log_cysc_cleaned", "bmi")),
    protein   = list(auc = res.auc.proteins, cox = res.cox.proteins,    recon = recon.proteins,
                     auc.id = "exposure",     cox.id = "short_name_new", recon.id = "short_name",
                     highlight = c("tshb", "igfbp6", "nbl1", "gdf15"))
  ),
  p.thresh = p.thresh, case.floor = 100L)

###############################################
####       figures for the manuscript      ####
###############################################

#--------------------------------#
##--       AUC prediction     --##
#--------------------------------#

## import function to do so
source("../functions/plot_regression_scatter.R")

## open pdf
pdf("../graphics/Summary.AUC.prediction.Delphi-2M.20260618.pdf", width = 6.3, height = 3.15)
## 2 x 4 grid
par(mfrow = c(2,4), mar=c(2.5,2.5,.5,.5), lwd=.5, cex.axis = .7, cex.lab = .7, tck = -.01, mgp = c(.8,.2,0))

## --> four top hits from SuSie and selected <-- ##

## what to plot
tmp.plt <- data.table(variable = c("n.sig.bonf", "median.estimate", "Long_standing_illness__disability_or_infirmityYes", "Falls_in_the_last_yearMore_than_one_fall",
                                   "log_crp_cleaned", "chol_cleaned", "whr"),
                      label = c("PCorNet - Sign. partners", "PCorNet - median r", "Long-standing illness", "Falls in in last year",
                                "C-reative protein", "Total cholesterol", "Waist-to-hip ratio"))

for(j in 1:nrow(tmp.plt)){
  ## example
  plot_regression_scatter(
    feature      = tmp.plt$variable[j],
    outcome      = "auc",
    results      = res.auc.pred,
    model_choice = "linear",
    dt           = auc.pred.dat,
    ylim = c(-3.5,3.5),
    add_stats = F,
    xlab = paste(tmp.plt$label[j], "\n(RINT)"),
    ylab = "Delphi-2M AUC (RINT)",
    covariate_handling = "marginal",
    covariates   = ~ log10(n_diseased),
    outcome_rint = TRUE,
    main = NA
  )
}

## --> genetics adjusted for other variables <-- ##

## example
plot_regression_scatter(
  feature      = "mean_rg.parsed",
  outcome      = "auc",
  results      = res.auc.pred.adjusted,
  model_choice = "linear",
  dt           = auc.pred.dat,
  ylim = c(-3.5,3.5),
  add_stats = F,
  xlab = "Average genetic correlation\n(RINT)",
  ylab = "Delphi-2M AUC (RINT)",
  covariate_handling = "marginal",
  covariates   = ~ log10(n_diseased) + n.sig.bonf + Long_standing_illness__disability_or_infirmityYes + median.estimate + Falls_in_the_last_yearMore_than_one_fall,
  outcome_rint = TRUE,
  main = NA
)

## close device
dev.off()

#----------------------------------#
##--  enhanced predicitive cap. --##
#----------------------------------#

## open pdf
pdf("../graphics/Summary.AUC.improvement.Delphi-2M.20260623.pdf", width = 6.3, height = 6.3)
## adopt paramaters
par(mar=c(2.5,2.5,.5,.5), lwd=.5, cex.axis = .7, cex.lab = .7, tck = -.01, mgp = c(.8,.2,0), yaxs = "r")

## more plots on a figure
layout(matrix(c(1:3,4,4,4), 3, 2), widths = c(.4,.6))

## --> recontruction by category <-- ##

## empty plot
plot(c(0,1), c(0,.6), type = "n", xaxt = "n", yaxt = "n", 
     xlab = "Ordered risk factors", 
     ylab = "Reconstruction by embeddings [R2] (within sex)")
## add axis
axis(2, lwd=.5)

## biomarker: create temporary file, order, and add dots, annotate top four
tmp <- res.improve$reconstruction[ feature.type == "biomarker"]
tmp <- tmp[order(recon.r2)]
## add position in the plot
tmp[, plt.pos := seq(0, 1, length.out = nrow(tmp))]
## place equidistant
points(tmp$plt.pos, tmp$recon.r2, cex=.5, pch = 21,
       lwd = .3, col = adjustcolor("#E8788A", alpha.f = .5))
## add top three
tmp <- tail(tmp, 3)
text(tmp$plt.pos, tmp$recon.r2, cex=.7, pos=2, 
     labels = sapply(tmp$feature, function(x) lab.set[ short_name == x]$label),
     xpd = NA)

## Proteins: create temporary file, order, and add dots, annotate top four
tmp <- res.improve$reconstruction[ feature.type == "protein"]
tmp <- tmp[order(recon.r2)]
## add position in the plot
tmp[, plt.pos := seq(0, 1, length.out = nrow(tmp))]
## place equidistant
points(tmp$plt.pos, tmp$recon.r2, cex=.5, pch = 21,
       lwd = .3, col = adjustcolor("#5B8DB8", alpha.f = .5))
## add top three
tmp <- tail(tmp, 3)
text(tmp$plt.pos, tmp$recon.r2, cex=.7, pos=4, 
     labels = toupper(tmp$feature),
     xpd = NA)


## add legend
legend("topleft", bty = "n", cex = .7, lty = 0,
       pch = 22, pt.lwd=.5, pt.cex = 1,
       ncol=2, legend = c("Biomarker", "Proteins"),
       pt.bg = c("#E8788A", "#5B8DB8"))

## --> Delta AUC above Delphi <-- ##

## empty plot
plot(c(-0.02,.2), c(0,14000), type = "n", xaxt = "n", yaxt = "n", 
     xlab = "Delta AUC over Delphi-2M (per pair)", 
     ylab = "(Bio)marker - disease pairs (sex-stratified)")
## add axis
axis(2, lwd = .5); axis(1, lwd = .5)

## add histogram
hist(res.improve$pairs$delta.auc.cv, add = T, breaks = 200,
     lwd=.3, border = "white", col = "grey80")
## add median
abline(v=median(res.improve$pairs$delta.auc.cv), lwd = 1, col = "#C44E52")
## get plotting coordinates
pm <- par("usr")
## add AUC≥0.02 and text
abline(v=0.02, lwd=.5, lty=2)
text(0.02, pm[4]-(pm[4]-pm[3])*.1, pos = 4, cex = .7,
     labels = paste0(sprintf("%.2f", nrow(res.improve$pairs[ delta.auc.cv >= 0.02])/nrow(res.improve$pairs)*100), 
                     "% improve >=0.02"))

## --> Absorption and pleiotropy <-- ##

## aggregate per-marker pleiotropy / absorption summary (panel c)
M          <- res.improve$features[, .(pleiotropy = mean(pleiotropy, na.rm = T), gain = mean(gain.med, na.rm = T),
                                       residual = mean(residual.med, na.rm = T)), by = .(feature, feature.type)]

## prepare for LOESS fit
d          <- M[is.finite(pleiotropy) & pleiotropy > 0 & is.finite(gain) & is.finite(residual)][order(pleiotropy)]
lx         <- log10(d$pleiotropy)
fg         <- lowess(lx, d$gain, f = .6); fr <- lowess(lx, d$residual, f = .6)
fr.y       <- approx(fr$x, fr$y, xout = fg$x, rule = 2)$y
xg         <- 10^fg$x

## empty plot
plot(NA, log = "x", xlim = range(d$pleiotropy), ylim = c(0, max(fg$y) * 1.5),
     xlab = "Number of diseases associated (Cox-model)", 
     ylab = "AUC gain", main = "", xaxt = "n", yaxt = "n")
## add axis
axis(1, lwd=.5); axis(2, lwd = .5)

## polygon describing what is absorbed
polygon(c(xg, rev(xg)), c(fg$y, rev(fr.y)), col = adjustcolor("#C44E52", 0.16), border = NA)

## add individual estimates
lines(xg, fg$y, col = "#4C72B0", lwd = 1)
lines(xg, fr.y, col = "#C44E52", lwd = 1)

## add legend
legend("topleft", legend = c("Average gain over age + sex", "Average residual over Delphi-2M"), 
       lwd = 1, lty = 1,
       col = c("#4C72B0", "#C44E52"), bty = "n", cex = .7)

## annotate
text(10^quantile(lx, 0.5), max(fg$y) * 0.42, "absorbed by\nDelphi-2M", col = "#9A3B40", cex = .7, font = 3)

## add GDF-15 and CRP
ann        <- M[toupper(feature) %in% c("GDF15", "LOG_CRP_CLEANED")]
lab.map    <- c(GDF15 = "GDF-15", LOG_CRP_CLEANED = "CRP")
## add differences
arrows(ann$pleiotropy, ann$gain, ann$pleiotropy, ann$residual, length = .05, lwd=.7, col = "grey80")
points(ann$pleiotropy, ann$gain, pch = 21, bg = "#DD8452", col = "grey80", cex = .6, lwd=.3)
points(ann$pleiotropy, ann$residual, pch = 21, bg = "#DD8452", col = "black", cex = .6, lwd=.3)
text(ann$pleiotropy, ann$gain, lab.map[toupper(ann$feature)], pos = 2, cex = .7, col = "grey20")

## --> top remaining (AUC≥0.05) <-- ##

par(mar=c(2.5,20,.5,.5), yaxs = "i")

## what to plot
tmp <- unique(res.improve$disease.best[ best.delta >= .05, .(event, best.marker, marker.class, class.label, sex, best.delta, icd)])
## refine by taking whatever is largest across sexes
tmp <- tmp[ order(event, best.delta)]
tmp[, ind := 1:.N, by = c("event")]
tmp <- tmp[ ind == 1]
## define plotting order
tmp <- tmp[ order(marker.class, -best.delta)]
tmp[, plt.srt := 1:.N]
## get statistics across sexes (if available)
tmp <- merge(tmp[, .(event, best.marker, marker.class, class.label, plt.srt)],
             res.improve$pairs,
             by.x = c("event", "best.marker"),
             by.y = c("event", "feature"),
             all.x = T)
## some lack case numbers in the opposite sex

## empty plot
plot(NA, xlim = range(tmp[, .(lci.cv, uci.cv)]), ylim = rev(c(.5, max(tmp$plt.srt)+.5)),
     xlab = "Delta AUC over Delphi-2M (best marker per disease)", 
     ylab = "", xaxt = "n", yaxt = "n")
## add axis
axis(1, lwd=.5)

## plotting coordinates
pm <- par("usr")

# ## divide by marker class
# rect(pm[1], tmp[, .(min.start = min(plt.srt)), by = "marker.class"]$min.start-.5, 
#      pm[2], tmp[, .(max.start = max(plt.srt)), by = "marker.class"]$max.start+.5,
#      border = NA, col = c("white", "grey90"))
## divide by marker class
rect(pm[1], 1:max(tmp$plt.srt)-.5, 
     pm[2], 1:max(tmp$plt.srt)+.5,
     border = NA, col = c("white", "grey90"))


## add zero crossing
abline(v=0, lwd=.5, lty=2)

## add confidence intervals: women, "Female" = "#D95F02"
arrows(tmp[ sex == "women"]$lci.cv, tmp[ sex == "women"]$plt.srt-.25,
       tmp[ sex == "women"]$uci.cv, tmp[ sex == "women"]$plt.srt-.25,
       length = 0, lwd=.5, col = adjustcolor("#D95F02", .5))
## add point estimates: women
points(tmp[ sex == "women"]$delta.auc.cv, tmp[ sex == "women"]$plt.srt-.25,
       pch = 22, cex=.6, lwd=.5, bg = "#D95F02")

## add confidence intervals: men, "Male" = "#1B9E77"
arrows(tmp[ sex == "men"]$lci.cv, tmp[ sex == "men"]$plt.srt+.25,
       tmp[ sex == "men"]$uci.cv, tmp[ sex == "men"]$plt.srt+.25,
       length = 0, lwd=.5, col = adjustcolor("#1B9E77", .5))
## add point estimates: women
points(tmp[ sex == "men"]$delta.auc.cv, tmp[ sex == "men"]$plt.srt+.25,
       pch = 22, cex=.6, lwd=.5, bg = "#1B9E77")

## add labels
tmp <- unique(tmp[, .(best.marker, description, plt.srt, class.label)])
text(pm[1], tmp$plt.srt, cex=.6, pos = 2, xpd=NA,
     labels = paste(sapply(tmp$best.marker, function(x){
       if(x %in% lab.set$short_name){
         gsub("Forced expiratory volume in 1-second (FEV1), Best measure", "FEV1 (best)", lab.set[ x == short_name]$label, fixed = T)
       }else{
         lab.prot[ id == x]$Assay
       }
     }), 
     gsub("Date |first reported ", "", tmp$description), sep = "-")) 

## add categories
arrows(pm[2], tmp[, .(min.pos = min(plt.srt)), by = "class.label"]$min.pos-.5, 
       pm[1]-(pm[2]-pm[1])*2, tmp[, .(min.pos = min(plt.srt)), by = "class.label"]$min.pos-.5,
       lwd=.5, length=0,
       xpd = NA)

## add legend
legend("bottomright", bty = "n", cex = .7, lty = 0,
       pch = 22, pt.lwd=.5, pt.cex = 1,
       ncol=1, legend = c("Women", "Men"),
       pt.bg = c("#D95F02", "#1B9E77"))

## close device
dev.off()