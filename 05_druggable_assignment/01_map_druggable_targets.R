#######################################################
#### Test for possible drug targets of embQTLs     ####
#### Maik Pietzner                      25/03/2026 ####
#######################################################

rm(list=ls())
setwd("<path_to_file>")
options(stringsAsFactors = F)
load(".RData")

## --> packages needed <-- ##
require(data.table)
require(doMC)
require(arrow)

#########################################
####      import relevant files      ####
#########################################

## --> association results <-- ##

## import lead credible set variants
res.credible        <- fread("../../03_credible_set_variants/input/Lead.credible.set.variants.incl.MHC.Delphi.embeddings.UKB.GWAS.plus.catalog.gene.annotation.20260421.txt")

## --> gene annotations <-- ##

## import information on human genes
human.genes     <- rtracklayer::readGFF("../../03_credible_set_variants/input/Homo_sapiens.GRCh37.87.gtf.gz", filter = list(type=c("gene")))
human.genes     <- as.data.table(human.genes)
## parse some chromosome names: careful chromosome is stored as factor!!
human.genes[, CHROM := ifelse(seqid == "X", 23, as.numeric(as.character(seqid)))] ## assigns some chromosomes wrongly!!
human.genes     <- human.genes[ !is.na(gene_name)]
## drop some unspecific genes
sort(table(human.genes$gene_biotype))
human.genes     <- human.genes[ gene_biotype %in% c("protein_coding", "processed_transcript", "IG_J_gene", "IG_C_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene")]

#########################################
####   import OT drug information    ####
#########################################

#-------------------------#
##--  drug target info --##
#-------------------------#

## path to results
ot.drugs   <- dir("<path_to_file>")
ot.drugs   <- grep("part", ot.drugs, value=T)
ot.drugs   <- lapply(ot.drugs, function(x) read_parquet(paste0("<path_to_file>", x)))
ot.drugs   <- lapply(ot.drugs, as.data.table)
ot.drugs   <- rbindlist(ot.drugs)
## transform some columns to make life easier
ot.drugs[, synonyms := sapply(synonyms, function(x) paste(x, collapse = "|"))]
ot.drugs[, crossReferences := sapply(crossReferences, function(x) paste(x, collapse = "|"))]
ot.drugs[, tradeNames := sapply(tradeNames, function(x) paste(x, collapse = "|"))]
ot.drugs[, childChemblIds := sapply(childChemblIds, function(x) paste(x, collapse = "|"))]
ot.drugs[, linkedDiseases.rows := sapply(linkedDiseases.rows, function(x) paste(x, collapse = "|"))]
ot.drugs[, linkedTargets.rows := sapply(linkedTargets.rows, function(x) paste(x, collapse = "|"))]

## expand by target id
ot.drugs   <- ot.drugs[ linkedTargets.rows != ""]
## drop drugs with missing information on clinical phase
ot.drugs   <- ot.drugs[ !is.na(maximumClinicalTrialPhase)]
ot.drugs   <- rbindlist(lapply(1:nrow(ot.drugs), function(x) return(data.table(ot.drugs[x,], ensembl_gene_id = strsplit(ot.drugs$linkedTargets.rows[x], "\\|")[[1]]))))

## N.B.: linkedDiseases.rows contains mapping to EFO terms

## reduce to approved medications
ot.drugs   <- ot.drugs[ maximumClinicalTrialPhase == 4 & hasBeenWithdrawn == F]
## how many unique drugs
uniqueN(ot.drugs$name)
## n = 2370

## add gene names to map to gene annotations
ot.drugs   <- merge(ot.drugs, human.genes, by.x = "ensembl_gene_id", by.y = "gene_id")
## how many unique drugs
uniqueN(ot.drugs$name)
## n = 2369

##########################################
####    coloc for druggable targets   ####
##########################################

#-----------------------------------#
##-- prepare regions for testing --##
#-----------------------------------#

## import drug-target mapping to GWAS endpoints
gwas.drugs   <- fread("<path_to_file>")

## add positions in build38 to credible set variants
snp.mapping  <- fread("<path_to_file>")
## create comparable identifier
res.credible[, marker_name_hg19 := paste(ifelse(CHROM == 23, "X", CHROM), GENPOS, pmin(ALLELE0, ALLELE1), pmax(ALLELE0, ALLELE1), sep = "_")]
## add
res.credible <- merge(res.credible, snp.mapping[, .(marker_name_hg19, pos_hg38)], by = "marker_name_hg19", all.x = T)
## delete what is no longer needed
rm(snp.mapping); gc(reset=T)

## --> define druggable targets to test <-- ##

## helper function
expand_by_vector <- function(dt, vec, col.name) {
  ## replicate dt once per entry in vec and tag each copy
  res <- rbindlist(lapply(vec, function(x) {
    tmp           <- copy(dt)
    tmp[, (col.name) := x]
    tmp
  }))
  return(res)
}

## run for all embedding x file x region combinations
res.test         <- expand_by_vector(unique(gwas.drugs[, .(CHR_numeric, region.start.third.merge, region.end.third.merge, file)]), vec = paste0("emb", 1:120), col.name = "embedding")

## write a file set to implement colocalisation analysis: N.B.: coordinates are in build 38(!)
write.table(res.test, "input.drug.colocalisation.txt", sep = "\t", row.names = F, quote = F)

#-----------------------------------#
##--        collate results      --##
#-----------------------------------#

## import results: careful, some results are to be dropped due to poor drug - indication mapping
res.coloc <- dir("../output/")
res.coloc <- rbindlist(lapply(res.coloc, function(x) fread(paste0("../output/", x))), fill = T)

#-----------------------------------#
##-- merge with drug information --##
#-----------------------------------#

## add drug information
res.drugs  <- merge(res.coloc, unique(gwas.drugs[, !c("pos.hg38", "nsnps", "CHR_ENSEMBL", "BP", "A1", "A2", "ID_UCSC", "FRQ", "P", "N", "BETA", "SE")]),
                    by.x = c("disease", "CHR_ENSEMBL", "region_start.hg38", "region_end.hg38"),
                    by.y = c("file", "CHR_numeric", "region.start.third.merge", "region.end.third.merge"),
                    allow.cartesian = T)
## subset to strong findings
res.drugs[, emb.drug.region.id := paste(embedding, disease, CHR_ENSEMBL, region_start.hg38, region_end.hg38, sep = "$")]
## get findings meeting threshold
jj         <- unique(res.drugs[PP.H4.abf > .8 & R2.lead.embQTL > .8 & type.snp == "regional.lead.drug"]$emb.drug.region.id)
res.drugs  <- res.drugs[ emb.drug.region.id %in% jj]

## how many drugs
length(unique(res.drugs[ maxPhaseForIndication == 4]$name))

## look into some
View(res.drugs[ type.snp == "top.shared", .(name, gene_name, gene_name.complex, type, embedding, trait)])

## drop one mapping manually: drops unspecific cancer GWAS
res.drugs  <- res.drugs[ efoName != "cancer"]

## write to file
write.table(res.drugs, "Intermediate.Results.Drug.Coloc.Embeddings.20260609.txt", sep = "\t", row.names = F)

#############################################
####   align with embedding 'meanings'   ####
#############################################

## pull in association results, to assign each embedding a 'disease-orientation'
res.embedding.obs <- fread("../../02_association_analysis/input/Results.Embedding.associations.UKB.minimal.extensive.20260319.txt.gz")
## write minimal results
write.table(res.embedding.obs[ pval.minimal < .05/nrow(res.embedding.obs), .(embedding, label_new, beta.minimal, se.minimal, pval.minimal, category)],
            "Results.Embedding.associations.UKB.minimal.sig.results.20260609.txt", sep = "\t", row.names = F)

#----------------------------------------------#
##-- orient embeddings via UKB associations --##
#----------------------------------------------#

## create mapping: embedding effect <-> drug indication
write.table(data.table(efoName = unique(res.drugs$efoName)), "Drug.Indication.UKB.exposure.association.mapping.20260608.txt", sep = "\t", row.names = F)
## import to map association direction
tmp.mapping       <- fread("Drug.Indication.UKB.exposure.association.mapping.20260608.txt")

## add to the data set
res.drugs         <- merge(res.drugs, tmp.mapping, by = "efoName")

## add sign to druggable data as orientation
res.drugs         <- merge(res.drugs, res.embedding.obs[, .(embedding, short_name_new, label_new, beta.minimal, se.minimal, pval.minimal)],
                           by = c("embedding", "short_name_new", "label_new"))

## define consistency with druggable effect
res.drugs[, druggable_direction := sign(BETA.embedding * BETA.drug * beta.minimal) > 0]

## write updated file to text
write.table(res.drugs[ type.snp == "top.shared" ], "Results.Drug.Coloc.Embeddings.20260609.txt", sep = "\t", row.names = F)

## ordinal "level of support": complex (weakest) < ligand.receptor < same (canonical)
str.lev <- c("complex", "ligand.receptor", "same")
res.drugs[, str.rank := match(type, str.lev)]               ## 1..3 for ordering / aggregation
## convenience boolean of the curated direction flag (TRUE = consistent/beneficial)
res.drugs[, consistent := (druggable_direction == TRUE)]

## write to file
write.table(res.drugs, "Intermediate.Results.Drug.Coloc.Embeddings.20260609.txt", sep = "\t", row.names = F)

#############################################
####          reporting manuscript       ####
#############################################

#-----------------------------#
##--    overall summary    --##
#-----------------------------#

## number of embeddings implied
uniqueN(res.drugs$embedding)
## N = 98

## number of unique loci
nrow(unique(res.drugs[, .(CHR_ENSEMBL, region_start.hg38)]))
## N = 71

## how many unique genes
uniqueN(res.drugs$gene_name.complex)
## N = 55

## how many canonical drug targets
uniqueN(res.drugs$gene_name)
## N = 81

## drugs
uniqueN(res.drugs$name)
## N = 217

## how many indications
uniqueN(res.drugs$disease_new)
## N = 36

## how many embedding - locus associations
nrow(unique(res.drugs[, .(embedding, CHR_ENSEMBL, region_start.hg38)]))
## n = 411

## how many emeddings with three or more distinct loci
sort(table(unique(res.drugs[, .(embedding, gene_name.complex)])$embedding))

## drugs involved
View(unique(res.drugs[, .(gene_name.complex, name, embedding)]))
sort(table(unique(res.drugs[, .(name, embedding)])$name))

## how many embedding - target - indication associations
nrow(unique(res.drugs[, .(embedding, gene_name.complex, efoName)]))

## directional consistency:
## count effects by embedding x instrumented gene x indication
fnd <- res.drugs[, .(str.rank = max(str.rank),                       # strongest support at the locus
                     npos     = sum(druggable_direction == TRUE),
                     nneg     = sum(druggable_direction == FALSE)),
                 by = .(embedding, gene_name.complex, efoName)]       # one row per causal locus x indication
fnd[, dir := fcase(nneg == 0, "consistent",                    # all GWAS agree: beneficial
                   npos == 0, "inconsistent",                  # all GWAS agree: adverse
                   default  = "mixed")]                         # GWAS disagree for this indication

#############################################
####            prepare figure           ####
#############################################

#-- therapeutic area and gene count (panel a grouping) --#

## add therapeutic areas
res.drugs[, area := fcase(
  grepl("hyperchol|hyperlip|atheroscler|coronary|myocard|angina|hypertens|heart failure|cardiovascular", tolower(efoName)), "Cardiovascular",
  grepl("metabolic|obesity|diabetes|thyroid",                    tolower(efoName)), "Metabolic/Endocrine",
  grepl("asthma|pulmonary",                                      tolower(efoName)), "Respiratory",
  grepl("rheumatoid|colitis|crohn|gout|arthr|eczema",            tolower(efoName)), "Immune",
  grepl("alzheimer|anxiety|schizophren|smoking",                 tolower(efoName)), "Neurologic",
  grepl("cancer|carcinoma|melanoma",                             tolower(efoName)), "Oncology",
  default = "Other")]
area.lev <- c("Cardiovascular", "Metabolic/Endocrine", "Respiratory", "Immune", "Neurologic", "Oncology", "Other")
area.col <- setNames(c("#9e0142", "#d53e4f", "#f46d43","#66c2a5","#3288bd","#5e4fa2","#bdbdbd"), area.lev)

## Each gene gets ONE dominant area (its most frequent indication area) so an
## embedding's bar is a clean stack over areas rather than double-counting genes. - this is the canonical target!
gene.area <- res.drugs[, .N, by = .(gene_name.complex, area)][ order(gene_name.complex, -N), .(area = area[1]), by = gene_name.complex]
## distinct druggable loci per embedding x area
emb.area  <- merge(unique(res.drugs[, .(gene_name.complex, embedding)]), gene.area, by = "gene_name.complex")[,
                                                                                                              .(n.loci = uniqueN(gene_name.complex)), by = .(embedding, area)]
emb.area[, area := factor(area, levels = area.lev)]
## per-embedding totals + its dominant area (for ordering and the side strip)
emb.sum   <- emb.area[, .(tot.loci = sum(n.loci), dom.loci = max(n.loci)), by = embedding]
emb.sum   <- merge(emb.sum, emb.area[ order(embedding, -n.loci)][, .(dom = area[1]), by = embedding], by = "embedding")
emb.sum[, `:=`(dom.share = dom.loci/tot.loci, dom = factor(dom, levels = area.lev),
               emb.num = as.integer(gsub("emb","",embedding)))]
## order embeddings by dominant area, then by no. of loci; y = top-to-bottom row
setorder(emb.sum, dom, -tot.loci, emb.num); emb.sum[, y := nrow(emb.sum):1]
## y-extent of each area block (for the coloured side strip)
dom.blk <- emb.sum[, .(y.lo = min(y), y.hi = max(y)), by = dom]

#-- aggregation panel b --#

# ## Direction is gene-independent (it depends on BETA.embedding/BETA.drug/beta.minimal,
# ## none of which vary by the partner gene), so collapse to one row per FINDING
# ## (emb.drug.region.id), taking the strongest support level present.
# find <- merge(res.drugs[, .(str.rank = max(str.rank)), by = emb.drug.region.id],
#               unique(res.drugs[, .(emb.drug.region.id, consistent)]), by = "emb.drug.region.id")
# sup  <- dcast(find[, .N, by = .(str.rank, consistent)], str.rank ~ consistent, value.var = "N", fill = 0)
# setnames(sup, c("FALSE","TRUE"), c("inconsistent","consistent"), skip_absent = TRUE)
# for(cc in c("consistent","inconsistent")) if(!cc %in% names(sup)) sup[, (cc) := 0]

seg.col <- c(consistent = "#E69F00", inconsistent = "#0072B2", mixed = "#bdbdbd")
sup <- dcast(fnd[, .N, by = .(str.rank, dir)], str.rank ~ dir, value.var = "N", fill = 0)
for(cc in names(seg.col)) if(!cc %in% names(sup)) sup[, (cc) := 0]
M <- t(as.matrix(sup[ match(1:3, str.rank), names(seg.col), with = FALSE])); M[is.na(M)] <- 0

#-- summary figure --##

pdf("../graphics/Summary.druggable.loci.Delphi2M.embeddings.20260609.pdf", width = 6.3, height = 6.3)

## define layout of the plot
layout(matrix(c(1,1,1,2,3,4,5,5), nrow = 2, byrow = TRUE), heights = c(.4,.6), widths = c(.3,.3,.2,.2))

#-- a: fingerprint --------------------------------------------------------------

## plotting parameters
par(mar = c(1.5, 1.5, .5, .5), mgp = c(.6,0,0), tck = -.01, cex.axis = .7, lwd = .5, cex.lab = .7, bty = "l", xaxs = "i")

## define boundaries
n.e <- nrow(emb.sum); y.max <- max(emb.sum$tot.loci)

## empty plot
plot(0, 0, type="n", xlim=rev(c(.5,n.e+.5)), ylim=c(0, y.max), axes=F, 
     ylab="druggable loci recovered", xlab="")

## stacked bar per embedding: one rectangle per area, widths = n.loci in that area
for(i in 1:n.e){ 
  ## what element to plot
  e  <- emb.sum$embedding[i]
  ## starting parameters
  xx <- emb.sum$y[i]
  yl <- 0
  ## loop through each area
  for(a in area.lev){ 
    ## get the width
    w <- emb.area[embedding==e & area==a]$n.loci
    ## add if any
    if(length(w)&&w>0){ 
      ## draw rectangle
      rect( xx-.5, yl, xx + .5, yl+w, col=area.col[a], border="white", lwd=.3)
      ## update
      yl <- yl+w 
    } 
  } 
}

## plotting coordinates
pm <- par("usr")

## coloured strip on the far left marking each dominant-area block
rect(dom.blk$y.lo-.5, pm[3]-(pm[4]*pm[3])*.10,  dom.blk$y.hi+.5, pm[3]-(pm[4]*pm[3])*.05, 
     col=area.col[as.character(dom.blk$dom)], border=NA, xpd = NA)
axis(2, lwd=.5)

## add labels
mtext("98 embeddings (grouped by dominant axis)", side=1, line=0, cex=.5)
## add legend
legend("topright", bty="n", cex=.6, inset=c(.01,.02), fill=area.col, border=NA, legend=area.lev)

#-- b: support x direction ------------------------------------------------------

## define colours: green = consistent (beneficial), red = inconsistent (potential adverse)
seg.col <- c(consistent = "#E69F00", inconsistent = "#0072B2", mixed = "#bdbdbd")

## stacked horizontal bars: consistent (green) + inconsistent (red), per support level
M       <- t(as.matrix(sup[ match(1:3, str.rank), names(seg.col), with = FALSE])); M[is.na(M)] <- 0

## empty plot
plot(c(.5, ncol(M)+.5), c(0, max(colSums(M))), type = "n", xlab = "", xlim=rev(c(.5, ncol(M)+.5)),
     ylab = "druggable findings (embedding \u00d7 causal gene \u00d7 indication)", xaxt = "n", yaxt = "n")

## add rectangle
rect(1:ncol(M)-.4, 0, 1:ncol(M)+.4, colSums(M), col = "#bdbdbd", border = "white", lwd=.3)
rect(1:ncol(M)-.4, 0, 1:ncol(M)+.4, colSums(M[1:2,]), col = "#0072B2", border = "white", lwd=.3)
rect(1:ncol(M)-.4, 0, 1:ncol(M)+.4, M[1,], col = "#E69F00", border = "white", lwd=.3)

## add axis
axis(2, lwd=.5); axis(1, at =1:3, labels = str.lev, lwd=.5)

## associated legend
legend("topleft", bty="n", cex=.6, 
       fill=seg.col, 
       border=NA,
       legend=c("consistent","inconsistent", "mixed across GWAS"), 
       ncol=2)

#-- c/d: effect across embeddings at druggable loci -----------------------------

## adopt graphical parameters
par(mar = c(2.5, 3.5, 1.5, .5), bty = "o", mgp = c(1.5,.2,0), xaxs = "i", yaxs = "r")

## PCSK9 -- via rs142130958 and hypercholesterolemia_phecode_272_11_2024_verma_a_science_gcst90479924_multiple

## temporary data to plot
tmp <- unique(res.drugs[ SNP == "rs142130958" & type.snp == "top.shared" & disease == "hypercholesterolemia_phecode_272_11_2024_verma_a_science_gcst90479924_multiple", 
                         .(embedding, BETA.embedding, SE.embedding, BETA.drug, beta.minimal, druggable_direction)])
## create direction adjused beta -> embedding lowering effect
tmp[, BETA.display := -1*(BETA.embedding * sign(BETA.drug) * sign(beta.minimal))]
## order
tmp <- tmp[ order(-BETA.display)]
tmp[, ci.l := BETA.display - 1.96 * SE.embedding ]
tmp[, ci.u := BETA.display + 1.96 * SE.embedding ]

## empty plot
plot(c(min(c(tmp$ci.l, 0)), max(c(tmp$ci.u, 0))), c(.5, nrow(tmp)), type = "n", ylab = "", xaxt = "n", yaxt = "n",
     xlab = "Drug-oriented genetic effect on embedding\n(per risk allele) \u00b1 95% CI")
## zero crossoing and axis
abline(v=0, lwd=.5); axis(1, lwd=.5)

## add confidence interval
arrows(tmp$ci.l, 1:nrow(tmp), tmp$ci.u, 1:nrow(tmp), lwd=.5, length = 0, 
       col = ifelse(tmp$druggable_direction, "#E69F00", "#0072B2"))
## add point estimates
points(tmp$BETA.display, 1:nrow(tmp), pch = 22, col = "white",
       bg = ifelse(tmp$druggable_direction, "#E69F00", "#0072B2"), 
       cex = .8)

## store plotting coordinates
pm <- par("usr")

## add embeddings
text(pm[1], 1:nrow(tmp), cex=.7, labels = tmp$embedding,
     pos = 2, xpd = NA)

## add text
mtext("PCSK9-inhibitors and hypercholesterolemia", side=3, line=.4, adj=0, cex=.6, font=1)


## TSLP -- via rs1837253 and asthma_2024_verma_a_science_gcst90480249_multiple

## temporary data to plot
tmp <- unique(res.drugs[ SNP == "rs1837253" & type.snp == "top.shared" & disease == "asthma_2024_verma_a_science_gcst90480249_multiple", 
                         .(embedding, BETA.embedding, SE.embedding, BETA.drug, beta.minimal, druggable_direction)])
## create direction adjused beta -> embedding lowering effect
tmp[, BETA.display := -1*(BETA.embedding * sign(BETA.drug) * sign(beta.minimal))]
## order
tmp <- tmp[ order(-BETA.display)]
tmp[, ci.l := BETA.display - 1.96 * SE.embedding ]
tmp[, ci.u := BETA.display + 1.96 * SE.embedding ]

## empty plot
plot(c(min(c(tmp$ci.l, 0)), max(c(tmp$ci.u, 0))), c(.5, nrow(tmp)), type = "n", ylab = "", xaxt = "n", yaxt = "n",
     xlab = "Drug-oriented genetic effect on embedding\n(per risk allele) \u00b1 95% CI")
## zero crossoing and axis
abline(v=0, lwd=.5); axis(1, lwd=.5)

## add confidence interval
arrows(tmp$ci.l, 1:nrow(tmp), tmp$ci.u, 1:nrow(tmp), lwd=.5, length = 0, 
       col = ifelse(tmp$druggable_direction, "#E69F00", "#0072B2"))
## add point estimates
points(tmp$BETA.display, 1:nrow(tmp), pch = 22, col = "white",
       bg = ifelse(tmp$druggable_direction, "#E69F00", "#0072B2"), 
       cex = .8)

## store plotting coordinates
pm <- par("usr")

## add embeddings
text(pm[1], 1:nrow(tmp), cex=.7, labels = tmp$embedding,
     pos = 2, xpd = NA)

## add text
mtext("Tezepelumab (TSLP) and asthma", side=3, line=.4, adj=0, cex=.6, font=1)


#-- e: most druggable embedding -----------------------------

## adopt graphical parameters
par(mar = c(2.5, 9, 1.5, .5))

## Embedding 92 -- across 12 loci

## temporary data to plot
tmp <- unique(res.drugs[ embedding == "emb92", 
                         .(embedding, BETA.embedding, SE.embedding, BETA.drug, beta.minimal, druggable_direction, gene_name, gene_name.complex, PP.H4.abf, name)])
## retain only top association by locus
tmp <- tmp[ order(gene_name.complex, -druggable_direction, -PP.H4.abf)]
tmp[, ind := 1:.N, by = "gene_name.complex"]
tmp <- tmp[ ind == 1]

## create direction adjused beta -> embedding lowering effect
tmp[, BETA.display := -1*(BETA.embedding * sign(BETA.drug) * sign(beta.minimal))]
## order
tmp <- tmp[ order(-BETA.display)]
tmp[, ci.l := BETA.display - 1.96 * SE.embedding ]
tmp[, ci.u := BETA.display + 1.96 * SE.embedding ]

## empty plot
plot(c(min(c(tmp$ci.l, 0)), max(c(tmp$ci.u, 0))), c(.5, nrow(tmp)), type = "n", ylab = "", xaxt = "n", yaxt = "n",
     xlab = "Drug-oriented genetic effect on embedding\n(per risk allele) \u00b1 95% CI")
## zero crossoing and axis
abline(v=0, lwd=.5); axis(1, lwd=.5)

## add confidence interval
arrows(tmp$ci.l, 1:nrow(tmp), tmp$ci.u, 1:nrow(tmp), lwd=.5, length = 0, 
       col = ifelse(tmp$druggable_direction, "#E69F00", "#0072B2"))
## add point estimates
points(tmp$BETA.display, 1:nrow(tmp), pch = 22, col = "white",
       bg = ifelse(tmp$druggable_direction, "#E69F00", "#0072B2"), 
       cex = .8)

## store plotting coordinates
pm <- par("usr")

## add embeddings
text(pm[1], 1:nrow(tmp), cex=.6, font = 3,
     labels = paste(tmp$gene_name.complex, "-\n(", tmp$name, "-",tmp$gene_name,")"), 
     pos = 2, xpd = NA)

## add text
mtext("Druggable loci - embedding 92", side=3, line=.4, adj=0, cex=.6, font=1)

## close device
dev.off()

#############################################
####    benchmark against ICD-10 coloc   ####
#############################################

#-------------------------#
##--  import findings  --##
#-------------------------#

## import results accordingly
res.coloc.icd10 <- fread("<path_to_file>")

## subset to findings with strong evidence
res.coloc.icd10 <- res.coloc.icd10[ PP.H4.abf >= .8 & R2.leads >= .8]

#-------------------------------------#
##--  harmonise to common columns   --#
#-------------------------------------#

## helper function
usplit          <- function(x, sep) { v <- unique(trimws(unlist(strsplit(as.character(x), sep)))); v[ v != "" & !is.na(v)] }

## EMBEDDINGS -> standard columns; drop the unspecific cancer GWAS as in the paper
E <- data.table(
  inst.gene  = res.drugs$gene_name.complex,
  canon.gene = res.drugs$gene_name,                       ## may carry '|' for complexes
  drug       = res.drugs$name,
  indication = res.drugs$efoName,
  locus      = paste(res.drugs$CHR_ENSEMBL, res.drugs$region_start.hg38, sep = ":"),
  phase      = suppressWarnings(as.integer(res.drugs$maxPhaseForIndication)))

## ICD-10 -> standard columns (already exploded per target x drug)
I <- data.table(
  inst.gene  = res.coloc.icd10$inference.gene,
  canon.gene = res.coloc.icd10$gene.single,
  drug       = res.coloc.icd10$name,
  indication = res.coloc.icd10$efoName,
  locus      = paste(res.coloc.icd10$chr, res.coloc.icd10$region_start.hg38, sep = ":"),
  phase      = suppressWarnings(as.integer(res.coloc.icd10$maxPhaseForIndication)))

#-------------------------------------#
##--        set builders            --#
#-------------------------------------#

## returns the unique member set for one level, optionally approved-only (phase 4)
memb <- function(dt, level, phase4 = FALSE){
  if(phase4) dt <- dt[ phase == 4]
  x <- switch(level,
              locus       = dt$locus,
              inst.gene   = dt$inst.gene,
              canon.gene  = unlist(lapply(dt$canon.gene, usplit, sep = "[|,]")),
              drug        = unlist(lapply(dt$drug,       usplit, sep = "[|]")),
              indication  = dt$indication,
              target.ind  = paste(dt$inst.gene, dt$indication, sep = " @ "))
  sort(unique(x[ !is.na(x) & x != ""]))
}

venn <- function(level, phase4 = FALSE){
  e <- memb(E, level, phase4); i <- memb(I, level, phase4)
  data.table(level = level,
             n.embeddings   = length(e),
             n.icd10        = length(i),
             both           = length(intersect(e, i)),
             icd10.only     = length(setdiff(i, e)),
             embeddings.only= length(setdiff(e, i)))
}

#-------------------------------------#
##--     the contrast table         --#
#-------------------------------------#

levs <- c("locus", "inst.gene", "canon.gene", "drug", "indication", "target.ind")
lab  <- c(locus="druggable loci (region)", inst.gene="instrumented genes",
          canon.gene="canonical targets", drug="drugs", indication="indications (EFO)",
          target.ind="target x indication")

contrast.all    <- rbindlist(lapply(levs, venn, phase4 = FALSE))
contrast.phase4 <- rbindlist(lapply(levs, venn, phase4 = TRUE))
contrast.all[,    level := lab[level]]
contrast.phase4[, level := lab[level]]

cat("\n================ ICD-10 WGS vs EMBEDDINGS (all phases) ================\n")
print(contrast.all)
cat("\n================ approved drugs only (max phase = 4) ================\n")
print(contrast.phase4)

#-------------------------------------#
##--  who is unique to each approach --#
#-------------------------------------#

## the manuscript lists: what the direct ICD-10 GWAS add, and what embeddings add
emb.only.gene <- setdiff(memb(E,"inst.gene"), memb(I,"inst.gene"))
icd.only.gene <- setdiff(memb(I,"inst.gene"), memb(E,"inst.gene"))
emb.only.drug <- setdiff(memb(E,"drug"),      memb(I,"drug"))
icd.only.drug <- setdiff(memb(I,"drug"),      memb(E,"drug"))

fwrite(data.table(inst.gene = emb.only.gene), "output/genes.embeddings.only.20260717.txt", sep="\t", row.names=F, quote=F)
fwrite(data.table(inst.gene = icd.only.gene), "output/genes.icd10.only.20260717.txt",      sep="\t", row.names=F, quote=F)
fwrite(data.table(drug      = emb.only.drug), "output/drugs.embeddings.only.20260717.txt", sep="\t", row.names=F, quote=F)
fwrite(data.table(drug      = icd.only.drug), "output/drugs.icd10.only.20260717.txt",      sep="\t", row.names=F, quote=F)

#-------------------------------------#
##--   manuscript-ready sentences   --#
#-------------------------------------#

g <- contrast.all[ level == "instrumented genes"]
d <- contrast.all[ level == "drugs"]

cat("\n---------------- for the manuscript ----------------\n")
cat(sprintf("Direct ICD-10 WGS colocalisation implicated %d instrumented target genes (%d approved-drug targets),\n",
            g$n.icd10, contrast.phase4[level=="instrumented genes"]$n.icd10))
cat(sprintf("of which %d (%.0f%%) were also recovered through the embedding GWAS.\n",
            g$both, 100 * g$both / max(g$n.icd10, 1)))
cat(sprintf("Embeddings additionally implicated %d target genes and %d drugs NOT seen by any direct ICD-10 GWAS;\n",
            g$embeddings.only, d$embeddings.only))
cat(sprintf("conversely %d genes and %d drugs were recovered ONLY by the direct ICD-10 GWAS.\n",
            g$icd10.only, d$icd10.only))
cat("----------------------------------------------------\n")

# Direct ICD-10 WGS colocalisation implicated 133 instrumented target genes (81 approved-drug targets),
# of which 45 (34%) were also recovered through the embedding GWAS.
# Embeddings additionally implicated 10 target genes and 12 drugs NOT seen by any direct ICD-10 GWAS;
# conversely 88 genes and 258 drugs were recovered ONLY by the direct ICD-10 GWAS.

#-------------------------------------#
##--           add column           --#
#-------------------------------------#

## add column to res.drugs
res.drugs[, drug.covered.by.single.icd10.gwas := name %in% I$drug]
res.drugs[, locus.covered.by.single.icd10.gwas := gene_name.complex %in% I$inst.gene]
## write to file
write.table(res.drugs, "Intermediate.Results.Drug.Coloc.Embeddings.20260721.txt", sep = "\t", row.names = F)

