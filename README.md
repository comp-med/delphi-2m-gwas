# Genetic decoding reveals druggable biology implicitly learned by a medical-history foundation model

This repository contains the analysis code accompanying the manuscript.

> Zeng W, Kohleick L, Beuchel C, Zoodsma M, Koprulu M, Carrasco Zanini J, Wild B, Langenberg C, Pietzner M. *Genetic decoding reveals druggable biology implicitly learned by a medical-history foundation model.* (2026).

---

#### Note

*The provided scripts are not designed to work out of the box. Analyses were run across a SLURM high-performance computing cluster (BIH, Berlin) using a mixture of native genomics binaries, a singularity container for R, and conda environments for Python and Snakemake, and they depend on UK Biobank individual-level data that cannot be redistributed. The scripts are intended to document the analytical steps used to generate the results reported in the manuscript. All file paths have been replaced with `<path_to_file>`, or with an informative placeholder where the program being called matters (`<path_to_regenie>`, `<path_to_plink2>`, `<path_to_container>`, `<path_to_ukb_bgen>`). Individual-level UK Biobank data are available to bona fide researchers on application via [ukbiobank.ac.uk](https://www.ukbiobank.ac.uk/); this research was conducted under application no. 44448.*

*Candidate effector genes were assigned using an ensemble variant-to-gene classifier that is developed and reported in separate work; that code is therefore not included here. The scripts in `04_credible_set_annotation/` import the resulting per-locus assignments and document every downstream step that used them.*

*Delphi-2M itself was not developed by us. We re-implemented it following the repository released by the original authors ([gerstung-lab/Delphi](https://github.com/gerstung-lab/Delphi), MIT licence). Only the files we modified relative to that repository are included here, under `01_reimplementation_Delphi-2M/delphi_modified/`, at their upstream-relative paths.*

---

### Project structure

| Script | Description |
| --- | --- |
| **01 — Re-implementation of Delphi-2M** | |
| [`01_reimplementation_Delphi-2M/01_train_Delphi.sh`](01_reimplementation_Delphi-2M/01_train_Delphi.sh) | Launch training of Delphi-2M on the UK Biobank ICD-10 first-occurrence data using the upstream `train.py` on a single CUDA device |
| [`01_reimplementation_Delphi-2M/02_embedding_generation.py`](01_reimplementation_Delphi-2M/02_embedding_generation.py) | Use the trained model as an encoder: mean-pool the final-layer hidden states across all non-padding tokens per participant to yield the 120-dimensional participant-level embeddings, and export the linear predictors used downstream |
| [`01_reimplementation_Delphi-2M/delphi_modified/data/example_ukb_to_bin.py`](01_reimplementation_Delphi-2M/delphi_modified/data/example_ukb_to_bin.py) | Convert UK Biobank first-occurrence fields, sex and lifestyle tokens into the binary trajectory format expected by Delphi-2M; 80/20 train/validation split (modified from upstream) |
| [`01_reimplementation_Delphi-2M/delphi_modified/config/train_delphi.py`](01_reimplementation_Delphi-2M/delphi_modified/config/train_delphi.py) | Training configuration (12 layers, block size 48, batch size 128, seed 42) used for the re-implementation (modified from upstream) |
| [`01_reimplementation_Delphi-2M/delphi_modified/evaluate_additional_metrics.py`](01_reimplementation_Delphi-2M/delphi_modified/evaluate_additional_metrics.py) | Extend the upstream AUC-only evaluation with precision, recall, F1, AUPRC and balanced accuracy, retaining disease chunking, sex stratification and age-bracket evaluation (modified from upstream) |
| **02 — Genome-wide association testing** | |
| [`02_GWAS/01_REGENIE_step1_emb.sh`](02_GWAS/01_REGENIE_step1_emb.sh) | REGENIE step 1 (whole-genome ridge regression) across all 120 embeddings in the combined-sex sample |
| [`02_GWAS/02_REGENIE_step2_emb.sh`](02_GWAS/02_REGENIE_step2_emb.sh) | REGENIE step 2 single-variant association testing per chromosome against imputed genotypes, combined sexes; chromosome 23 handled as X |
| [`02_GWAS/03_REGENIE_step1_emb_female.sh`](02_GWAS/03_REGENIE_step1_emb_female.sh) | REGENIE step 1, female-only sample |
| [`02_GWAS/04_REGENIE_step2_emb_female.sh`](02_GWAS/04_REGENIE_step2_emb_female.sh) | REGENIE step 2, female-only sample |
| [`02_GWAS/05_REGENIE_step1_emb_male.sh`](02_GWAS/05_REGENIE_step1_emb_male.sh) | REGENIE step 1, male-only sample |
| [`02_GWAS/06_REGENIE_step2_emb_male.sh`](02_GWAS/06_REGENIE_step2_emb_male.sh) | REGENIE step 2, male-only sample |
| **03 — Signal selection and statistical fine-mapping** | |
| [`03_signal_selection/01_proc_gwas.sh`](03_signal_selection/01_proc_gwas.sh) | Post-process REGENIE output: group results by embedding, intersect variants with the published UK Biobank whole-genome-sequencing GWAS by rsID, filter and sort, define regional lead signals with a dedicated extended-MHC rule on chromosome 6, and write the region files used for fine-mapping |
| [`03_signal_selection/02_run_smk_process_regions.sh`](03_signal_selection/02_run_smk_process_regions.sh) | Launch the Snakemake workflow that prepares regional summary statistics and genotype dosages per embedding |
| [`03_signal_selection/03_run_smk_finemapping.sh`](03_signal_selection/03_run_smk_finemapping.sh) | Launch the Snakemake fine-mapping workflow across all embeddings and regions |
| [`03_signal_selection/04_finemapping_workflow.smk`](03_signal_selection/04_finemapping_workflow.smk) | Snakemake workflow definition dispatching SuSiE per embedding–region pair |
| [`03_signal_selection/05_run_susie.R`](03_signal_selection/05_run_susie.R) | Fine-map one region with `susieR`, sweeping the number of causal signals and retaining the model with the most credible sets; exclude credible sets whose lead variants are in LD (r²>0.25); refit a joint model on the remaining lead variants and keep only signals that stay genome-wide significant with concordant direction and a joint estimate within 25% of the marginal estimate; export marginal and joint statistics, credible-set membership, PIP and LD with the lead variant |
| **04 — Credible-set annotation, effector genes and enrichment** | |
| [`04_credible_set_annotation/01_annotate_credible_set_variants.R`](04_credible_set_annotation/01_annotate_credible_set_variants.R) | Head script for locus definition and annotation: LD-clump credible-set lead variants into independent loci via a graph-component approach, fold in regional sentinels and the MHC, harmonise statistics across embeddings, assign effector genes from variant-to-gene evidence with a nearest-gene fallback, lift over to build 38 and intersect with the NHGRI-EBI GWAS Catalog via LD proxies, and run EFO-term and ICD-10 enrichment testing with redundancy pruning |
| [`04_credible_set_annotation/02_extract_proxies.sh`](04_credible_set_annotation/02_extract_proxies.sh) | SLURM array job computing LD proxies (r²>0.1, 3 Mb window) for all credible-set variants in an unrelated European-ancestry subsample using PLINK2 |
| [`04_credible_set_annotation/03_clump_gwas_signals.sh`](04_credible_set_annotation/03_clump_gwas_signals.sh) | SLURM array job clumping GWAS Catalog variants into loci with PLINK2 for the enrichment analysis |
| [`04_credible_set_annotation/04_wgs_icd10_lookup.sh`](04_credible_set_annotation/04_wgs_icd10_lookup.sh) | SLURM array job looking up credible-set variants across 763 published UK Biobank whole-genome-sequencing ICD-10 GWAS |
| [`04_credible_set_annotation/obtain_snps.sh`](04_credible_set_annotation/obtain_snps.sh) | Extract genotype dosages for a supplied variant list from the UK Biobank BGEN files |
| **05 — Druggable target assignment and colocalisation** | |
| [`05_druggable_assignment/01_map_druggable_targets.R`](05_druggable_assignment/01_map_druggable_targets.R) | Head script for the drug analysis: parse Open Targets drug, target and indication annotations, restrict to approved medications, map targets to genes and to embedding loci, define the drug-proxying regions to test, collate colocalisation results, orient embeddings against their UK Biobank trait associations to classify effects as consistent or adverse, and benchmark the recovered drug set against colocalisation with single ICD-10 GWAS |
| [`05_druggable_assignment/02_submit_drug_coloc.sh`](05_druggable_assignment/02_submit_drug_coloc.sh) | SLURM array job running colocalisation with one task per embedding–drug-locus pair |
| [`05_druggable_assignment/03_run_drug_coloc.R`](05_druggable_assignment/03_run_drug_coloc.R) | Run `coloc` for a single embedding–drug-locus pair, extract the regional summary statistics, compute LD and draw stacked locus-compare plots for strong colocalisations |
| [`05_druggable_assignment/02_submit_drug_coloc.sh`](05_druggable_assignment/02_submit_drug_coloc.sh) | SLURM array job running colocalisation with one task per genomic region |
| [`05_druggable_assignment/03_run_drug_coloc.R`](05_druggable_assignment/03_run_drug_coloc.R) | Run `coloc` for all drug loci within one genomic region against all embeddings, gated by a regional signal check (p<10-5) before the expensive merge and colocalisation steps; compute LD and draw stacked locus-compare plots for strong colocalisations |
| [`05_druggable_assignment/obtain_ld_matrix.sh`](05_druggable_assignment/obtain_ld_matrix.sh) | Compute the LD matrix among variants in a region from the UK Biobank BGEN files |
| **06 — Predictive augmentation of Delphi-2M** | |
| [`06_predictive_augmentation/01_data_prep_collate_results.R`](06_predictive_augmentation/01_data_prep_collate_results.R) | Head script for the prediction analyses: prepare UK Biobank phenotypes, biomarkers and first-occurrence outcomes for survival analysis, drop prevalent and early events, build the partial-correlation disease network summaries, collate Cox and ΔAUC results, relate discrimination gain to how well each marker is reconstructed from the embeddings, and generate the reporting tables and figures |
| [`06_predictive_augmentation/02_submit_cox.sh`](06_predictive_augmentation/02_submit_cox.sh) | SLURM array job for the Cox models across outcomes |
| [`06_predictive_augmentation/03_run_cox_model.R`](06_predictive_augmentation/03_run_cox_model.R) | Fit Cox proportional-hazards models for one outcome across all candidate markers, with and without the Delphi-2M linear predictor |
| [`06_predictive_augmentation/04_collate_delta_auc_input.R`](06_predictive_augmentation/04_collate_delta_auc_input.R) | Build the outcome × marker × population parameter file that drives the ΔAUC job array |
| [`06_predictive_augmentation/04_submit_parcor.sh`](06_predictive_augmentation/04_submit_parcor.sh) | SLURM job for the partial-correlation disease network |
| [`06_predictive_augmentation/05_compute_partial_correlation_disease.R`](06_predictive_augmentation/05_compute_partial_correlation_disease.R) | Compute the partial-correlation network among lifetime ICD-10 codes, conditioning each disease pair on age, sex and all other diseases via the precision matrix, after pruning near-duplicate diseases (|r|>0.7) to keep the matrix invertible |
| [`06_predictive_augmentation/06_submit_delta_auc.sh`](06_predictive_augmentation/06_submit_delta_auc.sh) | SLURM array job for the ΔAUC analysis |
| [`06_predictive_augmentation/07_run_delta_auc_cv.R`](06_predictive_augmentation/07_run_delta_auc_cv.R) | Cross-fitted ΔAUC for one outcome, restricted to continuous markers; this is the estimate reported in the manuscript |
| [`06_predictive_augmentation/08_submit_r2_embeddings.sh`](06_predictive_augmentation/08_submit_r2_embeddings.sh) | SLURM job for the marker-reconstruction analysis |
| [`06_predictive_augmentation/09_compute_r2_embeddings.R`](06_predictive_augmentation/09_compute_r2_embeddings.R) | Quantify how much of each risk factor and biomarker can be reconstructed from the 120 embeddings, within sex and pooled |
| [`06_predictive_augmentation/10_submit_cox_proteins.sh`](06_predictive_augmentation/10_submit_cox_proteins.sh) | SLURM array job for the Cox models in the Olink proteomic subcohort |
| [`06_predictive_augmentation/11_run_cox_proteins.R`](06_predictive_augmentation/11_run_cox_proteins.R) | Cox models for one outcome using plasma proteins as exposures |
| [`06_predictive_augmentation/12_submit_delta_auc_proteins.sh`](06_predictive_augmentation/12_submit_delta_auc_proteins.sh) | SLURM array job for the protein ΔAUC analysis |
| [`06_predictive_augmentation/13_run_delta_auc_cv_proteins.R`](06_predictive_augmentation/13_run_delta_auc_cv_proteins.R) | Cross-fitted ΔAUC for plasma proteins in the Olink subcohort, with a global age- and sex-adjusted fallback for sparse outcomes |
| [`06_predictive_augmentation/14_submit_r2_embeddings_proteins.sh`](06_predictive_augmentation/14_submit_r2_embeddings_proteins.sh) | SLURM job for protein reconstruction from the embeddings |
| [`06_predictive_augmentation/15_compute_r2_embeddings_proteins.R`](06_predictive_augmentation/15_compute_r2_embeddings_proteins.R) | Quantify how much of each plasma protein can be reconstructed from the embeddings, mirroring script 09 in the Olink subcohort |
| **Functions** | |
| [`functions/run_cox_models.R`](functions/run_cox_models.R) | `run.cox()` — parallel Cox proportional-hazards models across a set of exposures, with automatic detection of categorical exposures, rank-inverse-normal transformation of continuous exposures, optional recruitment-centre handling (fixed effect, frailty or stratification) and optional adjustment for the Delphi-2M linear predictor |
| [`functions/delphi_delta_auc_cv.R`](functions/delphi_delta_auc_cv.R) | Cross-fitted ΔAUC with out-of-sample prediction within 5-year-age × sex strata, a global age- and sex-adjusted fallback when strata collapse, DeLong comparison of paired AUCs and reporting of the component AUCs |
| [`functions/improvement_report.R`](functions/improvement_report.R) | Assemble the reporting analysis linking marker reconstruction, absorbed predictive value and residual gain, run jointly over biomarkers and plasma proteins |
| [`functions/reconstruction_vs_headroom_perpair.R`](functions/reconstruction_vs_headroom_perpair.R) | Per-pair comparison of how well a marker is reconstructed from the embeddings against how much discrimination it adds, with a case-count floor, plus the corresponding scatter plot |
| [`functions/regression_adaptive.R`](functions/regression_adaptive.R) | `regression_analysis()` — unified regression across multiple features and outcomes with automatic outcome-type detection (linear, logistic, negative binomial or multinomial), optional rank-inverse-normal transformation, and linear plus spline models returned in long format |
| [`functions/plot_regression_scatter.R`](functions/plot_regression_scatter.R) | Scatter plot of a feature against a continuous outcome with the fitted linear or spline line overlaid, re-fitting through the same pipeline so the plot matches the reported estimates |
| [`functions/map_closest_gene.R`](functions/map_closest_gene.R) | `map_nearest_gene()` — assign the nearest gene to each variant, optionally keeping ties |
| [`functions/annotate_gwas_catalog.R`](functions/annotate_gwas_catalog.R) | Annotate a variant set with GWAS Catalog associations, in parallel, allowing for LD proxies above a supplied r² threshold |
| [`functions/gene_enrichment_hpa.R`](functions/gene_enrichment_hpa.R) | Multi-method tissue and cell-type expression enrichment for a gene list against any long-format expression table (for example Human Protein Atlas consensus RNA), combining Fisher's exact, rank-based and distributional tests |
| [`functions/pathway_pruning.R`](functions/pathway_pruning.R) | Prune redundant pathways within each embedding by greedy set cover, skipping candidates whose gene sets exceed a Jaccard overlap threshold with an already selected pathway |
| [`functions/prune_consistent_by_factor.R`](functions/prune_consistent_by_factor.R) | Prune enrichment results within strata using a globally consistent trait priority, so the same trait is retained across embeddings |
| [`functions/prune_with_tracking.R`](functions/prune_with_tracking.R) | Prune overlapping enrichment results for one trait while recording which terms were collapsed into each retained row |
| [`functions/fn_drug_coloc.R`](functions/fn_drug_coloc.R) | Shared functions for the embedding × drug colocalisation, including the cheap regional signal check used to skip regions before the expensive colocalisation step |
| [`functions/extract_genomic_region.R`](functions/extract_genomic_region.R) | Extract a genomic region from formatted parquet summary statistics using polars (provided by Carl Beuchel) |
| [`functions/extract_genomic_data.R`](functions/extract_genomic_data.R) | As above, extended to extract by variant list as well as by coordinate range |
| [`functions/import_ld_matrix_regenie.R`](functions/import_ld_matrix_regenie.R) | `get.corr.sq.matrix()` — read the compressed binary LD matrix written by REGENIE together with its variant list |
| [`functions/get_dosages.R`](functions/get_dosages.R) | `get_dosages()` — call the regional dosage extraction, read the resulting dosage file, transpose it to one row per participant and rebuild the allele information needed to realign effect estimates |
| [`functions/function_wakefield.R`](functions/function_wakefield.R) | `doWakefield()` — fallback used when SuSiE fails to converge for a region: derive approximate Bayes factors and posterior inclusion probabilities from the marginal statistics, form the 95% credible set, add LD to the lead variant, and apply the same joint-model filter (genome-wide significance, concordant direction, joint estimate within 25% of the marginal) as the SuSiE path |
| [`functions/plot_locus_compare.R`](functions/plot_locus_compare.R) | Stacked locus-zoom plot comparing embedding and drug-locus association statistics, coloured by LD with the lead variant |

---

### Software

| Software | Version | Purpose |
| --- | --- | --- |
| Python | 3.11 | Delphi-2M re-implementation and embedding generation |
| PyTorch | — | Delphi-2M model training and inference |
| R | 4.3.2 | All statistical analyses downstream of the GWAS |
| REGENIE | 4.0 | Genome-wide association testing and LD computation |
| PLINK2 | 2.00a | LD proxies, clumping and genotype extraction |
| BEDtools | 2.30.0 | Definition and merging of genomic regions |
| bgenix / qctool | — | Extraction of genotype dosages from BGEN files |
| Snakemake | — | Orchestration of the regional fine-mapping workflow |
| susieR | 0.14.2 | Statistical fine-mapping (`susie_rss`) |
| corrcoverage | — | Approximate Bayes factors and credible sets for the Wakefield fallback |
| coloc | 5.2.3 | Bayesian colocalisation |
| survival | — | Cox proportional-hazards models |
| pROC | — | AUC estimation and DeLong comparison of paired AUCs |
| data.table, arrow, polars | — | Data handling and parquet access |
| igraph | — | Graph-component approach to locus definition |
| gprofiler2 | — | Pathway, tissue and cell-type enrichment |

---

### Licence

Released under the MIT Licence; see [LICENSE](LICENSE). Code under `01_reimplementation_Delphi-2M/delphi_modified/` is modified from [gerstung-lab/Delphi](https://github.com/gerstung-lab/Delphi) and remains subject to its original MIT Licence.

---

### Citation

Zeng W, Kohleick L, Beuchel C, Zoodsma M, Koprulu M, Carrasco Zanini J, Wild B, Langenberg C, Pietzner M. *Genetic decoding reveals druggable biology implicitly learned by a medical-history foundation model.* (2026).
