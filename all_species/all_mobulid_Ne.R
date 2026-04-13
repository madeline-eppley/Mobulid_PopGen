### all mobulid species Ne / effective population size analysis ####
# M. Eppley, v1. current version uploaded to github 4/13/2026

outdir <- "/Users/madelineeppley/Desktop/manta26pub"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
setwd(outdir)

pop_colors <- c("IND"="#6DBAA4","MEX"="#ED9A6C","NIC"="#8C9FCB","ECU"="#DA8EC0","BYC"="#B8B8B8")
discrete_palette <- c("#264653","#2a9d8f","#8ab17d","#f4a261","#e76f51")

# i'm going to run a test on the thurstoni data first
vcf_thur <- read.vcfR("/Users/madelineeppley/Desktop/manta/finalvcfs/thurstoni_final_v2.vcf")
popmap <- data.frame(
  sample = c("mthur_446","IN_10_MB_B","IN_12_MB_B",
             "MT_07_ELP","MT_04_ELP","MT_02_ELP","MT_05_ELP",
             "NIC_6003_E","NIC_6048","NIC_6003_F","NIC_6016_1",
             "MT_01_EC","MT_05_EC","MT_04_EC_conc","MT_07_EC_conc",
             "MT_10_EC","MT_09_EC_conc","MT_08_EC_conc","MT_03_EC",
             "BYC_RMO_40","BYC_RMO_47","BYC_RMO_48","BYC_RMU_16",
             "BYC_RMU_17","BYC_RMU_18","BYC_RMU_26","BYC_DESC_60","BYC_DESC_62"),
  pop = c(rep("IND",3),rep("MEX",4),rep("NIC",4),rep("ECU",8),rep("BYC",9)))

gl_thur <- vcfR2genlight(vcf_thur)
pop(gl_thur) <- popmap$pop[match(indNames(gl_thur),popmap$sample)]
n_snps_all <- nrow(vcf_thur@fix)

library(dartR)
gl.LDNe(
  gl_thur,
  outfile = "thurstoni_genepopLD.txt",
  outpath = tempdir(),
  neest.path = "/Users/madelineeppley/Desktop/manta26pub/NeEstimator",
  critical = 0,
  singleton.rm = TRUE,
  mating = "random",
  plot.out = TRUE,
  plot_theme = theme_dartR(),
  plot_colors_pop = discrete_palette,
  save2tmp = FALSE)

## this code worked, so i'm going to save it for future uses, but we don't have enough samples per-pop to run the Ne analyis
## going to make an all mobulids vcf, and the populations will be species, then re-run the Ne analysis on all

## all mobulids analysis
library(vcfR)
library(adegenet)

vcf <- read.vcfR("/Users/madelineeppley/Desktop/manta/finalvcfs/all_mobulids_final.vcf")
vcf_samples <- colnames(vcf@gt)[-1]
nrow(vcf@fix) # 22771 SNPs
length(vcf_samples) # 96 individuals

# sanity check the removed individuals (these were removed in vcftools filtering)
removed <- c("IN_4_MB", "REV_10_MB", "PER_DPPW28_MM", "mthur_446",
             "MT_04_ELP", "MT_02_ELP", "MT_05_ELP", "NIC_6003_F",
             "NIC_6016_1", "MT_09_EC_conc", "MT_08_EC_conc", "MT_03_EC",
             "BYC_RMO_40")
still_present <- intersect(removed, vcf_samples)
if (length(still_present) == 0) {
  cat("PASS: All 13 removed individuals are absent from VCF\n\n")} else {
  cat("WARNING: These removed individuals are still in VCF:", still_present, "\n\n")}

# pop map with all species
popmap <- data.frame(
  sample = c(
    "BYC_RMB_01","IN_1_MB","IN_2_MB","IN_5_MB","IN_6_MB","IN_7_MB",
    "PER_001_MB","PER_003_MB","PER_004_MB","PER_005_MB","PER_006_MB",
    "PER_007_MB","PER_DZW81_4_MB",
    "REV_13_MB","REV_14_MB_B","REV_15_MB","REV_17_MB","REV_18_MB",
    "REV_19_MB","REV_20_MB_B",
    "mu_11","mu_12","mu_21","mu_3","mu_4","mu_7",
    "NIC_2728_01","NIC_2740_01",
    "PER_DCLCW130_MM","PER_DCLCW87_MM","PER_DPPW019_MM","PER_DPPW98_MM",
    "IN_10_MB_B","IN_12_MB_B",
    "MT_07_ELP",
    "NIC_6003_E","NIC_6048",
    "MT_01_EC","MT_05_EC","MT_04_EC_conc","MT_07_EC_conc","MT_10_EC",
    "BYC_RMO_47","BYC_RMO_48","BYC_RMU_16","BYC_RMU_17","BYC_RMU_18",
    "BYC_RMU_26","BYC_DESC_60","BYC_DESC_62",
    "BYC_RMM_10","BYC_RMM_11","BYC_RMM_12","BYC_RMM_13","BYC_RMM_14",
    "BYC_RMM_20","BYC_RMM_21","BYC_RMM_22","BYC_RMM_24","BYC_RMM_25",
    "BYC_RMM_31","BYC_RMM_32","BYC_RMB_34","BYC_RMM_35","BYC_RMB_36",
    "BYC_RMM_38","BYC_RMU_41","BYC_RMU_42","BYC_RMM_44","BYC_RMM_50",
    "BYC_RMM_51","BYC_DESC_54","BYC_DESC_55","BYC_RMM_58","BYC_DESC_63",
    "BYC_DESC_64","BYC_RMM_02","BYC_RMM_03","BYC_RMM_05",
    "BYCI_RMM_68","BYCI_RMM_66","BYCI_DESC_67","BYCI_RMM_72",
    "BYC_RMB_57","BYC_RMO_45","BYC_RMT_27","BYC_RMT_29","BYC_RMT_46",
    "BYC_RMM_30","BYC_RMT_07","BYC_RMT_28","BYC_RMT_49",
    "BYC_RMT_04","BYC_RMT_06","BYCI_RMT_69","BYCI_RMT_71"),
  pop = c(rep("birostris",20), rep("munkiana",12), rep("thurstoni",18),
          rep("mobular",33), rep("tarapacana",13)))

popmap <- popmap[popmap$sample %in% vcf_samples, ]
print(table(popmap$pop))

# check for any inds not in pop map above
unmatched <- setdiff(vcf_samples, popmap$sample)
if (length(unmatched) > 0) {
  cat("WARNING - samples in VCF but not in popmap:", unmatched, "\n\n")} else {
  cat("PASS: All VCF samples matched in popmap\n\n")}

# per-individual missingness
gt_matrix <- extract.gt(vcf, element = "GT")
miss_per_ind <- apply(gt_matrix, 2, function(x) sum(is.na(x)) / length(x))
miss_df <- data.frame(
  sample = names(miss_per_ind),
  missing = round(miss_per_ind, 4),
  species = popmap$pop[match(names(miss_per_ind), popmap$sample)])
miss_df <- miss_df[order(miss_df$species, -miss_df$missing), ]

# check
print(miss_df, row.names = FALSE)
max(miss_df$missing) # 0.3621, filtering cutoff at 40% worked
miss_df$sample[which.max(miss_df$missing)] # most missing ind is REV_18_MB
mean(miss_df$missing) # 0.0744 great

# per species missingness summary
print(tapply(miss_df$missing, miss_df$species, function(x) 
  round(c(mean = mean(x), median = median(x), max = max(x), n = length(x)), 4)))
# mostly thurstoni and munkiana are the higher missingness species
# thurstoni median 0.1184, mean 0.1382
# munkiana median 0.1054, mean 0.1118

# per locus missingness summary
miss_per_locus <- apply(gt_matrix, 1, function(x) sum(is.na(x)) / length(x))
print(summary(miss_per_locus))
sum(miss_per_locus > 0.2) # 1 locus
sum(miss_per_locus > 0.5) # 0, filtering worked


# loci shared across species
# at each locus, check which species have at least one genotyped individual
for (sp in c("birostris", "munkiana", "thurstoni", "mobular", "tarapacana")) {
  sp_samples <- popmap$sample[popmap$pop == sp]
  sp_gt <- gt_matrix[, sp_samples, drop = FALSE]
  n_loci <- sum(apply(sp_gt, 1, function(x) any(!is.na(x))))
  cat("  ", sp, ":", n_loci, "of", nrow(gt_matrix), 
      "(", round(100 * n_loci / nrow(gt_matrix), 1), "%)\n")}

# data looks really nice and complete
# birostris : 22771 of 22771 (100 %)
# munkiana : 22067 of 22771 (96.9 %)
# thurstoni : 22771 of 22771 (100 %)
# mobular : 22771 of 22771 (100 %)
# tarapacana : 22318 of 22771 (98 %)

# count loci present by species
loci_by_species <- sapply(c("birostris","munkiana","thurstoni","mobular","tarapacana"),
                          function(sp) {
                            sp_samples <- popmap$sample[popmap$pop == sp]
                            apply(gt_matrix[, sp_samples, drop = FALSE], 1, function(x) any(!is.na(x)))})
n_species_per_locus <- rowSums(loci_by_species)
print(table(n_species_per_locus))
# 1157 loci have 4 species
# 21614 loci have all 5 species

## all Mobulids effective population size
library(dartR)
library(ggplot2)

outdir <- "/Users/madelineeppley/Desktop/manta26pub"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

vcf <- read.vcfR("/Users/madelineeppley/Desktop/manta/finalvcfs/all_mobulids_final.vcf")
nrow(vcf@fix)
ncol(vcf@gt) - 1

# popmap
popmap <- data.frame(
  sample = c(
    "BYC_RMB_01","IN_1_MB","IN_2_MB","IN_5_MB","IN_6_MB","IN_7_MB",
    "PER_001_MB","PER_003_MB","PER_004_MB","PER_005_MB","PER_006_MB",
    "PER_007_MB","PER_DZW81_4_MB",
    "REV_13_MB","REV_14_MB_B","REV_15_MB","REV_17_MB","REV_18_MB",
    "REV_19_MB","REV_20_MB_B",
    "mu_11","mu_12","mu_21","mu_3","mu_4","mu_7",
    "NIC_2728_01","NIC_2740_01",
    "PER_DCLCW130_MM","PER_DCLCW87_MM","PER_DPPW019_MM","PER_DPPW98_MM",
    "IN_10_MB_B","IN_12_MB_B",
    "MT_07_ELP",
    "NIC_6003_E","NIC_6048",
    "MT_01_EC","MT_05_EC","MT_04_EC_conc","MT_07_EC_conc","MT_10_EC",
    "BYC_RMO_47","BYC_RMO_48","BYC_RMU_16","BYC_RMU_17","BYC_RMU_18",
    "BYC_RMU_26","BYC_DESC_60","BYC_DESC_62",
    "BYC_RMM_10","BYC_RMM_11","BYC_RMM_12","BYC_RMM_13","BYC_RMM_14",
    "BYC_RMM_20","BYC_RMM_21","BYC_RMM_22","BYC_RMM_24","BYC_RMM_25",
    "BYC_RMM_31","BYC_RMM_32","BYC_RMB_34","BYC_RMM_35","BYC_RMB_36",
    "BYC_RMM_38","BYC_RMU_41","BYC_RMU_42","BYC_RMM_44","BYC_RMM_50",
    "BYC_RMM_51","BYC_DESC_54","BYC_DESC_55","BYC_RMM_58","BYC_DESC_63",
    "BYC_DESC_64","BYC_RMM_02","BYC_RMM_03","BYC_RMM_05",
    "BYCI_RMM_68","BYCI_RMM_66","BYCI_DESC_67","BYCI_RMM_72",
    "BYC_RMB_57","BYC_RMO_45","BYC_RMT_27","BYC_RMT_29","BYC_RMT_46",
    "BYC_RMM_30","BYC_RMT_07","BYC_RMT_28","BYC_RMT_49",
    "BYC_RMT_04","BYC_RMT_06","BYCI_RMT_69","BYCI_RMT_71"),
  pop = c(rep("birostris",20), rep("munkiana",12), rep("thurstoni",18),
          rep("mobular",33), rep("tarapacana",13)))



# convert to genlight and assign pops
gl_all <- vcfR2genlight(vcf)
pop(gl_all) <- popmap$pop[match(indNames(gl_all), popmap$sample)]

print(table(pop(gl_all)))

# run Ne estimation single call, estimates Ne per population
ne_results <- gl.LDNe(
  gl_all,
  outfile = "all_mobulids_genepopLD.txt",
  outpath = outdir,
  neest.path = "/Users/madelineeppley/Desktop/manta26pub/NeEstimator",
  critical = 0,
  singleton.rm = TRUE,
  mating = "random",
  plot.out = TRUE,
  plot_theme = theme_dartR(),
  plot_colors_pop = discrete_palette,
  save2tmp = FALSE)

ne_results
saveRDS(ne_results, file = file.path(outdir, "all_mobulids_Ne_results.rds"))
