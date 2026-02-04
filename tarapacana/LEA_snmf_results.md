```R
###################
# LEA snmf structure
###################

# convert vcf to geno format for LEA
vcf2geno("/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_4pop.vcf",
         "/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_4pop.geno")

# run snmf
project_4pop <- snmf("/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_4pop.geno",
                     K = 1:5, entropy = TRUE, repetitions = 10, project = "new")

print(summary(project_4pop))
plot(project_4pop, col = "blue", pch = 19, cex = 1.2, main = "All SNPs: M. tarapacana")
```

<img width="1470" height="1296" alt="image" src="https://github.com/user-attachments/assets/fc586a51-2b24-4ca7-a7ef-fc12d92c9124" />


```R
###################
# pcadapt outliers
###################
geno_tar <- read.pcadapt("/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_4pop.vcf", type = "vcf")
obj <- pcadapt(geno_tar, K = 2)

plot(obj, option = "manhattan")
pvals <- obj$pvalues
alpha <- 0.01  
outliers <- which(pvals < alpha)
length(outliers) #1644

# outlier info
out_meta <- vcf_tar@fix[outliers, c("CHROM", "POS", "ID")]
out_meta <- as.data.frame(out_meta, stringsAsFactors = FALSE)
out_meta$LocusID <- paste0(out_meta$CHROM, ":", out_meta$POS)

write.csv(out_meta, "/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_4pop_outlier_snps_meta.csv",
          row.names = FALSE)

# vcf with outliers only
vcf_tar_outliers <- vcf_tar
vcf_tar_outliers@fix <- vcf_tar@fix[outliers, , drop = FALSE]
vcf_tar_outliers@gt <- vcf_tar@gt[outliers, , drop = FALSE]

write.vcf(vcf_tar_outliers,
          "/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_4pop_outliers.vcf.gz")

# then unzip
system("gunzip -f /Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_4pop_outliers.vcf.gz")

# run SNMF on only outliers
vcf2geno("/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_4pop_outliers.vcf",
         "/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_4pop_outliers.geno")

project_outliers <- snmf("/Users/madelineeppley/Desktop/manta/tarapacana/tarapacana_4pop_outliers.geno",
                         K = 1:5, entropy = TRUE, repetitions = 10, project = "new")

print(summary(project_outliers))
plot(project_outliers, col = "blue", pch = 19, cex = 1.2, main = "Outlier SNPs: M. tarapacana")
```

<img width="1470" height="1296" alt="image" src="https://github.com/user-attachments/assets/17a32aec-c88d-4063-b7e7-db333624fe9f" />


