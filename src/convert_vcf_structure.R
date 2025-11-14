r
library(adegenet)
library(vcfR)

# Read VCF file
vcf <- read.vcfR("input.vcf")

# Convert to genlight object
gl <- vcfR2genlight(vcf)

# Convert to STRUCTURE format
genind_obj <- gi2genpop(gl)
# Export to STRUCTURE format using genind2structure() or write.structure()
