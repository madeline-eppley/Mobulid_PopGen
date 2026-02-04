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
