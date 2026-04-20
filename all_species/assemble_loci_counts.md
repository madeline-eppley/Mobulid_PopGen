## STACKS loci param opt counts figures
### for all species

To scp the files from the `/opt` directory on explorer to local  `/Desktop/manta` folder, I used this, but changed the species name directory for each species.

```bash
scp eppley.m@login.explorer.northeastern.edu:/projects/gatins/2025_Mobulid/mobular/opt/loci_counts.csv ~/Desktop/manta
```

### then in R, I made one script to create all supp. figure plots

```R
## assemble STACKS parameter space supp figure plots
# M. Eppley, v1 uploaded to github 4/20/26

# libs
library(ggplot2)
library(dplyr)

# one general function to apply to all species since output format is the same
create_param_plot <- function(csv_path, species_name, optimal_M, optimal_n, output_path) {
  df <- read.csv(csv_path)
  df$M <- as.numeric(gsub("m3_M(\\d+)_n\\d+", "\\1", df$params))
  df$n <- as.numeric(gsub("m3_M\\d+_n(\\d+)", "\\1", df$params))
  # optimal loci count
  optimal_loci <- df %>% filter(M == optimal_M, n == optimal_n) %>% pull(loci)
  # plot
  p <- ggplot(df, aes(x=factor(M), y=loci, fill=factor(n))) +
    geom_bar(stat="identity", position="dodge", width=0.7, color="black", linewidth=0.3) +
    geom_text(aes(label=format(loci, big.mark=",")), 
              position=position_dodge(width=0.7), 
              vjust=-0.5, size=2.5) +
    scale_fill_brewer(palette="Blues", name="n value") +
    labs(title=paste0(species_name, " STACKS parameter space (r80 loci)"),
         x="M: mismatches within individuals", 
         y="Number of polymorphic loci (r80)") +
    theme_classic() +
    theme(
      plot.title = element_text(face="bold", size=14),
      axis.text = element_text(size=11),
      axis.title = element_text(size=12),
      legend.position = "right"
    ) +
    scale_y_continuous(labels=scales::comma, limits=c(0, max(df$loci)*1.1))
  
  ggsave(output_path, p, width=8, height=5, dpi=300)
  cat("\n=== ", species_name, " ===\n")
  cat("Optimal parameters: M=", optimal_M, ", n=", optimal_n, " with ", 
      format(optimal_loci, big.mark=","), " loci\n", sep="")
  df %>% arrange(desc(loci)) %>% head() %>% print()
  return(p)}

# plot for tarapacana
tarapacana_plot <- create_param_plot(
  csv_path = "~/Desktop/manta/tarapacana_loci_counts.csv",
  species_name = "M. tarapacana",
  optimal_M = 2,
  optimal_n = 2,
  output_path = "~/Desktop/manta/tarapacana_param_opt.png")

# plot for birostris
birostris_plot <- create_param_plot(
  csv_path = "~/Desktop/manta/birostris_loci_counts.csv",
  species_name = "M. birostris",
  optimal_M = 2,
  optimal_n = 2,
  output_path = "~/Desktop/manta/birostris_param_opt.png")

# plot for munkiana
munkiana_plot <- create_param_plot(
  csv_path = "~/Desktop/manta/munkiana_loci_counts.csv",
  species_name = "M. munkiana",
  optimal_M = 3,
  optimal_n = 4,
  output_path = "~/Desktop/manta/munkiana_param_opt.png")

# plot for thurstoni
thurstoni_plot <- create_param_plot(
  csv_path = "~/Desktop/manta/thurstoni_loci_counts.csv",
  species_name = "M. thurstoni",
  optimal_M = 2,
  optimal_n = 3,
  output_path = "~/Desktop/manta/thurstoni_param_opt.png")

# plot for mobular
mobular_plot <- create_param_plot(
  csv_path = "~/Desktop/manta/mobular_loci_counts.csv",
  species_name = "M. mobular",
  optimal_M = 2,
  optimal_n = 2,
  output_path = "~/Desktop/manta/mobular_param_opt.png")

# display check
print(tarapacana_plot)
print(birostris_plot)
print(munkiana_plot)
print(thurstoni_plot)
print(mobular_plot)
```
