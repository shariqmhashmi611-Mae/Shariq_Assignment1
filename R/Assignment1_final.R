##***************************
## Assignment 1, Trichinella nativa comparing with Trichinella Spiralis
##
## Jordyn Anene  
## Luis Humberto Villasana Avila  
## Shariq Mumtaz Hashmi
## Sumayyah Khan
##
## 2026-11-02
##
## Secondary contributor Adhithya Dhivakar revised code
## Main edits made: 
##   1. Replaced the duplicated nativa/spiralis steps with two functions,
##      called load_te_data() and filter_dense_te(). The contig cutoff is now the
##      min_repeats argument. Also merged the two identical non-TE
##      lists into ls_no_te and updated the data paths to match the file
##      names in the repository.
##   2. Added run_wilcox() and a results table called df_wilcox that collects the
##      four Part 2 Wilcoxon tests with Holm-adjusted p-values.
##   3. Added a sensitivity analysis that runs the Part 2
##      tests again at cutoffs of 20-100 repeats and a figure.
##
## Revised 2026-11-23
##***************************

## Checking if this works in Github


## Cleaning the environment to make sure that the code can run on its own ----
rm(list = ls())

## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
library(Biostrings)
library(broom)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())


## Preparing CSV files ----
# Function to load a RepeatMasker CSV and fix the bracketed qry_left / in_repeat_end columns
load_te_data = function(path) {
  df = read.csv(path)
  df = df |>
    mutate(qry_left2 = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")), .after = qry_left) |>
    mutate(in_repeat_end2 = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")), .after = in_repeat_end)
  return(df)
}

# Original (unfiltered) tables for both species
df_te_nat_org   = load_te_data("../data/GCA_001447565.2.csv")
df_trich_te_org = load_te_data("../data/GCA_001447595.1.csv")



## Removing non-TE elements and classifying matching classes into either NA transposons or retrotransposons----

# List of non-TE elements (same list used for both species)
ls_no_te = c("Simple_repeat", "simple_repeat", "Low_complexity", "Satellite", "rRNA", "snRNA", "ARTEFACT", "Unknown")

# Updating datasets with only TEs
df_te_nat = df_te_nat_org |>
  filter(matching_class %notin% ls_no_te)

df_trich_te = df_trich_te_org |>
  filter(matching_class %notin% ls_no_te)

# Making sure there are not any non-TE left in Nativa
df_te_nat$matching_class %in% ls_no_te |>
  sum()
df_te_nat |>
  filter(matching_class %notin% ls_no_te) |> View()

# Making sure there are not any non-TE left in Spiralis
df_trich_te$matching_class %in% ls_no_te |>
  sum()

# To observe what are the TE elements present in Nativa
df_te_nat$matching_class |> unique()

# To observe what are the TE elements present in Spiralis
df_trich_te$matching_class |> unique()

# To identify how many unique qry_ids we have in Nativa
df_te_nat$qry_id |> unique() |> length()
#=============> we have 438 unique contigs

# To identify how many unique qry_ids we have in Nativa
df_trich_te$qry_id |> unique() |> length()
#=============> we have 315 unique contigs 

# To obtain a list of the dna transposons for T. nativa 
ls_dna_qry = df_te_nat |> 
  filter(str_starts(matching_class, "DNA/")) |>
  pull(qry_id) |> unique()
ls_dna_qry

# To obtain a list of the dna transposons for T. spiralis 
ls_ts_dna_qry = df_trich_te |>
  filter(str_starts(matching_class, "DNA/|RC")) |>
  pull(qry_id) |> unique()
ls_ts_dna_qry

# To obtain a list of a all the retrotransposons in Nativa
ls_retro_qry = df_te_nat |>
  filter(str_starts(matching_class, "PLE|PLE/|LTR/|LTR|LINE|LINE/|SINE/")) |>
  pull(qry_id) |> unique()
ls_retro_qry

# To obtain a list of a all the retrotransposons in Spiralis
ls_ts_retro_qry = df_trich_te |>
  filter(str_starts(matching_class, "PLE|PLE/|LTR/|LTR|LINE|LINE/|SINE/")) |>
  pull(qry_id) |> unique()
ls_ts_retro_qry

# Checking that the number of DNA transposons and retrotransposons match in Nativa
Biostrings::setdiff(ls_dna_qry, ls_retro_qry) |> length()
Biostrings::setdiff(ls_retro_qry, ls_dna_qry) |> length()

# Checking that the number of DNA transposons and retrotransposons match in Spiralis
Biostrings::setdiff(ls_ts_dna_qry, ls_ts_retro_qry) |> length()
Biostrings::setdiff(ls_ts_retro_qry, ls_ts_dna_qry) |> length()


## Creating different data frames according to the type of TE ----
#Data frame with only DNA transponsos in Nativa
df_dna_only = df_te_nat |>
  filter(qry_id %in% Biostrings::setdiff(ls_dna_qry, ls_retro_qry))
df_dna_only

#Data frame with only DNA transponsos in Spiralis
df_ts_dnaonly = df_trich_te |>
  filter(qry_id %in% Biostrings::setdiff(ls_ts_dna_qry, ls_ts_retro_qry))
df_ts_dnaonly

#Data frame with only retrotransponsos in Nativa
df_retro_only = df_te_nat |>
  filter(qry_id %in% Biostrings::setdiff(ls_retro_qry, ls_dna_qry))
df_retro_only

#Data frame with only retrotransponsos in spiralis
df_ts_retroonly = df_trich_te |>
  filter(qry_id %in% Biostrings::setdiff(ls_ts_retro_qry, ls_ts_dna_qry))
df_ts_retroonly

# Data frame with both DNA transponsos and retrotransponsos in Nativa
df_bothTE_nat = df_te_nat |>
  filter(qry_id %in% Biostrings::intersect(ls_dna_qry, ls_retro_qry))
df_bothTE_nat

# Data frame with both DNA transponsos and retrotransponsos in Spiralis
df_ts_both = df_trich_te |>
  filter(qry_id %in% Biostrings::intersect(ls_ts_dna_qry, ls_ts_retro_qry))
df_ts_both


## Finding the inflection point for large contigs ----

# For T.Spiralis
# Creating a plot to determine the inflection point to identify the size of large contigs in spiralis
df_trich_te_org |> 
  group_by(qry_id) |> 
  summarise(count = n()) |> 
  arrange(desc(count)) |> 
  pull(count) |> 
  plot()

# Using log10 to make more clear where the inflection point is located in spiralis
df_trich_te_org |> 
  group_by(qry_id) |> 
  summarise(count = n()) |> 
  arrange(desc(count)) |> 
  pull(count) |> 
  log10() |>
  plot()
#================> We determine that large contigs are above 50 in width


# Determine the inflection point to identify the size of large contigs in Nativa
df_te_nat_org |>
  group_by(qry_id) |>
  summarise(count = n()) |>
  arrange(desc(count)) |>
  pull(count) |>
  plot()

# Using log10 to make more clear where the inflection point is located in Nativa
df_te_nat_org |>
  group_by(qry_id) |>
  summarise(count = n()) |>
  arrange(desc(count)) |>
  pull(count) |>
  log10() |>
  plot()
#================> We determine that large contigs are above 50 in width

#================> Cutoff for repeat-dense contigs: >50 repeats per contig (both species).
#================> the curve bends from a gentle to a steep decline at ~50
#================> Repeats per contig (all repeat classes >>>> org tables) were sorted and plotted on a log10 scale;

# Function: keep TEs on repeat-dense contigs and label Species and TE_type
# df_org = the original, unfiltered table (contigs are counted before non-TEs are removed)
filter_dense_te = function(df_org, species, min_repeats = 50) {
  ls_dense = df_org |>
    group_by(qry_id) |>
    summarise(count = n()) |>
    arrange(desc(count)) |>
    filter(count > min_repeats) |>
    pull(qry_id)

  df_dense = df_org |>
    filter(matching_class %notin% ls_no_te) |>
    filter(qry_id %in% ls_dense) |>
    mutate(Species = species,
           TE_type = if_else(str_starts(matching_class, "DNA/|RC"), "DNA Transposons", "Retrotransposons"))
  return(df_dense)
}

df_dense_nat = filter_dense_te(df_te_nat_org, "T.nativa")
df_dense_ts  = filter_dense_te(df_trich_te_org, "T.spiralis")

# Number of contigs kept in each species (contigs that still have TE rows after non-TE removal)
df_dense_nat$qry_id |> unique() |> length()
#=================> 225 contigs (268 had >50 repeats before non-TE removal)
df_dense_ts$qry_id |> unique() |> length()
#=================> 186 contigs (220 had >50 repeats before non-TE removal)

#=================>  The same cutoff is used for both species so the comparison is consistent.
#=================>  Note: we haven't used FASTA files, contigs are ranked by repeat count, not by length, so "repeat dense" rather than "largest".

nrow(df_te_nat)
nrow(df_dense_nat)

nrow(df_trich_te)
nrow(df_dense_ts)

#====> Dense contigs hold 823 of 1319 TEs in T.Nativa and 838 of 1037 TEs in T.Spiralis.
  
# Identifying the amount of DNA Transposons and Retrotransposons in both Nativa and Spiralis to check values
df_dense_nat |>
  count(TE_type)

df_dense_ts |>
  count(TE_type)

df_dense_nat$qry_id |> unique()
df_dense_ts$qry_id |> unique()
  
df_dense_nat |>
  filter(TE_type == "Retrotransposons") |>
  count(matching_class)  

df_dense_ts |>
  filter(TE_type == "Retrotransposons") |>
  count(matching_class)

# Combine both species into one table
df_both_species = bind_rows(df_dense_nat, df_dense_ts)

nrow(df_both_species)
df_both_species |> count(Species, TE_type)


## Determine the timings of the TEs ----
# Do DNA transposons and retrotransposons differ in age, and is that pattern the same in T.nativa and T.spiralis?
# Do older TEs appear more decayed? Use the relationship between perc_div and perc_del

# Exploratory: Checking for TE divergence
df_both_species |> 
  ggplot(aes(x = perc_div, fill = TE_type)) +
  geom_histogram(binwidth = 1, color = "white") +
  labs(x = "Divergence from reference sequence %", y = "Number of TE copies", title = "Distribution of TE divergence by type and species") +
  facet_grid(TE_type ~ Species, scales = 'free_y') +
  theme_light() +
  scale_fill_brewer(palette = "Pastel2") +
  theme(legend.position = "none")

ggsave("../Output/Fig4_TE_Species_dist.png", width = 10, height = 6)


#==> T. nativa: no impossible values (all between ~2% and ~37% divergence); data passes the check.
#==> Retrotransposons: one main peak at ~20%; a few very old copies (~37%) separated from the main distribution,
#    and a thin tail of younger copies (2–8%) too small to indicate a recent burst. Outliers kept as genuine values.
#==> DNA transposons: patchy distribution (n = 85), a few younger copies at ~9–10%; gaps reflect small sample size.

# Plotting perc_div (substitution, as a proxy for age) against perc_del (deletion) for both species.
# To examine whether older TE copies carry more deletions, and whether this relationship differs between DNA transposons and Retrotransposon, and T.spiralis species and T.nativa species.
  
df_both_species |>
  ggplot(aes(x = perc_div, y = perc_del, colour = TE_type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = lm) +
  facet_wrap(~ Species) +
  labs(x = "Divergence from the reference sequence in %", y = "Deletion in %", color = "TE type")

ggsave("../Output/fig_div_vs_del.png", width = 10, height = 6)
  
# In both species, older retrotransposon copies show slightly more deletion, whereas DNA transposons show no clear relationship; the overlapping confidence bands indicate the differences between TE types and species are weak.

# Part 2
# Are DNA transposons and retrotransposons different ages, and is that difference the same in T.native and T.spiralis?

#Statistic analysis: Using wilcoxon ranked sum tested, based on the histogram we generated we noticed that the distibution of DNA TEs is regular. With wilcoxon ranked sum we will check whether the divergence differ between the two classes of TE (DNA vs Retrotransposon).
  
help("wilcox.test")

wilcox.test(perc_div ~ TE_type, data = df_both_species |> filter(Species == "T.nativa"))
#In T.nativa, divergence did not differ between DNA transposons and retrotransposons (Wilcoxon rank-sum test, W = 29899, p = 0.48; medians ≈ 21% for both).

wilcox.test(perc_div ~ TE_type, data = df_both_species |> filter(Species == "T.spiralis"))
    
# DNA transposons and retrotransposons did not differ in divergence in T.nativa, but DNA transposons were slightly more diverged in T.spiralis (W = 41606, p = 0.038). The T. spiralis difference is small (~1.4 percentage points), so the TE-type age pattern is broadly similar between species.      #=======> The small T.spiralis difference (p = 0.038) does not survive correction for running two tests, so it should be treated as a tentative trend, not a firm result.
      
#  Comparing divergence between species within each TE type, to test directly whether TE timing differs between T.nativa and T.spiralis.
      
wilcox.test(perc_div ~ Species, data = df_both_species |> filter(TE_type == "Retrotransposons"))

wilcox.test(perc_div ~ Species, data = df_both_species |> filter(TE_type == "DNA Transposons"))

# Run a Wilcoxon test and return the result as a one-row table
# group_col = the column to compare (e.g. "TE_type" or "Species")
run_wilcox = function(df, group_col, label) {
  test = wilcox.test(reformulate(group_col, response = "perc_div"), data = df)
  tidy(test) |>
    mutate(comparison = label) |>
    select(comparison, statistic, p.value)
}

# Four comparisons from Part 2, collected in one table
df_wilcox = bind_rows(
  run_wilcox(filter(df_both_species, Species == "T.nativa"), "TE_type", "T.nativa: DNA vs Retro"),
  run_wilcox(filter(df_both_species, Species == "T.spiralis"), "TE_type", "T.spiralis: DNA vs Retro"),
  run_wilcox(filter(df_both_species, TE_type == "Retrotransposons"), "Species", "Retro: nativa vs spiralis"),
  run_wilcox(filter(df_both_species, TE_type == "DNA Transposons"), "Species", "DNA: nativa vs spiralis")
) |>
  mutate(p_holm = p.adjust(p.value, method = "holm"))

df_wilcox

#==> After Holm correction across the four comparisons, none is significant (adjusted p = 0.119 at lowest), so the T. spiralis DNA vs retro difference is a tentative trend.
      
df_median_div = df_both_species |>
  group_by(Species, TE_type) |>
  summarise(med = median(perc_div), n = n())

## Sensitivity analysis: do the Part 2 conclusions depend on the 50-repeat cutoff? ----

# Same as run_wilcox(), but returns NA instead of an error if a group is too small at a strict cutoff
run_wilcox_safe = function(df, group_col, label) {
  tryCatch(run_wilcox(df, group_col, label),
           error = function(e) tibble(comparison = label, statistic = NA_real_, p.value = NA_real_))
}

# Rebuild the combined table at one cutoff and rerun the four Part 2 comparisons
run_cutoff = function(min_repeats) {
  df = bind_rows(
    filter_dense_te(df_te_nat_org, "T.nativa", min_repeats),
    filter_dense_te(df_trich_te_org, "T.spiralis", min_repeats)
  )
  bind_rows(
    run_wilcox_safe(filter(df, Species == "T.nativa"), "TE_type", "T.nativa: DNA vs Retro"),
    run_wilcox_safe(filter(df, Species == "T.spiralis"), "TE_type", "T.spiralis: DNA vs Retro"),
    run_wilcox_safe(filter(df, TE_type == "Retrotransposons"), "Species", "Retro: nativa vs spiralis"),
    run_wilcox_safe(filter(df, TE_type == "DNA Transposons"), "Species", "DNA: nativa vs spiralis")
  ) |>
    mutate(p_holm = p.adjust(p.value, method = "holm"),
           min_repeats = min_repeats,
           n_TEs = nrow(df),
           .before = 1)
}

df_sensitivity = map_dfr(c(20, 30, 50, 75, 100), run_cutoff)
df_sensitivity

# Adjusted p-values across cutoffs (dashed line = 0.05)
df_sensitivity |>
  ggplot(aes(x = min_repeats, y = p_holm, colour = comparison)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(x = "Minimum repeats per contig (cutoff)", y = "Holm-adjusted p-value",
       colour = "Comparison", title = "Sensitivity of TE divergence tests to the contig cutoff")

ggsave("../Output/fig_cutoff_sensitivity.png", width = 10, height = 6)

#==> No comparison is significant after Holm correction at any cutoff from 20 to 100 (lowest adjusted p = 0.119), so the Part 2 conclusions do not depend on the cutoff of 50. The raw p < 0.05 results at cutoff 50 do not hold at neighbouring cutoffs.

#=====> Medians were calculated per species and TE type to report the size of the age difference alongside the Wilcoxon p-values, since a p-value alone shows whether groups differ, not by how much. ANd to add captions to the Boxplot.
      
df_both_species |> filter(Species == "T.nativa") |>
  ggplot(aes(x = TE_type, y = perc_div, fill = TE_type)) +
  scale_fill_brewer(palette = "Pastel2") +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.5, colour = "grey30") +
  theme_light()
  
df_both_species |> filter(Species == "T.spiralis") |>
  ggplot(aes(x = TE_type, y = perc_div, fill = TE_type)) +
  scale_fill_brewer(palette = "Pastel2") +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.5, colour = "grey30") +
  theme_light()
  
df_both_species |>
  ggplot(aes(x = TE_type, y = perc_div, fill = TE_type)) +
  geom_violin(alpha = 0.4) +
  scale_fill_brewer(palette = "Pastel2") +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.1, alpha = 0.2, size = 0.8, colour = "grey30") +
  facet_wrap(~ Species) +
  labs(x = "TE type", y = "Divergence from reference sequence (%)", title = "TE divergence by type and species") +
  labs(caption = "T. nativa: DNA 21.2% (n = 85) vs Retro 21.3% (n = 738), Wilcoxon p = 0.48\nT. spiralis: DNA 22.8% (n = 100) vs Retro 21.4% (n = 738), Wilcoxon p = 0.038\nBetween species: Retro p = 0.76, DNA p = 0.030") +
    theme(legend.position = "none")
    
ggsave("../Output/fig_div_TEtype.png", width = 10, height = 6)
    
#Histograms showed irregular distributions for DNA transposons and unequal group sizes, so a non-parametric Wilcoxon test was used instead of a t-test.

#==> Median divergence is around 21% for both TE types in T.nativa, so DNA transposons and retrotransposons are of similar age. In T.spiralis, DNA transposons sit slightly higher (~23% vs ~21%), suggesting slightly older DNA TE copies.
#==> Retrotransposons look almost identical in both species (median ~21%, range ~3–38%), and span a wider range of ages than DNA TEs, i.e. activity over a longer period.
#==> All violins have a single bulge at ~20–25% and no peak at low divergence: no evidence of a recent burst of TE activity in either species; most copies are old.
#==> DNA groups are small (n = 85, 100) vs retros (n = 738 each).
  
# Part 3

#Do older TE copies tend to be shorter, and does that differ between DNA transposons and retrotransposons, and between T. nativa and T. spiralis?

df_both_species |>
  ggplot(aes(x = perc_div, y = qry_width, colour = TE_type)) +
  scale_y_log10() +
  geom_point(size = 0.8, alpha = 0.3) +
  geom_smooth(method = lm) +
  facet_grid(TE_type ~ Species, scales = "free_y") +
  labs(x = "Divergence from the reference sequence (%)",
       y = "TE copy length (bp, log scale)",
       title = "Older TE copies appear longer",
       caption = "Each point is one TE copy; lines are linear fits with 95% confidence bands.\nn: T. nativa DNA = 85, Retro = 738; T. spiralis DNA = 100, Retro = 738") +
  theme_light() +
  theme(legend.position = "none")
    
ggsave("../Output/fig_div_TEwidth_4.png", width = 10, height = 8)
    
#==> Most TE copies are short fragments (<200 bp), far below full length TE size; consistent with old, decayed TEs. Contrary to expectation, copy length increases with divergence, clearly for retrotransposons (both species, n = 738), weakly for DNA transposons (wide confidence bands, n = 85, 100).
#==> Likely bias: highly diverged copies must be longer to pass RepeatMasker's score threshold, so old short fragments are under-detected. Length vs divergence therefore reflects detectability as well as decay.
    
# Part 4: Do contigs in T.nativa and T.spiralis differ in their composition of TE families?
# PCA on TE family copy counts per contig (one row per contig, one column per family),
# coloured by species, to see whether contigs from the two species carry similar or different TE family mixes.

# Creating a wide table to run PCA
    
df_pca_wide = df_both_species |>
  group_by(qry_id, Species, matching_class) |>
  summarise(count = n()) |>
  pivot_wider(names_from = matching_class,
              values_from = count,
              values_fill = 0)

dim(df_pca_wide)

help("prcomp")
    
pca_te = df_pca_wide |>
  ungroup() |>
  select(-qry_id, -Species) |>
  prcomp(scale. = TRUE)

summary(pca_te)

# # Table of PC scores for each contig (from pca_te$x), with species added, for plotting PC1 vs PC2.

help("as_tibble")

df_pca_score = pca_te$x |>
  as_tibble() |>
  mutate(Species = df_pca_wide$Species)

df_pca_score
    
#Plotting
df_pca_score |>
  ggplot(aes(x = PC1, y = PC2, colour = Species)) +
  geom_point(alpha = 0.5) +
stat_ellipse(level = 0.95, segments = 51, type = "t") +
  labs(x = "PC1 = 8.6% of variance", y = "PC2 = 7.1% of variance", title = "Contigs from both species have similar TE composition") +
  labs( caption = "Each point is one contig; ellipses cover ~95% of each species contigs.\nPCA on scaled TE family counts (23 families).") +
  theme_light()
  
pca_te$rotation[ ,1:2]
  
ggsave("../Output/PCA_species.png", width = 8, height = 6)
