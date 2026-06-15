###############################################################################
#Inter and intraindividual variation analysis
###############################################################################
#Load Required Packages
install.packages("BiocManager")
BiocManager::install("HMP16SData")
library(HMP16SData)
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(tidyr)
set.seed(123)  #for reproducibility
#Load HMP V3 -V5 dataset:
se <- V35()
#Convert to phyloseq:
ps <- as_phyloseq(se)
ps
#Check the available body subsite names:
unique(sample_data(ps)$HMP_BODY_SUBSITE)
#Selection of 5 subsites
target_sites <- c(
  "Saliva",
  "Stool",
  "Anterior Nares",
  "Right Antecubital Fossa",
  "Vaginal Introitus"
)

ps_sub <- subset_samples(
  ps,
  HMP_BODY_SUBSITE %in% target_sites
)

ps_sub
#Keep only individuals with ≥2 visits:
meta <- as.data.frame(sample_data(ps_sub))

rsid_counts <- table(meta$RSID)

multi_ids <- names(rsid_counts[rsid_counts >= 2])

ps_multi <- subset_samples(
  ps_sub,
  RSID %in% multi_ids
)

ps_multi
#Find individuals common to ALL 5 subsites
meta <- as.data.frame(sample_data(ps_multi))

rsid_list <- meta |>
  dplyr::group_by(HMP_BODY_SUBSITE) |>
  dplyr::summarise(rsids = list(unique(RSID)), .groups = "drop")
rsid_list
#Take the intersection:
common_rsids <- Reduce(intersect, rsid_list$rsids)
length(common_rsids)
#Subset to ONLY those common individuals
ps_common <- subset_samples(
  ps_multi,
  RSID %in% common_rsids
)
#Verify
sample_data(ps_common) |>
  as.data.frame() |>
  dplyr::group_by(HMP_BODY_SUBSITE) |>
  dplyr::summarise(n_individuals = dplyr::n_distinct(RSID))
#Bray–Curtis distance
bray_common <- phyloseq::distance(ps_common, method = "bray")
#PERMANOVA — individual discrimination
meta_df <- data.frame(sample_data(ps_common))
class(meta_df)
head(meta_df$RSID)

adonis_individual <- adonis2(
  bray_common ~ RSID,
  data = meta_df,
  permutations = 999
)

adonis_individual
#PERMANOVA per subsite
library(vegan)

sites <- c(
  "Saliva",
  "Stool",
  "Anterior Nares",
  "Right Antecubital Fossa",
  "Vaginal Introitus"
)

permanova_site <- data.frame()

for (site in sites) {
  
  ps_site <- subset_samples(
    ps_common,
    HMP_BODY_SUBSITE == site
  )
  
  # Bray–Curtis for this subsite
  bray_site <- phyloseq::distance(ps_site, method = "bray")
  
  # Metadata (FORCED base data.frame)
  meta_site <- data.frame(sample_data(ps_site))
  
  ad <- adonis2(
    bray_site ~ RSID,
    data = meta_site,
    permutations = 999
  )
  
  permanova_site <- rbind(
    permanova_site,
    data.frame(
      Site = site,
      R2 = ad$R2[1],
      p_value = ad$`Pr(>F)`[1]
    )
  )
}

permanova_site
#Inter-individual distances
compute_inter <- function(ps_obj) {
  
  bray <- phyloseq::distance(ps_obj, "bray")
  bray_mat <- as.matrix(bray)
  meta <- data.frame(sample_data(ps_obj))
  
  rsids <- unique(meta$RSID)
  
  inter <- c()
  
  for (i in 1:(length(rsids) - 1)) {
    for (j in (i + 1):length(rsids)) {
      
      s1 <- rownames(meta[meta$RSID == rsids[i], ])
      s2 <- rownames(meta[meta$RSID == rsids[j], ])
      
      d <- bray_mat[s1, s2]
      inter <- c(inter, as.vector(d))
    }
  }
  
  data.frame(InterDist = inter)
}
#Inter-individual distances for all subsites
sites <- c(
  "Saliva",
  "Stool",
  "Anterior Nares",
  "Right Antecubital Fossa",
  "Vaginal Introitus"
)

inter_all <- data.frame()

for (site in sites) {
  
  ps_site <- subset_samples(
    ps_common,
    HMP_BODY_SUBSITE == site
  )
  
  inter_site <- compute_inter(ps_site)
  inter_site$Site <- site
  
  inter_all <- rbind(inter_all, inter_site)
}
#check
table(inter_all$Site)
summary(inter_all$InterDist)
#Inter-individual Variation plot: 
ggplot(inter_all,
       aes(x = Site, y = InterDist, fill = Site)) +
  geom_boxplot(outlier.alpha = 0.25) +
  theme_bw() +
  labs(
    title = "Inter-individual microbiome variation across forensic-relevant body subsites",
    y = "Inter-individual Bray–Curtis distance",
    x = "Body subsite"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

#Intra-individual distance:
compute_intra <- function(ps_obj) {
  
  bray <- phyloseq::distance(ps_obj, "bray")
  bray_mat <- as.matrix(bray)
  meta <- data.frame(sample_data(ps_obj))
  
  rsids <- unique(meta$RSID)
  intra <- c()
  
  for (id in rsids) {
    samples <- rownames(meta[meta$RSID == id, ])
    
    if (length(samples) > 1) {
      d <- bray_mat[samples, samples]
      intra <- c(intra, d[upper.tri(d)])
    }
  }
  
  data.frame(IntraDist = intra)
}
#Compute intra-individual distances for all subsites:
sites <- c(
  "Saliva",
  "Stool",
  "Anterior Nares",
  "Right Antecubital Fossa",
  "Vaginal Introitus"
)

intra_all <- data.frame()

for (site in sites) {
  
  ps_site <- subset_samples(
    ps_common,
    HMP_BODY_SUBSITE == site
  )
  
  intra_site <- compute_intra(ps_site)
  intra_site$Site <- site
  
  intra_all <- rbind(intra_all, intra_site)
}
table(intra_all$Site)
summary(intra_all$IntraDist)
#(plot)Temporal stability of microbiome across body sites.
library(ggplot2)

ggplot(intra_all,
       aes(x = Site, y = IntraDist, fill = Site)) +
  geom_boxplot(outlier.alpha = 0.25) +
  theme_bw() +
  labs(
    title = "Temporal stability of microbiomes across body subsites",
    y = "Intra-individual Bray–Curtis distance",
    x = "Body subsite"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

#Intra and inter-individual Plot:
combined_df <- rbind(
  data.frame(
    Distance = intra_all$IntraDist,
    Type = "Intra-individual",
    Site = intra_all$Site
  ),
  data.frame(
    Distance = inter_all$InterDist,
    Type = "Inter-individual",
    Site = inter_all$Site
  )
)

ggplot(combined_df,
       aes(x = Site, y = Distance, fill = Type)) +
  geom_boxplot(position = position_dodge(0.8),
               outlier.alpha = 0.3) +
  theme_bw() +
  labs(
    title = "Intra- vs inter-individual microbiome variation across body subsites",
    y = "Bray–Curtis distance",
    x = "Body subsite"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  "Figure_intra.png",
  width = 10,
  height = 5,
  dpi = 600
)

############################################################
# Genus-level composition of the salivary microbiome (Inter-individual genus level)
# (49 unique individuals from the common forensic cohort)
############################################################

library(phyloseq)
library(ggplot2)

# Saliva samples from the common 49-individual cohort
ps_saliva <- subset_samples(
  ps_common,
  HMP_BODY_SUBSITE == "Saliva"
)

# Verify
nsamples(ps_saliva)
length(unique(sample_data(ps_saliva)$RSID))

# Metadata
meta <- data.frame(sample_data(ps_saliva))
meta$SampleID <- sample_names(ps_saliva)

# Retain one sample per individual
# Visit 1 preferred; Visit 2 used if Visit 1 unavailable

meta_unique <- meta[order(meta$RSID, meta$VISITNO), ]
meta_unique <- meta_unique[!duplicated(meta_unique$RSID), ]

# Verify
nrow(meta_unique)
length(unique(meta_unique$RSID))
table(meta_unique$VISITNO)

# Create 49-individual phyloseq object
ps_saliva_49 <- prune_samples(
  meta_unique$SampleID,
  ps_saliva
)

# Verify
nsamples(ps_saliva_49)
length(unique(sample_data(ps_saliva_49)$RSID))

# Collapse to genus level
ps_genus_49 <- tax_glom(
  ps_saliva_49,
  taxrank = "GENUS"
)

# Relative abundance
ps_rel_49 <- transform_sample_counts(
  ps_genus_49,
  function(x) x / sum(x)
)

# Top 10 genera
genus_means <- sort(
  taxa_sums(ps_rel_49) / nsamples(ps_rel_49),
  decreasing = TRUE
)

top10_taxa <- names(genus_means)[1:10]

# Save for Figure 2
top10_taxa_49 <- top10_taxa

# Keep top 10 genera
ps_top10_49 <- prune_taxa(
  top10_taxa_49,
  ps_rel_49
)

# Convert for plotting
df_49 <- psmelt(ps_top10_49)

df_49$Sample <- factor(
  df_49$Sample,
  levels = unique(df_49$Sample)
)

# Plot
p1 <- ggplot(
  df_49,
  aes(
    x = Sample,
    y = Abundance,
    fill = GENUS
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.9
  ) +
  labs(
    title = "Genus-level composition of the salivary microbiome",
    x = "Individuals",
    y = "Relative abundance"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(face = "bold")
  )

p1
#Repeat for other subsites to generate the plot 
############################################################
# Salivary microbiome composition in male and female participants
############################################################

library(phyloseq)
library(ggplot2)

# Create genus list from Figure 1
top10_genera <- unique(df_49$GENUS)

# Full saliva dataset
ps_saliva_full <- subset_samples(
  ps,
  HMP_BODY_SUBSITE == "Saliva"
)

# Metadata
meta_full <- data.frame(sample_data(ps_saliva_full))

set.seed(123)

# Select 87 females and 87 males
female_rsids <- unique(
  meta_full$RSID[
    meta_full$SEX == "Female"
  ]
)

male_rsids <- unique(
  meta_full$RSID[
    meta_full$SEX == "Male"
  ]
)

sel_female <- sample(
  female_rsids,
  87
)

sel_male <- sample(
  male_rsids,
  87
)

selected_rsids <- c(
  sel_female,
  sel_male
)

# Keep selected participants
ps_saliva_sex <- subset_samples(
  ps_saliva_full,
  RSID %in% selected_rsids
)

# One sample per participant
meta_sex <- data.frame(sample_data(ps_saliva_sex))
meta_sex$SampleID <- sample_names(ps_saliva_sex)

meta_sex <- meta_sex[
  order(meta_sex$RSID, meta_sex$VISITNO),
]

meta_sex <- meta_sex[
  !duplicated(meta_sex$RSID),
]

ps_saliva_sex174 <- prune_samples(
  meta_sex$SampleID,
  ps_saliva_sex
)

# Verify
nsamples(ps_saliva_sex174)
table(sample_data(ps_saliva_sex174)$SEX)

# Genus level
ps_genus_sex <- tax_glom(
  ps_saliva_sex174,
  taxrank = "GENUS"
)

# Relative abundance
ps_rel_sex <- transform_sample_counts(
  ps_genus_sex,
  function(x) x / sum(x)
)

# Convert to dataframe
df_sex <- psmelt(ps_rel_sex)

# Keep only the genera present in Figure 1
df_sex <- subset(
  df_sex,
  GENUS %in% top10_genera
)

# Preserve genus order from Figure 1
df_sex$GENUS <- factor(
  df_sex$GENUS,
  levels = top10_genera
)

# Preserve sample order
df_sex$Sample <- factor(
  df_sex$Sample,
  levels = unique(df_sex$Sample)
)

# Plot
p2 <- ggplot(
  df_sex,
  aes(
    x = Sample,
    y = Abundance,
    fill = GENUS
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.9
  ) +
  facet_wrap(
    ~SEX,
    scales = "free_x",
    nrow = 1,
    labeller = as_labeller(
      c(
        Female = "(A) Female",
        Male = "(B) Male"
      )
    )
  ) +
  labs(
    title = "Salivary microbiome composition in male and female participants",
    x = "Individuals",
    y = "Relative abundance"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(face = "bold")
  )

p2
###############################################################################
#WILCOXON RESULTS
###############################################################################
sites <- unique(inter_all$Site)

wilcox_results <- data.frame()

for(site in sites){
  
  inter_vals <- subset(inter_all, Site == site)$InterDist
  intra_vals <- subset(intra_all, Site == site)$IntraDist
  
  wt <- wilcox.test(inter_vals, intra_vals)
  
  wilcox_results <- rbind(
    wilcox_results,
    data.frame(
      Site = site,
      P_value = wt$p.value,
      Median_Inter = median(inter_vals),
      Median_Intra = median(intra_vals)
    )
  )
}

wilcox_results
# Abundance for all 49 individuals
ps_saliva_sex <- subset_samples(
  ps_saliva,
  !is.na(SEX)
)
