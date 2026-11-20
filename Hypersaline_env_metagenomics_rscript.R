######################################
#LOAD LIBRARIES
######################################
library(tidyverse)
library(vegan)
library(reshape2)
#######################################
#Read Metadata
#######################################
metadata <- read.csv(
  "metadata.csv",
  stringsAsFactors = FALSE
)

head(metadata)
######################################
#Read all bracken files
#####################################
files <- list.files(
  pattern = "*.tabular",
  full.names = TRUE
)
#Combine all samples
all_data <- lapply(files, function(x){
  
  df <- read.delim(
    x,
    header = TRUE
  )
  
  df$Sample <- gsub(
    ".tabular",
    "",
    basename(x),
    fixed = TRUE
  )
  
  return(df)
  
}) %>%
  bind_rows()
# Add Environment Information
all_data <- left_join(
  all_data,
  metadata,
  by = "Sample"
)
#Create Species-level dataset
species_data <- all_data %>%
  filter(taxonomy_lvl == "S")
#Create Sample x Species Matrix
species_matrix <- species_data %>%
  select(
    Sample,
    name,
    fraction_total_reads
  ) %>%
  pivot_wider(
    names_from = name,
    values_from = fraction_total_reads,
    values_fill = 0
  )
#Create Community Matrix
community <- species_matrix[, -1]

rownames(community) <- species_matrix$Sample
dim(community)

##############################
#Unique Species Plot
###############################
#Check truly unique species plot
species_env <- all_data %>%
  group_by(name, Environment) %>%
  summarise(
    abundance = mean(fraction_total_reads),
    .groups = "drop"
  )

species_occurrence <- species_env %>%
  group_by(name) %>%
  summarise(
    n_env = n()
  )

unique_species <- species_env %>%
  inner_join(
    species_occurrence %>%
      filter(n_env == 1),
    by = "name"
  )
#Inspect
unique_species %>%
  count(Environment)
unique_species %>%
  arrange(Environment, desc(abundance))
#Generate top 10 unique species per environment
top_unique_species <- unique_species %>%
  group_by(Environment) %>%
  slice_max(
    abundance,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()
#Plot
library(ggplot2)

library(ggplot2)

ggplot(
  top_unique_species,
  aes(
    x = reorder(name, abundance),
    y = abundance,
    fill = Environment
  )
) +
  
  geom_col() +
  
  coord_flip() +
  
  facet_wrap(
    ~Environment,
    scales = "free_y",
    ncol = 2
  ) +
  
  theme_bw() +
  
  labs(
    title = "Environment-Specific Microbial Signature Taxa",
    x = "Species",
    y = "Mean Relative Abundance"
  ) +
  
  theme(
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 10,
      colour = "black",
      face = "plain"
    ),
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    
    legend.position = "none"
  )

##########################################
#Alpha Diversity Plot And Wilcoxon test
###########################################
#Create sample x species mtrix
library(dplyr)
library(tidyr)

species_matrix <- all_data %>%
  select(Sample, name, fraction_total_reads) %>%
  pivot_wider(
    names_from = name,
    values_from = fraction_total_reads,
    values_fill = 0
  )
#Create community Matrix
community <- species_matrix[, -1]

rownames(community) <- species_matrix$Sample
#calculate Aplha Diversity
library(vegan)

shannon <- diversity(
  community,
  index = "shannon"
)

simpson <- diversity(
  community,
  index = "simpson"
)

richness <- specnumber(community)
#Create Diversity Dataframe
alpha_div <- data.frame(
  Sample = rownames(community),
  Shannon = shannon,
  Simpson = simpson,
  Richness = richness
)
#Add Environment Information
metadata <- metadata %>%
  distinct(Sample, Environment)

alpha_div <- left_join(
  alpha_div,
  metadata,
  by = "Sample"
)
head(alpha_div)
#Plot Shannon Diversity
library(ggplot2)

ggplot(alpha_div,
       aes(
         x = Environment,
         y = Shannon,
         fill = Environment
       )) +
  
  geom_boxplot(
    alpha = 0.8
  ) +
  
  geom_jitter(
    width = 0.1,
    size = 2
  ) +
  
  theme_bw(base_size = 14) +
  
  labs(
    x = "",
    y = "Shannon Diversity Index"
  ) +
  
  theme(
    legend.position = "none"
  )
#Kruskal-Wallis test
kruskal.test(
  Shannon ~ Environment,
  data = alpha_div
)
########################################
#BETA DIVERSITY PLOT
########################################
#=========================================
# PCoA ANALYSIS
#=========================================

library(vegan)
library(ggplot2)
library(dplyr)

# Species abundance matrix
species_matrix <- all_data %>%
  filter(taxonomy_lvl == "S") %>%
  select(Sample, name, fraction_total_reads) %>%
  tidyr::pivot_wider(
    names_from = name,
    values_from = fraction_total_reads,
    values_fill = 0
  )

species_mat <- as.data.frame(species_matrix)

rownames(species_mat) <- species_mat$Sample
species_mat$Sample <- NULL

# Match metadata order
metadata2 <- metadata[
  match(rownames(species_mat), metadata$Sample),
]

# Bray-Curtis distance
bray_dist <- vegdist(
  species_mat,
  method = "bray"
)

# PCoA
pcoa <- cmdscale(
  bray_dist,
  k = 2,
  eig = TRUE
)

# Variance explained
var_exp <- round(
  100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]),
  2
)

# Create dataframe
pcoa_df <- data.frame(
  Sample = rownames(species_mat),
  PCoA1 = pcoa$points[,1],
  PCoA2 = pcoa$points[,2]
)

pcoa_df <- left_join(
  pcoa_df,
  metadata2,
  by = "Sample"
)

# Plot
ggplot(
  pcoa_df,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = Environment
  )
) +
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  stat_ellipse(
    aes(fill = Environment),
    geom = "polygon",
    alpha = 0.15,
    colour = NA
  ) +
  labs(
    title = "PCoA of Hypersaline Microbial Communities",
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)")
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )
ggsave(
  "Figurebeta_diversity.tiff",
  width = 8,
  height = 6,
  dpi = 600,
  compression = "lzw"
)
#PERMANOVA
adonis2(
  bray_dist ~ Environment,
  data = metadata,
  permutations = 999
)
####################################
#SPEARMAN CORRELATION ANALYSIS fromTop 20 species
####################################
#Calculate mean species abundance per environment
library(dplyr)

top20_species <- all_data %>%
  group_by(name) %>%
  summarise(
    MeanAbundance = mean(fraction_total_reads),
    .groups = "drop"
  ) %>%
  arrange(desc(MeanAbundance)) %>%
  slice(1:20)

top20_species
#Create Environment × top 20 Species matrix
species_env_top20 <- all_data %>%
  filter(name %in% top20_species$name) %>%
  group_by(Environment, name) %>%
  summarise(
    MeanAbundance = mean(fraction_total_reads),
    .groups = "drop"
  )
#COnvert to Matrix
library(tidyr)

top20_matrix <- species_env_top20 %>%
  pivot_wider(
    names_from = name,
    values_from = MeanAbundance,
    values_fill = 0
  )

env_names <- top20_matrix$Environment

top20_matrix <- as.data.frame(top20_matrix)

rownames(top20_matrix) <- env_names

top20_matrix$Environment <- NULL
#Calculate Spearman correlation matrix
cor_top20 <- cor(
  t(as.matrix(top20_matrix)),
  method = "spearman"
)
#Generate heatmap
library(pheatmap)
pheatmap(
  cor_top20,
  
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  
  clustering_method = "complete",
  
  color = colorRampPalette(
    c("steelblue","white","red")
  )(100),
  
  fontsize = 12,
  fontsize_row = 12,
  fontsize_col = 12,
  
  main = "Spearman Correlation Between Environments",
  
  border_color = NA
)
##################################
#Top microbial Species Across environment plot
#########################################
#Remove Homo sapiens
all_data_clean <- all_data %>%
  filter(name != "Homo sapiens")
#Calculate mean abundance per environment
library(dplyr)

species_env <- all_data_clean %>%
  group_by(Environment, name) %>%
  summarise(
    abundance = mean(fraction_total_reads),
    .groups = "drop"
  )
#Select top 20 species overall
top20_species <- species_env %>%
  group_by(name) %>%
  summarise(
    total_abundance = sum(abundance)
  ) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:20)
#Keep only these species
top_species_env <- species_env %>%
  filter(name %in% top20_species$name)
#Create heatmap matrix
library(tidyr)

heatmap_data <- top_species_env %>%
  pivot_wider(
    names_from = name,
    values_from = abundance,
    values_fill = 0
  )
#Convert to matrix
heatmap_matrix <- as.data.frame(heatmap_data)

rownames(heatmap_matrix) <- heatmap_matrix$Environment

heatmap_matrix$Environment <- NULL

heatmap_matrix <- as.matrix(heatmap_matrix)
#Environment annotation
annotation_row <- data.frame(
  Group = rownames(heatmap_matrix)
)

rownames(annotation_row) <- rownames(heatmap_matrix)
#PLOT
library(pheatmap)

pheatmap(
  heatmap_matrix,
  
  scale = "row",
  
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  annotation_row = annotation_row,
  
  fontsize_row = 12,
  fontsize_col = 10,
  
  color = colorRampPalette(
    c("steelblue", "ivory", "orangered")
  )(100),
  
  border_color = "grey70",
  
  main = "Top Microbial Species Across Different Environments"
)
#########################################################
## COG CATEGORY STACKED BARPLOT
#########################################################

library(tidyverse)
library(ggplot2)
library(vegan)
#Read the COG excel file
library(readxl)

cog <- read_excel("cog_category.xlsx")
#Check sample names
colnames(cog)
#Concert the matrix into one long table
library(tidyr)
library(dplyr)

long <- cog %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "COG"
  )
#remove unassigned COGs
long <- long %>%
  filter(
    !is.na(COG),
    COG != "-"
  )
#Split multiple COG categories
library(tidyr)

long <- long %>%
  mutate(COG = strsplit(COG, "")) %>%
  unnest(COG)
head(long)
#Add Metadata
metadata <- read.csv("metadata.csv", stringsAsFactors = FALSE)
long <- left_join(long, metadata, by = "Sample")
head(long)
#########################################################
## STEP 1: Remove Function Unknown (S)
#########################################################

long2 <- long %>%
  filter(COG != "S")

#########################################################
## STEP 2: Count each COG category
#########################################################

count_data <- long2 %>%
  group_by(Sample, Environment, COG) %>%
  summarise(Count = n(), .groups = "drop")

#########################################################
## STEP 3: Calculate Relative Abundance
#########################################################

count_data <- count_data %>%
  group_by(Sample) %>%
  mutate(RelAbundance = Count/sum(Count)*100)

#########################################################
## STEP 4: Average across samples
#########################################################

plot_data <- count_data %>%
  group_by(Environment, COG) %>%
  summarise(RelAbundance = mean(RelAbundance),
            .groups="drop")

#########################################################
## STEP 5: Order COG Categories
#########################################################

plot_data$COG <- factor(
  plot_data$COG,
  levels=c(
    "J","A","K","L","B",
    "D","Y","V","T","M",
    "N","Z","W","U","O",
    "C","G","E","F","H",
    "I","P","Q"
  )
)

#########################################################
## STEP 6: Plot
#########################################################

p <- ggplot(plot_data,
            aes(Environment,
                RelAbundance,
                fill=COG))+
  
  geom_bar(stat="identity",
           position="fill",
           width=0.8)+
  
  scale_y_continuous(
    labels=scales::percent_format(scale=100)
  )+
  
  labs(
    title="Functional Composition Based on COG Categories",
    x="Environment",
    y="Relative Abundance (%)",
    fill="COG"
  )+
  
  theme_bw(base_family="Times New Roman")+
  
  theme(
    
    plot.title=element_text(
      face="bold",
      size=18,
      hjust=0.5
    ),
    
    axis.title=element_text(
      face="bold",
      size=16
    ),
    
    axis.text=element_text(
      size=13,
      colour="black"
    ),
    
    legend.title=element_text(
      face="bold",
      size=14
    ),
    
    legend.text=element_text(
      size=12
    ),
    
    panel.grid.major=element_line(colour="grey90"),
    panel.grid.minor=element_blank()
    
  )

p
ggsave(
  "Figure_COG.tiff",
  plot=p,
  width=8,
  height=6,
  units="in",
  dpi=600,
  compression="lzw"
)
#############################################################################
#PHYLUM PLOT (four environments)
#############################################################################
#Read All phylum files
library(dplyr)
library(readr)
library(stringr)

files <- list.files(
  "C:/Users/Lenovo/OneDrive/Documents/Saline Soil Analysis/bracken/phylum_data",
  pattern = "phylum.tabular$",
  full.names = TRUE
)

phylum_data <- lapply(files, function(file){
  
  df <- read.delim(file)
  
  sample_name <- tools::file_path_sans_ext(basename(file))
  
  df$Sample <- sample_name
  
  if(grepl("Atacama", sample_name)){
    df$Environment <- "Atacama Desert"
  }
  
  if(grepl("Nebraska", sample_name)){
    df$Environment <- "Nebraska Salt Marsh"
  }
  
  if(grepl("Rann", sample_name)){
    df$Environment <- "Rann of Kutch"
  }
  
  if(grepl("Xinjiang", sample_name)){
    df$Environment <- "Xinjiang Saline Soil"
  }
  
  return(df)
  
}) %>%
  bind_rows()
#Remove non-microbial taxa
phylum_data <- phylum_data %>%
  filter(
    !name %in% c(
      "Annelida",
      "Apicomplexa",
      "Arthropoda",
      "Ascomycota",
      "Bacillariophyta",
      "Basidiomycota",
      "Chlorophyta",
      "Cnidaria",
      "Echinodermata",
      "Euglenida",
      "Glomeromycota",
      "Mollusca",
      "Nematoda",
      "Phaeophyceae",
      "Platyhelminthes",
      "Porifera",
      "Streptophyta",
      "Preplasmiviricota",
      "Uroviricota"
    )
  )
#Calculate Mean abundance per environment
phylum_summary <- phylum_data %>%
  group_by(Environment, name) %>%
  summarise(
    MeanAbundance = mean(fraction_total_reads),
    .groups = "drop"
  )
#Select top 10 phyla globally
top10_phyla <- phylum_summary %>%
  group_by(name) %>%
  summarise(
    TotalAbundance = sum(MeanAbundance)
  ) %>%
  arrange(desc(TotalAbundance)) %>%
  slice(1:10)
#prepare plotting data
plot_data <- phylum_summary %>%
  filter(name %in% top10_phyla$name)
#Convert to percentages
plot_data <- plot_data %>%
  group_by(Environment) %>%
  mutate(
    RelativeAbundance =
      MeanAbundance /
      sum(MeanAbundance) * 100
  )
#Create Plot
library(ggplot2)

ggplot(
  plot_data,
  aes(
    x = Environment,
    y = RelativeAbundance,
    fill = name
  )
) +
  
  geom_bar(
    stat = "identity",
    width = 0.8
  ) +
  
  theme_bw() +
  
  labs(
    title = "Taxonomic Composition at Phylum Level",
    x = "",
    y = "Relative Abundance (%)",
    fill = "Phylum"
  ) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      colour = "black"
    ),
    
    legend.title = element_text(
      face = "bold"
    )
  )
#GENUS PLOT
#Read All genus files
library(dplyr)
library(readr)
library(stringr)

files <- list.files(
  "C:/Users/Lenovo/OneDrive/Documents/Saline Soil Analysis/bracken/genus_data",
  pattern = "genus.tabular$",
  full.names = TRUE
)

genus_data <- lapply(files, function(file){
  
  df <- read.delim(file)
  
  sample_name <- tools::file_path_sans_ext(basename(file))
  
  df$Sample <- sample_name
  
  if(grepl("Atacama", sample_name)){
    df$Environment <- "Atacama Desert"
  }
  
  if(grepl("Nebraska", sample_name)){
    df$Environment <- "Nebraska Salt Marsh"
  }
  
  if(grepl("Rann", sample_name)){
    df$Environment <- "Rann of Kutch"
  }
  
  if(grepl("Xinjiang", sample_name)){
    df$Environment <- "Xinjiang Saline Soil"
  }
  
  return(df)
  
}) %>%
  bind_rows()
#Remove non-microbial taxa
genus_data <- genus_data %>%
  filter(
    !name %in% c(
      "Homo"
    )
  )
#Calculate Mean abundance per environment
genus_summary <- genus_data %>%
  group_by(Environment, name) %>%
  summarise(
    MeanAbundance = mean(fraction_total_reads),
    .groups = "drop"
  )
#Select top 10 genera globally
top10_genus <- genus_summary %>%
  group_by(name) %>%
  summarise(
    TotalAbundance = sum(MeanAbundance)
  ) %>%
  arrange(desc(TotalAbundance)) %>%
  slice(1:10)
#prepare plotting data
plot_data <- genus_summary %>%
  filter(name %in% top10_genus$name)
#Convert to percentages
plot_data <- plot_data %>%
  group_by(Environment) %>%
  mutate(
    RelativeAbundance =
      MeanAbundance /
      sum(MeanAbundance) * 100
  )
#Create Plot
library(ggplot2)

ggplot(
  plot_data,
  aes(
    x = Environment,
    y = RelativeAbundance,
    fill = name
  )
) +
  
  geom_bar(
    stat = "identity",
    width = 0.8
  ) +
  
  theme_bw() +
  
  labs(
    title = "Taxonomic Composition at Genus Level",
    x = "",
    y = "Relative Abundance (%)",
    fill = "Genus"
  ) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      colour = "black"
    ),
    
    legend.title = element_text(
      face = "bold"
    )
  )
#############################################################
#FUNCTIONAL PCoA PLOT
############################################################
#Read the Excel file
library(readxl)
cog <- read_xlsx("cog_category.xlsx")
#Convert to long format
long <- cog %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "COG"
  )
#Remove missing values
long <- long %>%
  filter(
    !is.na(COG),
    COG != "-"
  )
#Split multiple COG categories
long <- long %>%
  mutate(COG = strsplit(COG, "")) %>%
  unnest(COG)
#Read Metadata
metadata <- read.csv("metadata.csv",
                     stringsAsFactors = FALSE)
#Merge Metadata
long <- left_join(
  long,
  metadata,
  by = "Sample"
)
#Count COG categories
count_data <- long %>%
  group_by(
    Sample,
    Environment,
    COG
  ) %>%
  summarise(
    Count = n(),
    .groups = "drop"
  )
#Calculate relative abundance within each sample
count_data <- count_data %>%
  group_by(Sample) %>%
  mutate(
    RelativeAbundance =
      Count / sum(Count)
  ) %>%
  ungroup()
#Create Functional matrix
functional_matrix <- count_data %>%
  select(
    Sample,
    COG,
    RelativeAbundance
  ) %>%
  pivot_wider(
    names_from = COG,
    values_from = RelativeAbundance,
    values_fill = 0
  )
#Add metadata
functional_metadata <- metadata %>%
  filter(
    Sample %in%
      functional_matrix$Sample
  )

functional_matrix <- left_join(
  functional_metadata,
  functional_matrix,
  by = "Sample"
)
#Bray-Curtis distance
bray_dist <- vegdist(
  functional_matrix[, -(1:2)],
  method = "bray"
)
#PERMANOVA
adonis2(
  bray_dist ~ Environment,
  data = functional_metadata,
  permutations = 999
)
#PCoA
pcoa <- cmdscale(
  bray_dist,
  eig = TRUE,
  k = 2
)
#Percentage explained
variance <- round(
  100 *
    pcoa$eig /
    sum(pcoa$eig[pcoa$eig > 0]),
  2
)
#Plot data
pcoa_df <- data.frame(
  Sample = functional_matrix$Sample,
  Environment = functional_matrix$Environment,
  PCoA1 = pcoa$points[,1],
  PCoA2 = pcoa$points[,2]
)
#Plot
p <- ggplot(
  pcoa_df,
  aes(
    PCoA1,
    PCoA2,
    colour = Environment
  )
) +
  
  geom_point(size = 3.5) +
  
  stat_ellipse(
    type = "t",
    linewidth = 0.8
  ) +
  
  theme_bw(base_family = "Times New Roman") +
  
  labs(
    title = "Principal Coordinate Analysis of Functional Profiles",
    x = paste0(
      "PCoA1 (",
      variance[1],
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      variance[2],
      "%)"
    )
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    axis.text = element_text(
      size = 12,
      colour = "black"
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

p
ggsave(
  "Figure_PCoA_functional.tiff",
  plot=p,
  width=8,
  height=6,
  units="in",
  dpi=600,
  compression="lzw"
)
###############################################################
##AMR_metadata and amr table merging
###############################################################
library(readr)
library(dplyr)

amr <- read_tsv("AMR_master_table.tsv")

metadata <- read_csv("metadata.csv")

amr_final <- left_join(amr, metadata, by = "Sample")
#Import metadata and amr table 
library(readr)
library(dplyr)

# Read AMR table
amr <- read_tsv("AMR_master_table.tsv")

# Read metadata
metadata <- read_csv("metadata.csv")

# Merge by Sample
amr <- left_join(amr, metadata, by = "Sample")

# Check
head(amr)
str(amr)
#################################################
#Basic Summary
#################################################
nrow(amr)                      # Total ARG hits

length(unique(amr$Sample))     # Samples

length(unique(amr$`Gene symbol`))

length(unique(amr$Class))

length(unique(amr$Subclass))
################################################
#ARG Abundance per sample
################################################
sample_counts <-
  amr %>%
  count(Sample)

################################################
#Antibiotic Class distribution and plot
################################################
rm(list=ls())

##############################
# Packages
##############################

library(tidyverse)
library(viridis)
library(RColorBrewer)
library(pheatmap)

##############################
# Read data
##############################

amr <- read.delim("AMR_master_table.tsv",
                  sep="\t",
                  header=TRUE,
                  stringsAsFactors=FALSE)

metadata <- read.csv("metadata.csv",
                     stringsAsFactors=FALSE)

##############################
# Merge metadata
##############################

amr <- left_join(amr, metadata, by="Sample")

##############################
# Theme
##############################

theme_pub <- theme_bw(base_size=14)+
  theme(
    axis.text=element_text(size=12,color="black"),
    axis.title=element_text(size=14,face="bold"),
    legend.title=element_text(face="bold"),
    legend.text=element_text(size=11),
    strip.background=element_rect(fill="white"),
    panel.grid.minor=element_blank(),
    panel.grid.major=element_line(colour="grey90"),
    plot.title=element_text(face="bold",hjust=0.5,size=16)
  )

class_counts <-
  amr %>%
  count(Class)


################################################
#ARG in Environments
################################################
env_class <-
  amr %>%
  count(Environment,
        Class)

p3 <-
  ggplot(env_class,
         aes(Environment,
             n,
             fill=Class))+
  
  geom_bar(stat="identity")+
  
  labs(title="ARG Classes across Environments",
       x="Environment",
       y="ARG Hits")+
  
  theme_pub+
  
  theme(axis.text.x=
          element_text(angle=0,
                       hjust=0.5))

p3

ggsave("Figure3_Environment_ARGclasses.tiff",
       p3,
       dpi=600,
       width=9,
       height=6)