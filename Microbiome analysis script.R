# ============================================================
# Supplementary R Script: Microbiome Analysis
# ============================================================

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(ggrepel)
library(tidyverse)
library(UpSetR)
# ------------------------------------------------------------
# Function to read Bracken output
# ------------------------------------------------------------
read_bracken <- function(file, sample_name) {
  df <- read.table(file, header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
  
  colnames(df) <- c("fraction_total_reads", "new_est_reads", "added_reads",
                    "taxonomy_lvl", "taxonomy_id", "name")
  
  df <- df %>%
    filter(taxonomy_lvl == "G") %>%
    filter(!name %in% c("root", "cellular organisms")) %>%
    mutate(Sample = sample_name)
  
  return(df)
}

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
AF <- read_bracken("AF_report_bracken_genuses.txt", "AF")
AM <- read_bracken("AM_report_bracken_genuses.txt", "AM")
NAF <- read_bracken("NAF_report_bracken_genuses.txt", "NAF")
NAM <- read_bracken("NAM_report_bracken_genuses.txt", "NAM")

all_data <- bind_rows(AF, AM, NAF, NAM)

# ------------------------------------------------------------
# Calculate relative abundance
# ------------------------------------------------------------
all_data <- all_data %>%
  group_by(Sample) %>%
  mutate(RelAbundance = (new_est_reads / sum(new_est_reads)) * 100)

# ------------------------------------------------------------
# Metadata
# ------------------------------------------------------------
metadata <- data.frame(
  Sample = c("AF", "AM", "NAF", "NAM"),
  Group = c("Tobacco", "Tobacco", "Non-Tobacco", "Non-Tobacco")
)

# ------------------------------------------------------------
# Create abundance matrix
# ------------------------------------------------------------
abundance_matrix <- all_data %>%
  select(Sample, name, RelAbundance) %>%
  pivot_wider(names_from = name, values_from = RelAbundance, values_fill = 0)

# Convert to dataframe
abundance_matrix <- as.data.frame(abundance_matrix)

# Set correct rownames
rownames(abundance_matrix) <- abundance_matrix$Sample

# Remove Sample column
abundance_matrix$Sample <- NULL

# ------------------------------------------------------------
# Alpha Diversity (Shannon) & Wilcoxon test
# ------------------------------------------------------------
shannon <- diversity(as.matrix(abundance_matrix), index = "shannon")
alpha_df <- data.frame(
  Sample = rownames(abundance_matrix),
  Shannon = shannon
)
alpha_df <- merge(alpha_df, metadata, by = "Sample")
alpha_df
table(alpha_df$Group)

# Wilcoxon test
wilcox.test(Shannon ~ Group, data = alpha_df)

# Plot
ggplot(alpha_df, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 3) +
  theme_classic()

# ------------------------------------------------------------
# Beta Diversity (Bray-Curtis + PCoA + PERMANOVA)
# ------------------------------------------------------------
#PCoA plot data input:
AF  <- read.table("AF_genus.txt", header=TRUE, sep="\t", fill=TRUE)
AM  <- read.table("AM_genus.txt", header=TRUE, sep="\t", fill=TRUE)
NAF <- read.table("NAF_genus.txt", header=TRUE, sep="\t", fill=TRUE)
NAM <- read.table("NAM_genus.txt", header=TRUE, sep="\t", fill=TRUE)
#Select Columns
AF  <- AF[, c("name","new_est_reads")]
AM  <- AM[, c("name","new_est_reads")]
NAF <- NAF[, c("name","new_est_reads")]
NAM <- NAM[, c("name","new_est_reads")]
#Rename
colnames(AF)=c("Taxa","AF")
colnames(AM)=c("Taxa","AM")
colnames(NAF)=c("Taxa","NAF")
colnames(NAM)=c("Taxa","NAM")
#Merge
merged <- Reduce(function(x,y) merge(x,y,by="Taxa",all=TRUE),
                 list(AF,AM,NAF,NAM))
#Replace and Save
merged[is.na(merged)] <- 0
write.csv(merged,"merged_genus_matrix.csv",row.names=FALSE)

# Step 1: Bray-Curtis distance
otu <- merged[, -1]
rownames(otu) <- merged$Taxa
otu <- t(otu)
bray <- vegdist(otu, method = "bray")

# Step 2: PCoA
pcoa <- cmdscale(bray, eig = TRUE, k = 2)

# Step 3: Convert to dataframe
pcoa_df <- as.data.frame(pcoa$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")

# Step 4: Add sample names
pcoa_df$Sample <- rownames(pcoa_df)

# Step 5: Define groups (IMPORTANT: correct order)
pcoa_df$Group <- c("Tobacco", "Tobacco", "Non-Tobacco", "Non-Tobacco")
# Order = AF, AM, NAF, NAM

# Step 6: Plot
ggplot(pcoa_df, aes(PCoA1, PCoA2, color = Group)) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = Sample), size = 5) +
  scale_color_manual(values = c("blue", "red")) +
  theme_classic(base_size = 16) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "PCoA based on Bray-Curtis Dissimilarity",
    x = "PCoA1",
    y = "PCoA2"
  )

# PERMANOVA
adonis2(bray ~ Group, data = metadata, permutations = 999)

# ------------------------------------------------------------
# Differential Abundance (Log2FC)
# ------------------------------------------------------------
abund <- read.csv("relative_abundance_genus.csv", row.names = 1)
abund_filtered <- abund[apply(abund, 1, max) > 1, ]
mean_tob <- rowMeans(abund_filtered[, c("AF","AM")])

mean_non <- rowMeans(abund_filtered[, c("NAF","NAM")])

fc <- log2((mean_tob + 0.001) / (mean_non + 0.001))

fc_df <- data.frame(
  Taxa = rownames(abund_filtered),
  Log2FC = fc
)
top_pos <- fc_df %>%
  arrange(desc(Log2FC)) %>%
  head(7)

top_neg <- fc_df %>%
  arrange(Log2FC) %>%
  head(7)

fc_final <- rbind(top_pos, top_neg)

ggplot(fc_final,
       aes(x=reorder(Taxa, Log2FC),
           y=Log2FC,
           fill=Log2FC>0))+
  geom_bar(stat="identity")+
  coord_flip()+
  scale_fill_manual(values=c("#E64B35","#4DBBD5"))+
  theme_classic(base_size=14)+
  theme(legend.position="none")


# ------------------------------------------------------------
# SIMPER Analysis
# ------------------------------------------------------------
simper_res <- simper(abundance_matrix, metadata$Group)

simper_df <- summary(simper_res)[[1]]
simper_df <- as.data.frame(simper_df)

simper_top <- simper_df %>%
  arrange(desc(average)) %>%
  slice(1:10)

# Plot SIMPER
ggplot(simper_top, aes(x = reorder(rownames(simper_top), average), y = average)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_classic()

# Save outputs
write.csv(simper_top, "simper_top_taxa.csv")

# ------------------------------------------------------------
# RELATIVE ABUNDANCE GENUS PLOT
# ------------------------------------------------------------
# Read genus abundance table
AF  <- read.table("AF_genus.txt", header=TRUE, sep="\t")
AM  <- read.table("AM_genus.txt", header=TRUE, sep="\t")
NAF <- read.table("NAF_genus.txt", header=TRUE, sep="\t")
NAM <- read.table("NAM_genus.txt", header=TRUE, sep="\t")

# Keep only name + fraction columns
AF  <- AF[, c("name", "fraction_total_reads")]
AM  <- AM[, c("name", "fraction_total_reads")]
NAF <- NAF[, c("name", "fraction_total_reads")]
NAM <- NAM[, c("name", "fraction_total_reads")]

# Rename abundance columns
colnames(AF)[2]  <- "AF"
colnames(AM)[2]  <- "AM"
colnames(NAF)[2] <- "NAF"
colnames(NAM)[2] <- "NAM"

# Merge all samples
merged <- reduce(list(AF, AM, NAF, NAM), full_join, by="name")
merged[is.na(merged)] <- 0

# Select Top 10 Genera
top10 <- merged %>%
  mutate(mean_abundance = rowMeans(select(., AF:NAM))) %>%
  arrange(desc(mean_abundance)) %>%
  slice(1:10)

# Convert to long format
melted <- top10 %>%
  pivot_longer(cols=AF:NAM,
               names_to="Sample",
               values_to="RelativeAbundance")

# Plot
ggplot(melted, aes(x=Sample, y=RelativeAbundance*100, fill=name)) +
  geom_bar(stat="identity") +
  labs(title="Top 10 Genera Relative Abundance",
       x="Sample",
       y="Relative Abundance (%)",
       fill="Genera") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust=0.4, face="bold"),
    axis.text.x = element_text(face="bold"),
    legend.position = "right"
  )
# ------------------------------------------------------------
# RELATIVE ABUNDANCE SPECIES PLOT
# ------------------------------------------------------------
# Read species abundance table
AF  <- read.table("AF_species.txt", header=TRUE, sep="\t")
AM  <- read.table("AM_species.txt", header=TRUE, sep="\t")
NAF <- read.table("NAF_species.txt", header=TRUE, sep="\t")
NAM <- read.table("NAM_species.txt", header=TRUE, sep="\t")

# Keep only name + fraction columns
AF  <- AF[, c("name", "fraction_total_reads")]
AM  <- AM[, c("name", "fraction_total_reads")]
NAF <- NAF[, c("name", "fraction_total_reads")]
NAM <- NAM[, c("name", "fraction_total_reads")]

# Rename abundance columns
colnames(AF)[2]  <- "AF"
colnames(AM)[2]  <- "AM"
colnames(NAF)[2] <- "NAF"
colnames(NAM)[2] <- "NAM"

# Merge all samples
merged <- reduce(list(AF, AM, NAF, NAM), full_join, by="name")
merged[is.na(merged)] <- 0

# Select Top 10 Species
top10 <- merged %>%
  mutate(mean_abundance = rowMeans(select(., AF:NAM))) %>%
  arrange(desc(mean_abundance)) %>%
  slice(1:10)

# Convert to long format
melted <- top10 %>%
  pivot_longer(cols=AF:NAM,
               names_to="Sample",
               values_to="RelativeAbundance")

# Plot
ggplot(melted, aes(x=Sample, y=RelativeAbundance*100, fill=name)) +
  geom_bar(stat="identity") +
  labs(title="Top 10 Species Relative Abundance",
       x="Sample",
       y="Relative Abundance (%)",
       fill="Species") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust=0.4, face="bold"),
    axis.text.x = element_text(face="bold"),
    legend.position = "right"
  )
# ------------------------------------------------------------
# RELATIVE ABUNDANCE PHYLUM PLOT
# ------------------------------------------------------------
# Read Phyla abundance table
AF  <- read.table("AF_phylum.txt", header=TRUE, sep="\t")
AM  <- read.table("AM_phylum.txt", header=TRUE, sep="\t")
NAF <- read.table("NAF_phylum.txt", header=TRUE, sep="\t")
NAM <- read.table("NAM_phylum.txt", header=TRUE, sep="\t")

# Keep only name + fraction columns
AF  <- AF[, c("name", "fraction_total_reads")]
AM  <- AM[, c("name", "fraction_total_reads")]
NAF <- NAF[, c("name", "fraction_total_reads")]
NAM <- NAM[, c("name", "fraction_total_reads")]

# Rename abundance columns
colnames(AF)[2]  <- "AF"
colnames(AM)[2]  <- "AM"
colnames(NAF)[2] <- "NAF"
colnames(NAM)[2] <- "NAM"

# Merge all samples
merged <- reduce(list(AF, AM, NAF, NAM), full_join, by="name")
merged[is.na(merged)] <- 0

# Select Top 10 Phyla
top10 <- merged %>%
  mutate(mean_abundance = rowMeans(select(., AF:NAM))) %>%
  arrange(desc(mean_abundance)) %>%
  slice(1:10)

# Convert to long format
melted <- top10 %>%
  pivot_longer(cols=AF:NAM,
               names_to="Sample",
               values_to="RelativeAbundance")

# Plot
ggplot(melted, aes(x=Sample, y=RelativeAbundance*100, fill=name)) +
  geom_bar(stat="identity") +
  labs(title="Top 10 Phyla Relative Abundance",
       x="Sample",
       y="Relative Abundance (%)",
       fill="Phylum") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust=0.4, face="bold"),
    axis.text.x = element_text(face="bold"),
    legend.position = "right"
  )
# ------------------------------------------------------------
# RELATIVE ABUNDANCE FAMILY PLOT
# ------------------------------------------------------------
# Read species abundance table
AF  <- read.table("AF_family.txt", header=TRUE, sep="\t")
AM  <- read.table("AM_family.txt", header=TRUE, sep="\t")
NAF <- read.table("NAF_family.txt", header=TRUE, sep="\t")
NAM <- read.table("NAM_family.txt", header=TRUE, sep="\t")

# Keep only name + fraction columns
AF  <- AF[, c("name", "fraction_total_reads")]
AM  <- AM[, c("name", "fraction_total_reads")]
NAF <- NAF[, c("name", "fraction_total_reads")]
NAM <- NAM[, c("name", "fraction_total_reads")]

# Rename abundance columns
colnames(AF)[2]  <- "AF"
colnames(AM)[2]  <- "AM"
colnames(NAF)[2] <- "NAF"
colnames(NAM)[2] <- "NAM"

# Merge all samples
merged <- reduce(list(AF, AM, NAF, NAM), full_join, by="name")
merged[is.na(merged)] <- 0

# Select Top 10 Family
top10 <- merged %>%
  mutate(mean_abundance = rowMeans(select(., AF:NAM))) %>%
  arrange(desc(mean_abundance)) %>%
  slice(1:10)

# Convert to long format
melted <- top10 %>%
  pivot_longer(cols=AF:NAM,
               names_to="Sample",
               values_to="RelativeAbundance")

# Plot
ggplot(melted, aes(x=Sample, y=RelativeAbundance*100, fill=name)) +
  geom_bar(stat="identity") +
  labs(title="Top 10 Family Relative Abundance",
       x="Sample",
       y="Relative Abundance (%)",
       fill="Family") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust=0.4, face="bold"),
    axis.text.x = element_text(face="bold"),
    legend.position = "right"
  )
# ------------------------------------------------------------
# UpSet PLOT
# ------------------------------------------------------------

# Convert abundance to presence/absence
binary <- as.data.frame((abund > 0) * 1)

# Generate UpSet plot
upset(binary,
      sets = c("AF","AM","NAF","NAM"),
      keep.order = TRUE,
      order.by = "freq",
      main.bar.color = "black",
      sets.bar.color = "black")