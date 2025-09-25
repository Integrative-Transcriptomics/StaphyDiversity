rm(list = ls())
.libPaths( c( .libPaths(), "/home/wittepaz/R/x86_64-pc-linux-gnu-library/4.4/") )
library(ape)
library(tidyr)
library(dplyr)
library(ggtree)
library(ggnewscale)
library(ggplot2)

output_path <- "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/results/repo/StaphyDiversity/figures/Figure2/output/"
input_path <- "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/results/repo/StaphyDiversity/figures/Figure2/input/"

tree <- read.tree(paste(input_path, "no_paralogs_NJ_2.nwk", sep = ""))
antismash_overview <- read.csv(paste(input_path, "summary_with_clusterblast_with_percentage_clustercompare.csv", sep = ""), header = TRUE, sep = ";")

gtdb_data <- read.csv(paste(input_path, "gtdb-search.tsv", sep = ""), header = TRUE, sep = "\t")
# Split s__ and get the last part
gtdb_data$gtdb_species_name <- sapply(gtdb_data$gtdb_taxonomy, function(x) strsplit(x, "s__")[[1]][2])
gtdb_data <- subset(gtdb_data, select = c("accession", "gtdb_species_name"))

head(gtdb_data)
# gtdb_data$gtdb_species_name <- apply(gtdb_data, function(x) print(x))
sum(duplicated(gtdb_data$gtdb_species_name))
# trim every entry
antismash_overview <- as.data.frame(lapply(antismash_overview, function(x) trimws(x, which = c("both"))))

# for the alst three columns, we need to convert the values to numeric
antismash_overview$staphyloferrin.A <-as.numeric(antismash_overview$staphyloferrin.A)
antismash_overview$staphyloferrin.B <- as.numeric(antismash_overview$staphyloferrin.B)
antismash_overview$staphylopine<- as.numeric(antismash_overview$staphylopine)
antismash_overview$Others <- as.numeric(antismash_overview$unknown)

# merge gtdb_data with antismash_overview 
antismash_overview <- merge(antismash_overview, gtdb_data, by.x = "X", by.y = "accession", all.x = TRUE)

sorted_antismash_overview <- antismash_overview %>% arrange(number_of_records, staphyloferrin.A, staphyloferrin.B, staphylopine, Others, gtdb_species_name)
export_antismash_overview <- sorted_antismash_overview[c("gtdb_species_name", "number_of_records", "staphyloferrin.A", "staphyloferrin.B", "staphylopine", "Others")]

# replace Others column by nothing if 0, otherwise write True
export_antismash_overview$Others <- ifelse(export_antismash_overview$Others == 0, "", "True")
# sum over number of records
sum(as.numeric(export_antismash_overview$number_of_records))

antismash_overview_subset <- antismash_overview[c("gtdb_species_name","staphyloferrin.A", "staphyloferrin.B", "staphylopine")]
antismash_metadata <- antismash_overview_subset


rownames(antismash_metadata) <- antismash_metadata$gtdb_species_name

antismash_metadata <- antismash_metadata[,-1]

# replace complete if 100, incomplete if > 0, none if 0
antismash_metadata[antismash_metadata > 0 & antismash_metadata < 100] <- "Incomplete BGC"
antismash_metadata[antismash_metadata == 100] <- "Complete BGC"
antismash_metadata[antismash_metadata == 0] <- "Not found"

# get species with staphyloferrin.A incomplete
staphyloferrin_A_incomplete <- antismash_overview[antismash_overview$staphyloferrin.A < 100 & antismash_overview$staphyloferrin.A > 0,]
# remove NA values
staphyloferrin_A_incomplete <- staphyloferrin_A_incomplete[!is.na(staphyloferrin_A_incomplete$gtdb_species_name),]
# from the columns gtdb_species_name and genome create a named list
named_list <- setNames(as.list(antismash_overview$gtdb_species_name), antismash_overview$X)

# from antimsash_overview$genome remove the ending with .[0-9] and create a named list
antismash_overview$X <- gsub("\\.[0-9]", "", antismash_overview$X)
named_list <- setNames(as.list(antismash_overview$gtdb_species_name), antismash_overview$X)
rename_tree <- function(tree, named_list){
  for (i in 1:length(tree$tip.label)){
    tree_tip <- tree$tip.label[i]
    if (tree_tip %in% names(named_list)){
      tree$tip.label[i] <- named_list[[tree_tip]]
    }
  }
  return(tree)
}

# Create a reverse named_list
reverse_named_list <- setNames(as.list(antismash_overview$X), antismash_overview$gtdb_species_name)

# for Staphylococcus lugdunensis, manual inspection showed that the BGC was not identified as Staphyloferrin A but as other. 
# So set to Incomplete BGC for staphyloferrin.A
antismash_metadata["Staphylococcus lugdunensis", "staphyloferrin.A"] <- "Incomplete BGC"
# for Staphylococcus ureilyticus, manual inspection of a further accession code showed a complete BGC, so set to Complete BGC
antismash_metadata["Staphylococcus ureilyticus", "staphyloferrin.A"] <- "Complete BGC"

# rename tree
tree <- rename_tree(tree, named_list)
# root tree with Duncaniella sp023665865
tree_viz <- ggtree(root(tree,outgroup = "Duncaniella sp023665865", edgelabel = TRUE) ) + geom_tiplab(size = 2.5) 

colnames(antismash_metadata) <- c("sfa", "sbn", "cnt")

antismash_view <- gheatmap(tree_viz,
                antismash_metadata,  colnames_position = "top", font.size=2,  offset=0.0175,width=.15, colnames_angle=90,colnames_offset_y =4.5 ) +
                scale_fill_manual(values = c("darkblue", "lightblue", "white")) + guides(fill=guide_legend(title="ClusterCompare results"))


# For separation of the new heatmap
temp1 <- antismash_view + new_scale_fill()

# how many columns have a value > 0
antismash_overview$positive_bgcs <- rowSums(antismash_overview[,c("staphyloferrin.A", "staphyloferrin.B", "staphylopine", "Others")] > 0, na.rm = TRUE) 
antismash_overview$Others_mod <- ifelse(antismash_overview$positive_bgcs == antismash_overview$number_of_records, antismash_overview$Others, antismash_overview$Others + as.numeric(antismash_overview$number_of_records) - antismash_overview$positive_bgcs)
unknown_bgc <- antismash_overview[c("gtdb_species_name","Others_mod")]
rownames(unknown_bgc) <- unknown_bgc$gtdb_species_name
unknown_bgc <- subset(unknown_bgc, select = -c(gtdb_species_name))
# if unknown_bgc is 0, set to False
unknown_bgc[unknown_bgc > 0] <- "True"
unknown_bgc[unknown_bgc == 0] <- "False"

# For Staphylococcus lugdunensis, set Other to false since it was manually inspected and found to be Staphyloferrin A
unknown_bgc["Staphylococcus lugdunensis", "Others_mod"] <- "False"
# For Staphylococcus ureilyticus, set to False as it was also Staphyloferrin A after manual inspection
unknown_bgc["Staphylococcus ureilyticus", "Others_mod"] <- "False"

colnames(unknown_bgc) <- c("others")


temp0 <- gheatmap(temp1,unknown_bgc, colnames_position = "top", font.size=2,offset=0.0315, width=.05,
         colnames_angle=90, colnames_offset_y = 3.75) +
                scale_fill_manual(values = c("white", "darkgrey"), labels = c("", "Found")) + guides(fill=guide_legend(title="Other BGCs"))

tmp_0 <- temp0 + new_scale_fill()


mmseqs_overview <- read.csv(paste(input_path, "homology_results_staphyloccus_adapted_with_cntA_new.tsv", sep = ""), header = FALSE, sep = "\t")
colnames(mmseqs_overview) <- c("lipoprotein", "hit", "hit_locus_tag", "hit_genome_id", "bits", "identity", "evalue", "start_alignment", "end_alingment", "length_alignment", "differences")
# set identity to numeric
mmseqs_overview$identity <- as.numeric(as.character(mmseqs_overview$identity))
# if duplicate hit_locus_tag, keep the one with the highest identity
mmseqs_overview <- mmseqs_overview[order(mmseqs_overview$identity, decreasing = TRUE),]
mmseqs_overview <- mmseqs_overview[!duplicated(mmseqs_overview$hit_locus_tag),]


summarize_lipoprotein  <- function(lipoprotein){
  lipoprotein_hits <- mmseqs_overview[mmseqs_overview$lipoprotein == lipoprotein,]
  lipoprotein_hits <- lipoprotein_hits[lipoprotein_hits$identity > 0.5,]
  # group and get number of hits
  lipoprotein_hits_grouped <- lipoprotein_hits %>% group_by(hit_genome_id) %>% summarise(numHits = n(), homology = identity[which.max(identity)]) %>% arrange(desc(homology))
  #rename homology column to homology<lipoprotein>
  colnames(lipoprotein_hits_grouped)[3] <- paste0(lipoprotein, "Homology")
  colnames(lipoprotein_hits_grouped)[2] <- paste0(lipoprotein, "Hits")
  return(lipoprotein_hits_grouped)
}

sirA_hits_mmseqs <- summarize_lipoprotein("SirA")
htsA_hits_mmseqs <- summarize_lipoprotein("HtsA")
fhuD1_hits_mmseqs <- summarize_lipoprotein("FhuD1")
fhuD2_hits_mmseqs <- summarize_lipoprotein("FhuD2")
sstD_hits_mmseqs <- summarize_lipoprotein("SstD")
cntA_hits_mmseqs <- summarize_lipoprotein("CntA")


# join all dataframes
mmseqs_combined <- full_join(sirA_hits_mmseqs, htsA_hits_mmseqs, by = "hit_genome_id")
mmseqs_combined <- full_join(mmseqs_combined, fhuD1_hits_mmseqs, by = "hit_genome_id")
mmseqs_combined <- full_join(mmseqs_combined, fhuD2_hits_mmseqs, by = "hit_genome_id")
mmseqs_combined <- full_join(mmseqs_combined, sstD_hits_mmseqs, by = "hit_genome_id")
mmseqs_combined <- full_join(mmseqs_combined, cntA_hits_mmseqs, by = "hit_genome_id")

# from mmseqs_combined$hit_genome_id remove the ending with .[0-9] 
mmseqs_combined$hit_genome_id <- gsub("\\.[0-9]", "", mmseqs_combined$hit_genome_id)

# # create genome name column via named list
mmseqs_combined$genome <- named_list[mmseqs_combined$hit_genome_id]

 # get only columns numHitsHtsA and numHitsSirA
mmseqs_combined_df <- as.data.frame(mmseqs_combined[c("HtsAHomology", "SirAHomology", "CntAHomology", "FhuD1Homology", "FhuD2Homology", "SstDHomology" )], row.names = mmseqs_combined$genome)
rownames(mmseqs_combined_df) <- mmseqs_combined$genome
mmseqs_combined_df[is.na(mmseqs_combined_df)] <- 0

# Add missing genomes from antismash_overview
missing_genomes <- setdiff(antismash_overview$gtdb_species_name, rownames(mmseqs_combined_df))
# Add missing genomes with 0 values
for (genome in missing_genomes) {
  mmseqs_combined_df[genome, ] <- 0
}

# visualize with previous tree
library(ggnewscale)
colnames(mmseqs_combined_df) <- c("htsA", "sirA", "cntA", "fhuD1", "fhuD2", "sstD")
viz_mmseqs <- gheatmap(tmp_0,mmseqs_combined_df, low="white", high="darkgreen", colnames_position = "top", font.size=2,offset=0.0375, width=.2,
         colnames_angle=90, colnames_offset_y = 3.75) +
    scale_fill_gradient2(low="white", mid="white", high="darkgreen", na.value="grey", midpoint=0.5,name="Similarity to reference gene", breaks=c(0,0.5,0.75,1.0), minor_breaks=c(0,0.5,0.75)) +
         theme(plot.title = element_text(hjust = 0.5))

temp2 <- viz_mmseqs + new_scale_fill()


 # Color scale for number of hits
my_integer_breaks <- seq(0, 4, 1)
my_factor_breaks <- as.factor(my_integer_breaks) 
number_of_colors <- length(my_integer_breaks) 
# Generate the color palette
color_palette_generator <- colorRampPalette(c("white", "darkred"))
my_color_palette <- color_palette_generator(number_of_colors)
# Map each factor level to a color
color_values <- setNames(my_color_palette, my_factor_breaks)


mmseqs_hits_df <- as.data.frame(mmseqs_combined[c("HtsAHits", "SirAHits", "CntAHits", "FhuD1Hits", "FhuD2Hits", "SstDHits")], row.names = mmseqs_combined$genome)
rownames(mmseqs_hits_df) <- mmseqs_combined$genome
mmseqs_hits_df[is.na(mmseqs_hits_df)] <- 0 # No hits above the set threshold, hence set to 0
for (genome in missing_genomes) {
  mmseqs_hits_df[genome, ] <- 0
}

# Change columns to factors
mmseqs_hits_df <- mmseqs_hits_df %>% mutate(across(everything(), as.character))

colnames(mmseqs_hits_df) <- c("htsA", "sirA", "cntA", "fhuD1", "fhuD2", "sstD")
viz_mmseqs <- gheatmap(temp2,mmseqs_hits_df, low="white", high="darkred", colnames_position = "top", font.size=2,offset=0.0575, width=.2,
         colnames_angle=90, colnames_offset_y = 3.75) + 
         scale_fill_manual(values = color_values, name = "Number of putative homologs", breaks = my_factor_breaks) 


# Add overview of assembly level
contig_overview <- read.csv(paste(input_path, "contigs_number.csv", sep = ""), header = TRUE, sep = ",")
# Remove last two characters from contig_overview$accession
contig_overview$accession <- gsub("\\.[0-9]$", "", contig_overview$accession)
# Rename Accession code using named list
contig_overview$gtdb_species_name <- named_list[contig_overview$accession]
contig_overview_df <- as.data.frame(contig_overview$assembly_level)
rownames(contig_overview_df) <- contig_overview$gtdb_species_name
colnames(contig_overview_df) <- c("assembly_level")
# For ureyliticus, set to Complete Genome as it was manually inspected with another accession code
contig_overview_df["Staphylococcus ureilyticus", "assembly_level"] <- "Complete Genome"

tmp_1 <- viz_mmseqs + new_scale_fill()

viz_contigs <- gheatmap(tmp_1,contig_overview_df, colnames_position = "top", font.size=2,  offset=0.075,width=.05, colnames_angle=90,colnames_offset_y =9.5 ) +
                scale_fill_manual(limits = c("Complete Genome", "Scaffold", "Contig"), values = c("#a1d99b", "#fdae6b", "#fc9272")) + guides(fill=guide_legend(title="Assembly level"))

viz_contigs
# export figure
ggsave(paste(output_path,"heatmap_overview_contigs.pdf", sep = ""), plot = viz_contigs, width = 30, height = 21, units = "cm", dpi = 600)
# The figure needs to be edited in Inkscape to adapt the order of the legends, resort the columns and improve readability of text


# Create overview table
antismash_metadata$GTDB <- rownames(antismash_metadata)
unknown_bgc$GTDB <- rownames(unknown_bgc)
# add suffix _similarity to colnames
colnames(mmseqs_combined_df) <- paste0(colnames(mmseqs_combined_df), "_similarity")
mmseqs_combined_df$GTDB <- rownames(mmseqs_combined_df)
# add suffix _hits to colnames
colnames(mmseqs_hits_df) <- paste0(colnames(mmseqs_hits_df), "_hits")
mmseqs_hits_df$GTDB <- rownames(mmseqs_hits_df)
contig_overview_df$GTDB <- rownames(contig_overview_df)

# Join all by rownames
overview_table <- Reduce(function(x, y) merge(x, y, by = "GTDB"), list(antismash_metadata, unknown_bgc, mmseqs_combined_df, mmseqs_hits_df, contig_overview_df))
# add column accession code
overview_table$accession <- reverse_named_list[overview_table$GTDB]
# change to character
overview_table$accession <- as.character(overview_table$accession)
# Sort columns
overview_table <- overview_table[c("accession", "GTDB", "assembly_level", "sfa", "sbn", "cnt", "others", "htsA_similarity", "sirA_similarity", "cntA_similarity", "fhuD1_similarity", "fhuD2_similarity", "sstD_similarity", "htsA_hits", "sirA_hits", "cntA_hits", "fhuD1_hits", "fhuD2_hits", "sstD_hits")]
colnames(overview_table)[1] <- "Accession NCBI"
colnames(overview_table)[2] <- "GTDB species"
colnames(overview_table)[3] <- "Assembly Level"

# Export overview table
write.csv(overview_table, file = paste0(output_path, "overview_table.csv", sep = ""), row.names = FALSE)