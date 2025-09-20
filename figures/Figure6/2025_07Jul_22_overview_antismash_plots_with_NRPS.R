library(ComplexUpset)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)


# --- Configuration ---
adapt_data_mmseqs <- function(path_mmseqs) {
read_mmseqs <- read.table(path_mmseqs, header=FALSE)
colnames(read_mmseqs) <- c("query","target","fident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits","qcov")
head(read_mmseqs)
# from long to wide format
read_mmseqs_adapted <- read_mmseqs %>%
    filter(query %in% c("HtsA", "SirA", "FhuD1", "FhuD2", "SstD", "CntA")) %>%
    filter(fident > 0.5) %>%
    mutate(genome_id = sapply(target, function(x) {
        parts <- strsplit(x, "_", fixed = TRUE)[[1]]
        paste(head(parts, 2), collapse = "_")
    }))
head(read_mmseqs_adapted)

mmseqs_aggregated_number <- read_mmseqs_adapted %>%
    filter(fident > 0.5) %>%
    group_by(genome_id, query) %>%
    summarise(fident = n(), .groups = 'drop')

mmseqs_aggregated_max <- read_mmseqs_adapted %>%
    filter(fident > 0.5) %>%
    group_by(genome_id, query) %>%
    summarise(fident = max(fident), .groups = 'drop')
# Create a wide format dataframe
wide_mmseqs <- mmseqs_aggregated_max %>%
    select(genome_id, query, fident) %>%
    pivot_wider(names_from = query, values_from = fident) %>% 
    # replace all non-zero values with 1
    mutate(across(-genome_id, ~ ifelse(!is.na(.), 1, 0))) %>%
    # replace NA with 0
    mutate(across(-genome_id, ~ ifelse(is.na(.), 0, .))) 
return(wide_mmseqs)
}



prepare_data <- function(file_path, path_mmseqs = NULL, output_path=NULL, path_assemblies=NULL,max_number=1100, min_intersection=2) {
# Columns that might contain vector-like percentage strings
percentage_cols_to_process <- c("staphyloferrin A", "staphyloferrin B", "staphylopine")
# Column name for the number of records
records_col_name <- "number_of_records"
# Name to assign to the first column (which has an empty header in the file)
first_col_new_name <- "GCF_ID"


# --- Read and Prepare Data ---

# Read the CSV.
# sep=";" for semicolon separated values.
# header=TRUE because the first line is the header.
# stringsAsFactors=FALSE is good practice.
# check.names=FALSE preserves original column names, including the empty first one.
# na.strings handles empty strings as NA, if any, apart from vector notation.
tryCatch({
  df <- read.csv(file_path, sep = ";", header = TRUE, stringsAsFactors = FALSE, 
                 check.names = FALSE, na.strings = c("NA", ""))
}, error = function(e) {
  stop(paste("Error reading CSV file:", e$message, 
             "\nPlease ensure the file path is correct and the file is accessible."))
})


# Rename the first column (which was loaded with an empty name "")
if (ncol(df) > 0 && colnames(df)[1] == "") {
  colnames(df)[1] <- first_col_new_name
} else {
  warning("First column name was not empty as expected, or no columns found. Skipping rename of first column.")
}

# Ensure the 'number_of_records' column is numeric
if (!records_col_name %in% colnames(df)) {
  stop(paste("Critical column '", records_col_name, "' not found in the CSV.", sep=""))
}
df[[records_col_name]] <- as.integer(df[[records_col_name]])

# --- Process Rows and Columns ---

for (i in 1:nrow(df)) {
  for (col_name in percentage_cols_to_process) {
    # Check if the column exists in the dataframe
    if (!col_name %in% colnames(df)) {
      warning(paste("Column '", col_name, "' not found. Skipping processing for this column.", sep=""))
      next # Move to the next column in percentage_cols_to_process
    }

    current_cell_value <- df[i, col_name]

    # Check if the cell value is a character string and contains '['
    if (grepl("[", current_cell_value, fixed = TRUE)) {
      # Remove brackets and any spaces
      cleaned_value_str <- gsub("\\[|\\]|\\s+", "", current_cell_value)
      
      # Split by comma and convert to numeric
      # Ensure that empty strings (from "[]" or "[ ]") result in an empty numeric vector
      if (cleaned_value_str == "") {
        numbers_vec <- numeric(0)
      } else {
        numbers_str_vec <- strsplit(cleaned_value_str, ",", fixed = TRUE)[[1]]
        numbers_vec <- as.numeric(numbers_str_vec)
        numbers_vec <- numbers_vec[!is.na(numbers_vec)] # Remove any NAs from parsing non-numbers
      }

      if (length(numbers_vec) > 0) {
        # Store the maximum value
        df[i, col_name] <- max(numbers_vec, na.rm = TRUE)
        
        # Adjust 'number_of_records'
        # The adjustment is (length_of_vector - 1)
        # If length is 1, adjustment is 0. If length is 0, this logic is skipped.
        adjustment <- length(numbers_vec) - 1
        df[i, records_col_name] <- df[i, records_col_name] - adjustment
        
      } else {
        # If vector was empty or unparseable (e.g., "[]", "[abc]")
        df[i, col_name] <- NA 
      }
    }
  }
}
# --- Final Data Type Conversion ---
# Convert the processed percentage columns to numeric
for (col_name in percentage_cols_to_process) {
  if (col_name %in% colnames(df)) {
    df[[col_name]] <- as.numeric(as.character(df[[col_name]]))
  }
}

# --- Create 'Other' Column ---
df$Other <- ifelse(as.numeric(df$number_of_records) == as.numeric(as.logical(df["staphyloferrin A"] > 0)) + as.numeric(as.logical(df["staphyloferrin B"] > 0)) + as.numeric(as.logical(df$staphylopine > 0)), "FALSE", "TRUE")



upset_vector <- apply(df, 1, function(row) {
    return_vector <- c()
    for (col in c("staphyloferrin A", "staphyloferrin B", "staphylopine")) {
        row_col_value <- as.numeric(row[col])
        if (row_col_value > 0 && row_col_value < 100) {
            return_vector <- c(return_vector, paste0(col, "_incomplete"))
        } else if (row_col_value == 100) {
            return_vector <- c(return_vector, paste0(col, "_complete"))
        } 
    }
    if (row["Other"] == "TRUE") {
        return_vector <- c(return_vector, "Other BGC")
    }
    if (length(return_vector) == 0) {
        return_vector <- c("None")
    }
    return(return_vector)

})
names(upset_vector) <- df$GCF_ID

vector_as_df <- from_list_to_df(upset_vector, path_assemblies)
data_mmseqs <- adapt_data_mmseqs(path_mmseqs)
# merge the two dataframes
vector_as_df <- merge(vector_as_df, data_mmseqs, by.x = "row.names", by.y = "genome_id", all.x = TRUE)
vector_as_df <- vector_as_df %>%
    mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) 
    # drop the column none
if ("None" %in% colnames(vector_as_df)) {
   vector_as_df <- vector_as_df %>%
    select(-None)
}
# save the dataframe to a csv file
write.csv(vector_as_df, 
          file = output_path, 
          row.names = FALSE)

# sum of vector_as_df
 sum_cols <- colSums(vector_as_df[-1])
 # get those with a minimal value of 1
sum_cols_names <- names(sum_cols[sum_cols > max(c(min_intersection-1, 0))])
print(sum_cols_names)
# Create queries for the upset plot
# For column names with complete fill with a specific color
# For column names with incomplete fill with a different color
list_queries <- lapply(colnames(vector_as_df)[-1], function(x) {
  if (x %in% sum_cols_names) {
    # If the column is in the sum_cols_names, create a query
    # Check if the column name contains "_incomplete" or "_complete"
     if (grepl("_incomplete$", x)) {
    upset_query(set = x, fill = "#9FD0DF")  # for incomplete
  } else if (grepl("_complete$", x)) {
    upset_query(set = x, fill = "#000078")  #  for complete
  }
  else {
     NULL
  }
  } else {
    NULL  # Skip columns that do not meet the criteria
  }
 

})
list_queries_final <- list()
# Remove NULL queries
for (i in seq_along(list_queries)) {
  if (!is.null(list_queries[[i]])) {
    list_queries_final[[i]] <- list_queries[[i]]
  }
}
print(list_queries_final)
print(head(vector_as_df))
# Remove columns that are not in sum_cols_names
vector_as_df <- vector_as_df[, c("Row.names", sum_cols_names)]

# rename columns using a mapper, staphyloferrin A_complete -> sfa_complete
colnames(vector_as_df) <- gsub("staphyloferrin A_complete", "sfa_complete", colnames(vector_as_df))
colnames(vector_as_df) <- gsub("staphyloferrin A_incomplete", "sfa_incomplete", colnames(vector_as_df))
colnames(vector_as_df) <- gsub("staphylopine_complete", "cnt_complete", colnames(vector_as_df))
colnames(vector_as_df) <- gsub("staphylopine_incomplete", "cnt_incomplete", colnames(vector_as_df))
colnames(vector_as_df) <- gsub("staphyloferrin B_incomplete", "sbn_incomplete", colnames(vector_as_df))
colnames(vector_as_df) <- gsub("staphyloferrin B_complete", "sbn_complete", colnames(vector_as_df))
colnames(vector_as_df) <- gsub("HtsA", "htsA", colnames(vector_as_df))
colnames(vector_as_df) <- gsub("SirA", "sirA", colnames(vector_as_df))
colnames(vector_as_df) <- gsub("FhuD1", "fhuD1", colnames(vector_as_df))
colnames(vector_as_df) <- gsub("FhuD2", "fhuD2", colnames(vector_as_df))
colnames(vector_as_df) <- gsub("SstD", "sstD", colnames(vector_as_df))
colnames(vector_as_df) <- gsub("CntA", "cntA", colnames(vector_as_df))



return(
  vector_as_df           
)

}
from_list_to_df <- function(vector,path_assemblies){
        # 1. Get all unique values that will become the column headers
    all_values <- unique(unlist(vector))

    # 2. Use lapply to create a list of binary vectors
    binary_list <- lapply(vector, function(element) {
    # For each element, check which of the all_values are present
    # as.integer() converts TRUE/FALSE to 1/0
    as.integer(all_values %in% element)
    })

    # 3. Combine the list of vectors into a matrix
    binary_matrix <- do.call(rbind, binary_list)

    read_assemblies <- read.csv(path_assemblies, stringsAsFactors = FALSE, col.names = c("accession"))


    # 4. Assign the column names and convert to a data frame for nice printing
    colnames(binary_matrix) <- all_values


    
    binary_df_base <- as.data.frame(binary_matrix)
        # filter df to only include accessions that are in read_assemblies
    binary_df_base_filt <- binary_df_base[rownames(binary_df_base) %in% read_assemblies$accession, ]
    print(colnames(binary_df_base_filt))
    return(binary_df_base_filt)
}

input_path <- "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/results/repo/StaphyDiversity/figures/Figure6/input/"
output_all <- "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/results/repo/StaphyDiversity/figures/Figure6/output/"


capitis_path <- paste(input_path, "capitis/summary_with_clusterblast_with_percentage_clustercompare.csv", sep = "")
capitis_mmseqs_path <- paste(input_path, "capitis/alnResult_new.m8", sep = "")
output_path_capitis <- paste(output_all, "upset_capitis_ccompare.csv", sep = "")
path_assemblies <- paste(input_path, "capitis/accessions_filtered_2.csv", sep = "")
upset_capitis <- prepare_data(capitis_path, capitis_mmseqs_path, output_path_capitis, path_assemblies, 180,0)
png(paste(output_all, "upset_capitis_new_ccompare.png", sep = ""), width = 600, height = 600)
print(upset(
  upset_capitis,
  intersect = colnames(upset_capitis)[-1],
  min_size = 0,
  queries = list(
      upset_query(set = "sfa_complete", fill = "#000078"),
      upset_query(set = "cnt_complete", fill = "#000078"),
      upset_query(set = "cnt_incomplete", fill = "#9FD0DF")
      ),
      encode_sets=FALSE,  # for annotate() to select the set by name disable encoding
      set_sizes=(
          upset_set_size()
          + geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
          + expand_limits(y=180)
          + theme(axis.text.x=element_text(angle=90))
      )
))
dev.off()
svg(paste(output_all, "upset_capitis_new_ccompare.svg", sep = ""), width = 6, height = 6)
print(upset(
  upset_capitis,
  intersect = colnames(upset_capitis)[-1],
  min_size = 0,
  queries = list(
      upset_query(set = "sfa_complete", fill = "#000078"),
      upset_query(set = "cnt_complete", fill = "#000078"),
      upset_query(set = "cnt_incomplete", fill = "#9FD0DF")
      ),
      encode_sets=FALSE,  # for annotate() to select the set by name disable encoding
      set_sizes=(
          upset_set_size()
          + geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
          + expand_limits(y=180)
          + theme(axis.text.x=element_text(angle=90))
      )
))
dev.off()

lugdunensis_path <- paste(input_path, "lugdunensis/summary_with_clusterblast_with_percentage_clustercompare.csv", sep = "")
lugdunensis_mmseqs_path <- paste(input_path, "lugdunensis/alnResult_new.m8", sep = "")
output_path_lugdunensis <- paste(output_all, "upset_lugdunensis_ccompare.csv", sep = "")
path_assemblies_lugdunensis <- paste(input_path, "lugdunensis/accessions_filtered_2.csv", sep = "")
upset_lugdunensis <- prepare_data(lugdunensis_path, lugdunensis_mmseqs_path, output_path_lugdunensis, path_assemblies_lugdunensis, 150,0)
nrps_putatives <- read.csv(paste(input_path, "lugdunensis/NRPS_overview_files.txt", sep = ""), header=FALSE)
upset_lugdunensis$NRPS <- ifelse(upset_lugdunensis$Row.names %in% nrps_putatives$V1, 1, 0)

png(paste(output_all, "upset_lugdunensis_new_ccompare.png", sep = ""), width = 600, height = 600)
print(
upset(
  upset_lugdunensis,
  intersect = colnames(upset_lugdunensis)[-1],
  min_size = 0,
  queries = list(
      upset_query(set = "sfa_incomplete", fill = "#9FD0DF")
      ),
      encode_sets=FALSE,  # for annotate() to select the set by name disable encoding
      set_sizes=(
          upset_set_size()
          + geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
          + expand_limits(y=150)
          + theme(axis.text.x=element_text(angle=90))
      )
))
dev.off()

svg(paste(output_all, "upset_lugdunensis_new_ccompare.svg", sep = ""), width = 6, height = 6)
print(
upset(
  upset_lugdunensis,
  intersect = colnames(upset_lugdunensis)[-1],
  min_size = 0,
  queries = list(
      upset_query(set = "sfa_incomplete", fill = "#9FD0DF")
      ),
      encode_sets=FALSE,  # for annotate() to select the set by name disable encoding
      set_sizes=(
          upset_set_size()
          + geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
          + expand_limits(y=150)
          + theme(axis.text.x=element_text(angle=90))
      )
))
dev.off()

hominis_path <- paste(input_path, "hominis/summary_with_clusterblast_with_percentage_clustercompare.csv", sep = "")
hominis_mmseqs_path <- paste(input_path, "hominis/alnResult_new.m8", sep = "")
output_path_hominis <- paste(output_all, "upset_hominis_ccompare.csv", sep = "")
path_assemblies_hominis <- paste(input_path, "hominis/accessions_filtered_2.csv", sep = "")
upset_hominis <- prepare_data(hominis_path, hominis_mmseqs_path, output_path_hominis, path_assemblies_hominis, 300, 0)
png(paste(output_all, "upset_hominis_new_ccompare.png", sep = ""), width = 600, height = 600)
print(
upset(
  upset_hominis,
  intersect = colnames(upset_hominis)[-1],
  min_size = 0,
  queries = list(
      upset_query(set = "sfa_incomplete", fill = "#9FD0DF")
      ),
      encode_sets=FALSE,  # for annotate() to select the set by name disable encoding
      set_sizes=(
          upset_set_size()
          + geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
          + expand_limits(y=300)
          + theme(axis.text.x=element_text(angle=90))
      )
))
dev.off()

svg(paste(output_all, "upset_hominis_new_ccompare.svg", sep = ""), width = 6, height = 6)
plot(
  upset(
  upset_hominis,
  intersect = colnames(upset_hominis)[-1],
  min_size = 0,
  queries = list(
      upset_query(set = "sfa_incomplete", fill = "#9FD0DF")
      ),
      encode_sets=FALSE,  # for annotate() to select the set by name disable encoding
      set_sizes=(
          upset_set_size()
          + geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
          + expand_limits(y=300)
          + theme(axis.text.x=element_text(angle=90))
      )  
    )
)
dev.off()

epidermidis_path <- paste(input_path, "epidermidis/summary_with_clusterblast_with_percentage_clustercompare.csv", sep = "")
epidermidis_mmseqs_path <- paste(input_path, "epidermidis/alnResult_new.m8", sep = "")
output_path_epidermidis <- paste(output_all, "upset_epidermidis_ccompare.csv", sep = "")
path_assemblies_epidermidis <- paste(input_path, "epidermidis/accessions_filtered_2.csv", sep = "")
upset_epidermidis <- prepare_data(epidermidis_path, epidermidis_mmseqs_path, output_path_epidermidis, path_assemblies_epidermidis, 1500,2)

png(paste(output_all, "upset_epidermidis_new_ccompare.png", sep = ""), width = 1200, height = 800)
print(upset(
  upset_epidermidis,
  intersect = colnames(upset_epidermidis)[-1],
  min_size = 0,
  queries = list(
      upset_query(set = "sfa_complete", fill = "#000078"),
      upset_query(set = "sfa_incomplete", fill = "#9FD0DF"),
      upset_query(set = "cnt_complete", fill = "#000078"),
      upset_query(set = "cnt_incomplete", fill = "#9FD0DF")
      ),
      encode_sets=FALSE,  # for annotate() to select the set by name disable encoding
      set_sizes=(
          upset_set_size()
          + geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
          + expand_limits(y=1500)
          + theme(axis.text.x=element_text(angle=90))
      )
))
dev.off()

svg(paste(output_all, "upset_epidermidis_new_ccompare.svg", sep = ""), width = 12, height = 8)
print(upset(
  upset_epidermidis,
  intersect = colnames(upset_epidermidis)[-1],
  min_size = 0,
  queries = list(
      upset_query(set = "sfa_complete", fill = "#000078"),
      upset_query(set = "sfa_incomplete", fill = "#9FD0DF"),
      upset_query(set = "cnt_complete", fill = "#000078"),
      upset_query(set = "cnt_incomplete", fill = "#9FD0DF")
      ),
      encode_sets=FALSE,  # for annotate() to select the set by name disable encoding
      set_sizes=(
          upset_set_size()
          + geom_text(aes(label=after_stat(count)), hjust=1.1, stat='count')
          + expand_limits(y=1500)
          + theme(axis.text.x=element_text(angle=90))
      )
))
dev.off()

