# Allele  association  script 1.0.0 https://github.com/4C554D43
# Desc: Generates quantitative report on alleles associated with a phenotype (or the absence thereof)
# 
# EXPERIMENTAL BUILD. USE AT YOUR OWN RISK
# Version 1.0.0,  25-11-2023: Initial version
#


# ----- CONFIG -----
#directory for input files. 
input_dir = r"(path)" 
filename_phenotype = "pheno.csv"
filename_genotype =  "geno.csv"

output_dir = r"(path)"
output_name = "association_report" 
create_logfile = TRUE
create_dir = TRUE #create new folder for each run?
# ----- END CONFIG -----

#include and init
if (!requireNamespace("devtools", quietly = TRUE)) {install.packages("devtools")}
if (!requireNamespace("readxl", quietly = TRUE)) {install.packages("readxl")}
if (!requireNamespace("dplyr", quietly = TRUE)) {install.packages("dplyr")}
library(readxl)
library(dplyr)

# ----- functions -----
timedate_string <- function()
{
  current_datetime <- Sys.time() 
  formatted_datetime <- format(current_datetime, "%Y-%m-%d_%H-%M-%S")
  return (as.character(formatted_datetime))
}#end timedate
print_timestamped <- function(first, second = NULL, third = NULL, fourth = NULL, fifth = NULL, sixth = NULL, se = NULL, ei = NULL, ne = NULL)
{
  message <- paste(first, second, third, fourth, fifth, sixth,se, ei, ne, collapse = " ")
  if ("\n" %in% message) {
    cat(timedate_string(), ": ", message, "\n")
  } else {
    cat(timedate_string(), ": ", message, "\n", sep = "")
  }
}#end print
trim_columns <- function(mat) {
  keep_columns <- c()
  for (col_index in 1:ncol(mat)) 
  {
    if (length(unique(na.omit(mat[, col_index]))) > 1)    # Check if more than one unique non-missing value in column
      keep_columns <- c(keep_columns, col_index)
  }
  removed_columns <- ncol(mat) - length(keep_columns) # Determine the number of columns removed
  result_matrix <- mat[, keep_columns, drop = FALSE]  # Subset the matrix to keep only the selected columns
  #print and return
  print_timestamped("trim_column(): Columns removed:", removed_columns)
  print_timestamped("trim_column(): Columns kept:", length(keep_columns))
  return(result_matrix)
}#end trim

binary_encode_matrix <- function(mat) { #assign most common value 1, others with 0. Missing assigned 0
  for (col_index in 1:ncol(mat)) 
  {
    col_values <- mat[, col_index]
    value_counts <- table(na.omit(col_values))
    most_common_value <- as.numeric(names(value_counts[which.max(value_counts)]))
    mat[, col_index] <- ifelse(is.na(col_values), 1, ifelse(col_values == most_common_value, 1, 0))
  }
  return(mat)
} #end binary encode
print_phenostats <- function(mat)#prints phenotype input matrix stats
{
  num_true = 0
  num_false = 0
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      # Print the value and its position
      if (mat[i, j] == 1){
        num_true = num_true + 1
      }else{num_false = num_false + 1}
      print_timestamped("Value at position (", i, ",", j, "): ", mat[i, j])
    }
  }
  print_timestamped("Phenotype samples: ", nrow(severity_matrix), "; severe: ", num_true, "; mild: ", num_false)
}#end phenostats

# ----- main -----
#resolve I/O
if (create_dir == TRUE)
{
  output_dir = paste(output_dir, timedate_string(), "_Report")
  dir.create(output_dir)
  if (!dir.exists(output_dir))
  {
    print_timestamped("[ERR] Could not create output directory")
    stop()
  }
}
sink(file.path(output_dir, r"(log.txt)"), split = TRUE, append = TRUE)
Sys.Date()
Sys.time()
R.Version()
path_phenotype = file.path(input_dir, filename_phenotype) #sample severity matrix path) 
path_genotype = file.path(input_dir, filename_genotype) #sample data matrix path
if (!file.exists((path_phenotype))|| !file.exists((path_genotype)))
{
  print_timestamped("[ERR] An input file (severity, sample data or tree) was not found!")
  stop()
}






phenotype_matrix <- as.matrix(read.csv(path_phenotype, header = TRUE, row.names = 1, sep=';'))
genotype_matrix <- as.matrix(read.csv(path_genotype,header = TRUE, row.names = 1, sep=';'))

genotype_matrix_trimmed = trim_columns(genotype_matrix) #trim columns with same alleles

find_associated_genes <- function(target_phenotype, phenotype_matrix, genotype_matrix) {
  hits = 0
  misses = 0
  
  
  for (col in colnames(genotype_matrix)) {#traverse columns
    geno_col <- genotype_matrix[, col] #isolate column
    
    if (all(is.na(geno_col))) { #skip empty
      misses <- misses + 1
      next}
    
    checked_variations <- c() #hold already checked gene variations
    
    # Track gene variations and their counts for each phenotype
    variation_counts_1 <- table(geno_col[phenotype_matrix == 1 & !is.na(geno_col)])
    variation_counts_0 <- table(geno_col[phenotype_matrix == 0 & !is.na(geno_col)])
    if (target_phenotype == 0)
    {
      #------------------------------
      # Check variations associated with phenotype 0
      for (variation in names(variation_counts_0)) 
      {
        if (!is.na(variation_counts_0[variation]) && variation_counts_0[variation] > 1 && is.na(variation_counts_1[variation]) && !(variation %in% checked_variations)) 
        {
          print_timestamped("Gene set", col, "associated with phenotype 0 for variation", variation, "appeared", variation_counts_0[variation], "times")
          hits = hits + 1
          checked_variations <- c(checked_variations, variation)# Track the checked variation for this column
        }
      }
      if (all(is.na(phenotype_matrix[geno_col == 0])))  # Check if all values for phenotype 1 are missing in the current column
        misses <- misses + 1
      #---------------------------------------------------------------------
    }else if (target_phenotype == 1)
    {
      #------------------------------
      # Check variations associated with phenotype 1
      for (variation in names(variation_counts_1)) 
      {
        if (!is.na(variation_counts_1[variation]) && variation_counts_1[variation] > 1 && is.na(variation_counts_0[variation]) && !(variation %in% checked_variations)) 
        {
          print_timestamped("Gene set", col, "associated with phenotype 1 for variation", variation, "appeared", variation_counts_1[variation], "times")
          hits = hits + 1
          checked_variations <- c(checked_variations, variation)# Track the checked variation for this column
        }
      }
      if (all(is.na(phenotype_matrix[geno_col == 1])))  # Check if all values for phenotype 1 are missing in the current column
        misses <- misses + 1
      #---------------------------------------------------------------------
      
    }else
    {
      print_timestamped("Only supports binary phenotype")
      
    }
    

  }
  
  # Print the results
  print_timestamped("Checked pheno:", target_phenotype)
  print_timestamped("Number of columns hit:", hits)
  print_timestamped("Number of columns not hit:", misses)
}

find_associated_genes(0, phenotype_matrix, genotype_matrix_trimmed)
find_associated_genes(1, phenotype_matrix, genotype_matrix_trimmed)


sink()



