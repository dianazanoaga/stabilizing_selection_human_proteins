# load necessary libraries

library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(rbgen)
library(genio)
library(Matrix)
library(mbend)
library(gridExtra)
library(readxl)
library(ggplot2)

# Set the directory containing the CSV files
folder_path <- "/data/FAC/FBM/DBC/zkutalik/default_sensitive/old/UKBB/pheno"

# Get a list of all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Load file with refresh information about the files
corresp <- as.data.frame(read_excel("/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/correspondence_uk_pheno.xlsx"))
    
# Specify the column names you're interested in
column_of_interest <- c("2405", "2734") # respectively : number of children fathered, number of live births female
													   

# Initialize a data frame to store the results

new_pheno <- data.frame(id = character(),
                        file = character(),
                        col = character(),
                        stringsAsFactors = FALSE)


# Loop over each column of interest
for (i in 1:length(column_of_interest)) {
  
  print(i)

  # Loop over each CSV file
  for (file in csv_files) {
    
    # Read only the header of the CSV file
    data <- fread(file, nrows = 0, header = TRUE)
    
    # Get the column names
    column_names <- colnames(data)
    
    matching_indices <- grep(paste0("^", column_of_interest[i], "-"), column_names)
    
    if (length(matching_indices) > 0) {
      
	  # Save the column of interest, file, and matching indices to the new_pheno data frame
	  mom_pheno <- data.frame(id = column_of_interest[i],
                              file = file,
                              col = paste(column_names[matching_indices], collapse = "_"),
                              stringsAsFactors = FALSE)
	
      new_pheno <- rbind(new_pheno, mom_pheno)

    }

  }
}

# Save the numeric part of the files and check which one corresponds to refresh 2
new_pheno$numeric_part <- as.numeric(gsub("\\D", "", new_pheno$file))

# Filter rows where the numeric part is in the `refresh 2` column of corresp
filtered_data <- new_pheno[new_pheno$numeric_part %in% corresp$`refresh 2`, ]


#################### Analysis of "number of children - father" #2405  #####################################################

file_paths <- filtered_data[filtered_data$id == "2405", "file"] #take file path where I have found #2405
col <- filtered_data[filtered_data$id == "2405" & filtered_data$file == file_paths[1], "col"] #take all columns name I have found for the phenotype #2405

data_n_child_fathered <- as.data.frame(fread(file_paths[1], header = T, select = c("eid", unlist(strsplit(col, "_"))))) #access the file, select only eid and the columns related to #2405

# Put NA where values are less than 0 (-1 and -3, corresponding to "Do not know" and "Prefer not to answer") and take average across instances 

data_n_child_fathered[data_n_child_fathered < 0 ] <- NA
data_n_child_fathered$n.children.fathered_avg <- rowMeans(data_n_child_fathered[,c(2:ncol(data_n_child_fathered))], na.rm=T)

### this is to check individuals that have different values among instances ##############################################################
check <- data_n_child_fathered[
  (data_n_child_fathered$n.children.fathered_avg != data_n_child_fathered$`2405-0.0`) &
  (data_n_child_fathered$n.children.fathered_avg != data_n_child_fathered$`2405-1.0`) &
  (data_n_child_fathered$n.children.fathered_avg != data_n_child_fathered$`2405-2.0`) &
  (data_n_child_fathered$n.children.fathered_avg != data_n_child_fathered$`2405-3.0`),
]

check <- check %>%
  filter(rowSums(is.na(.)) != ncol(.))

eid_string <- paste(check$eid, collapse = ", ")
print("'eid' with different values among instances:")
print(eid_string)
print(paste0("In total: ",length(check$eid)))

#################################################################################################################################################


dd <- data_n_child_fathered[data_n_child_fathered$"n.children.fathered_avg" == max(data_n_child_fathered$"n.children.fathered_avg", na.rm =T),]
dd <- dd %>%
  filter(rowSums(is.na(.)) != ncol(.))

eid_string <- paste(dd$eid, collapse = ", ")
print("'eid' with max children 100:")
print(eid_string)

p1 <- ggplot(data_n_child_fathered, aes(x = n.children.fathered_avg)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  geom_vline(aes(xintercept = mean(data_n_child_fathered$"n.children.fathered_avg", na.rm = TRUE)), color = "red", linetype = "dashed", size = 1) +
  labs(title = paste0("Histogram of Value- (Mean : ", mean(data_n_child_fathered$"n.children.fathered_avg", na.rm = TRUE), " )"), x = "Value", y = "Frequency") +
  theme_minimal()

p2 <- ggplot(data_n_child_fathered, aes(y = n.children.fathered_avg)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  labs(title = "Boxplot of Value", y = "Value") +
  theme_minimal()

# Arrange the plots together
grid.arrange(p1, p2, nrow = 2)



#################### Analysis of "number of live births - mother" #2734   ##################################################


file_paths <- filtered_data[filtered_data$id == "2734", "file"]
col <- filtered_data[filtered_data$id == "2734" & filtered_data$"file" == file_paths[1], "col"]

data_n_live_births_fem <- as.data.frame(fread(file_paths[1], header = T, select = c("eid", unlist(strsplit(col, "_")))))

data_n_live_births_fem[data_n_live_births_fem == -3] <- NA
data_n_live_births_fem$n.live.births_avg <- rowMeans(data_n_live_births_fem[,c(2:ncol(data_n_live_births_fem))], na.rm=T)


check <- data_n_live_births_fem[
  (data_n_live_births_fem$n.live.births_avg != data_n_live_births_fem$`2734-0.0`) &
  (data_n_live_births_fem$n.live.births_avg != data_n_live_births_fem$`2734-1.0`) &
  (data_n_live_births_fem$n.live.births_avg != data_n_live_births_fem$`2734-2.0`) &
  (data_n_live_births_fem$n.live.births_avg != data_n_live_births_fem$`2734-3.0`),
]


check <- check %>%
  filter(rowSums(is.na(.)) != ncol(.))

eid_string <- paste(check$eid, collapse = ", ")
print("'eid' with different values among instances:")
print(eid_string)
print(paste0("In total: ",length(check$eid)))


dd <- data_n_live_births_fem[data_n_live_births_fem$"n.live.births_avg" == max(data_n_live_births_fem$"n.live.births_avg", na.rm =T),]
dd <- dd %>%
  filter(rowSums(is.na(.)) != ncol(.))

eid_string <- paste(dd$eid, collapse = ", ")
print("'eid' with max children:")
print(eid_string)


p1 <- ggplot(data_n_live_births_fem, aes(x = n.live.births_avg)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  geom_vline(aes(xintercept = mean(data_n_live_births_fem$"n.live.births_avg", na.rm = TRUE)), color = "red", linetype = "dashed", size = 1) +
  labs(title = paste0("Histogram of Value- (Mean : ", mean(data_n_live_births_fem$"n.live.births_avg", na.rm = TRUE), " )"), x = "Value", y = "Frequency") +
  theme_minimal()

p2 <- ggplot(data_n_live_births_fem, aes(y = n.live.births_avg)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  labs(title = "Boxplot of Value", y = "Value") +
  theme_minimal()

# Arrange the plots together
grid.arrange(p1, p2, nrow = 2)


### Save data ###################################################################################################################

save(data_n_child_fathered, data_n_live_births_fem, file = "/data/FAC/FBM/DBC/zkutalik/default_sensitive/mzanoaga/selection_proteins/00_file_preparation/data/pheno_n.children_M_and_F.RData")

print("END!")


