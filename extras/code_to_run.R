################################################################################
#
# Installing the package
#
################################################################################
devtools::install_github("Ckonsa/ClusterPatientsTrajectories")
library(ClusterPatientsTrajectories)



################################################################################
#
# Cluster patients in CLI
#
################################################################################
only_cluster_patients <- function() {
  # VARIABLES - CHANGE ACCORDING TO YOUR NEEDS
  file_path <- "./tmp/datasets/example_data.csv"
  records_range <- c(3, 13)
  trajectory_length <- 10
  clusters_amount <- 4
  punishments <- c(14, 11, 7, 4, 2)
  use_dtw <- FALSE

  data <- read.csv(file_path)
  patients <- get_filtered_patients(data, records_range)
  filtered_data <- data %>% filter(SUBJECT_ID %in% patients)

  patients_data <- collect_patient_data(filtered_data, patients, trajectory_length)
  compare_matrix <- compare_trajectories(patients_data, punishments, use_dtw)
  clustering_result <- cutree(flashClust(dist(compare_matrix)), clusters_amount)

  clusters <- clusters_with_ids(patients_data$IDs, clustering_result)
  return(clusters)
}
clusters <- only_cluster_patients()
print(clusters)



################################################################################
#
# Run the workflow in CLI
#
################################################################################
run_workflow <- function() {
  # VARIABLES - CHANGE ACCORDING TO YOUR NEEDS
  file_path <- "./tmp/datasets/example_data.csv"
  records_range <- c(3, 13)
  trajectory_length <- 10
  clusters_amount <- 4
  punishments <- c(14, 11, 7, 4, 2)
  use_dtw <- FALSE
  most_frequent_drug_count <- 5
  normalize <- FALSE


  data <- read.csv(file_path)
  patients <- get_filtered_patients(data, records_range)
  filtered_data <- data %>% filter(SUBJECT_ID %in% patients)

  patients_data <- collect_patient_data(filtered_data, patients, trajectory_length)
  compare_matrix <- compare_trajectories(patients_data, punishments, use_dtw)
  row_dendrogram <- flashClust(dist(compare_matrix))
  clustering_result <- pheatmap(compare_matrix, cutree_rows = clusters_amount, cluster_rows = row_dendrogram$order, cluster_cols = FALSE)

  clusters <- clusters_with_ids(patients_data$IDs, cutree(clustering_result$tree_row, clusters_amount))
  most_frequent_drugs <- get_most_frequent_values(data, "DRUG_SOURCE_VALUE", most_frequent_drug_count)
  plot_clusters_drugs_count(filtered_data, clusters, most_frequent_drugs, normalize)
  plot_elbow(compare_matrix)
  additional_clusters_data <- get_additional_clusters_data(filtered_data, clusters, patients_data)
  print(additional_clusters_data)
}
run_workflow()



################################################################################
#
# Run GUI
#
################################################################################
shiny::runApp("./shiny/app.R")
