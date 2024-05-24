library(Rcpp)
library(dplyr)
library(pheatmap)
library(gridExtra)
library(flashClust)
sourceCpp("./src/functions.cpp")

#' Get k most frequent values in a column.
#'
#' @param data: dataframe
#' @param column: string
#' @param k: integer, k > 0
#' @return vector containing k most frequent values in a column
get_most_frequent_values <- function(data, column, k) {
  most_frequent <- data %>%
                    plyr::count(column) %>%
                    arrange(desc(freq)) %>%
                    top_n(k) %>%
                    select(all_of(column))
  return(as.vector(unlist((most_frequent))))
}

#' Remove consecutive elements of the same value.
#'
#' @param elements: vector
#' @return vector without multiple same elements in a row
#' @examples
#' remove_consecutive_elements(c(1, 2, 1, 1, 3)) => c(1, 2, 1, 3)
remove_consecutive_elements <- function(elements) {
  indexes_to_remove <- c()
  for (i in seq(along = elements)[-1]) {
    # Remove an element, if the element is the same as the element before it.
    if (elements[i - 1] == elements[i]) {
      indexes_to_remove <- append(indexes_to_remove, i)
    }
  }
  if (length(indexes_to_remove) != 0) {
    elements <- elements[-indexes_to_remove]
  }
  return(elements)
}

#' Modify the vector to length of n. Remove last elements or add filler elements.
#'
#' @param vector: vector
#' @param n: integer, n > 0
#' @return vector with length of n
#' @examples
#' modify_vector_to_size_n(c(1, 2, 3), 4) => c(1, 2, 3, -1)
#' modify_vector_to_size_n(c(1, 2, 3), 2) => c(1, 2)
modify_vector_to_size_n <- function(vector, n) {
  if (length(vector) >= n) {
    return(vector[1:n])
  } else {
    return(append(vector, rep(c(-1), each = n - length(vector))))
  }
}

#' Get patients drug trajectory of length n.
#'
#' @param data: dataframe
#' @param id: integer
#' @param n: integer n > 0
#' @return vector that is patients drug trajectory that contains drug codes.
get_patient_drug_trajectory <- function(data, id, n) {
  filtered_data <- data %>%
    filter(SUBJECT_ID == id) %>%
    arrange(DRUG_EXPOSURE_START_DATE)
  drugs <- filtered_data$DRUG_SOURCE_VALUE
  trajectory <- modify_vector_to_size_n(remove_consecutive_elements(drugs), n)
  return(trajectory)
}

#' Split string according to indexes.
#'
#' @param string: string
#' @param indexes: numeric vector
#' @return vector that contains pieces of the string
#' @examples
#' split_into_pieces("abcdef", c(1, 2, 5)) => c("a", "bcd", "ef")
split <- function(string, indexes) {
  pieces <- c()
  for (i in seq(along = indexes)) {
    start <- indexes[i]
    if (length(indexes) > i) {
      end <- indexes[i + 1] - 1
    } else {
      end <- nchar(string)
    }
    pieces <- append(pieces, substring(string, start, end))
  }
  return(pieces)
}

#' Calculate distance between two ATC codes.
#'
#' The distance is calculated by given punishments.
#' For the same elements the distance is 0.
#' If either one of the codes is missing, the distance is 1.
#'
#' @param code1: string
#' @param code2: string
#' @param punishments: numeric vector
#' @return numeric distance between two codes
#' @examples
#' distance_between_codes("R03BA02", "R04BK10", c(5,4,3,2,1)) => 7
#' distance_between_codes("R03BA02", "R03BA02", c(5,4,3,2,1)) => 0
#' distance_between_codes("-1", "R03BA02", c(5,4,3,2,1)) => 1
distance_between_codes <- function(code1, code2, punishments) {
  if (code1 == code2) return(0)
  if (code1 == -1 || code2 == -1) return(1)

  # ATC codes are split into five
  split_indexes <- c(1, 2, 4, 5, 6)
  code1_pieces <- split(code1, split_indexes)
  code2_pieces <- split(code2, split_indexes)

  distance <- 0
  for (i in seq(along = code1_pieces)) {
    # Add to distance if pieces of codes are different from each other.
    if (code1_pieces[i] != code2_pieces[i]) {
      distance <- distance + punishments[i]
    }
  }
  return(distance)
}

#' Get IDs of patients whose records count is in range.
#'
#' @param data: dataframe
#' @param records_range: numeric vector with two values - min and max record count
#' @return vector that contains IDs of patients that fit into given range
get_filtered_patients <- function(data, records_range) {
  filtered_patients <- data %>%
    group_by(SUBJECT_ID) %>%
    dplyr::summarise(records_count = n()) %>%
    filter(records_count >= records_range[1]) %>%
    filter(records_count <= records_range[2])
  return(filtered_patients$SUBJECT_ID)
}

#' Get clusters with patients IDs.
#'
#' Two vectors given as arguments have to have the same length.
#' Every element in the second argument shows into which cluster a patient in
#' the first argument on the same index belongs.
#'
#' @param patients_ids: numeric vector
#' @param clusters_with_indexes: numeric vector
#' @return vector which has a list for every cluster. Lists contain patients IDs, who belong to the cluster.
clusters_with_ids <- function(patients_ids, clusters_with_indexes) {
  clusters <- vector(mode = "list", length = max(clusters_with_indexes))
  for (i in seq(along = clusters_with_indexes)) {
    cluster <- clusters_with_indexes[i]
    clusters[[cluster]] <- append(clusters[[cluster]], patients_ids[i])
  }
  return(clusters)
}

#' Create a plot which has a plot for each cluster that showcases occurrences of different drugs per patient.
#'
#' @param data: dataframe
#' @param clusters: vector of lists, every list represents one cluster with IDs
#' @param drugs: string vector which contains drug codes
#' @return plot of drug occurrences in clusters
plot_clusters_drugs_count <- function(data, clusters, drugs, norm = FALSE) {
  plots <- list()
  for (i in seq(along = clusters)) {
    cluster <- clusters[[i]]

    drugs_count <- countDrugsInCluster(data, cluster, drugs)
    if (norm) drugs_count <- normalize(drugs_count)

    title <- paste(c("Klaster", i), collapse = " ")

    # If all the values in matrix are the same, we have to add breaks manually
    breaks <- NA
    if (max(drugs_count) == min(drugs_count)) {
      breaks <- seq(0, 5, length.out = 100)
    }
    # If the cluster has only one element, the clustering can not be done
    if (length(cluster) == 1) {
      cluster_cols <- FALSE
    } else {
      row_dendrogram <- flashClust(dist(drugs_count))
      cluster_cols <- row_dendrogram$order
    }

    results <- pheatmap(drugs_count, cluster_rows = FALSE,
                        cluster_cols = cluster_cols, labels_row = drugs, main = title,
                        breaks = breaks)

    plots[[i]] <- results[[4]]
  }
  plot <- grid.arrange(arrangeGrob(grobs = plots, ncol = 2))
  return(plot)
}

#' Get the most frequent drug in the cluster.
#'
#' The most frequent drug per patient is taken into account.
#'
#' @param data: dataframe
#' @param cluster: numeric vector which has patients IDs
#' @return string of most frequent drug code in given cluster
get_most_frequent_drug <- function(data, cluster) {
  frequent_drugs <- c()
  for (i in seq(along = cluster)) {
    patient <- cluster[i]
    filtered_data <- data %>% filter(SUBJECT_ID == patient)
    frequent_drug <- filtered_data %>%
                      count(DRUG_SOURCE_VALUE) %>%
                      top_n(1, n) %>%
                      pull(DRUG_SOURCE_VALUE)
    frequent_drugs <- append(frequent_drugs, frequent_drug)
  }
  frequency_table <- table(frequent_drugs)
  return(names(frequency_table[which.max(frequency_table)]))
}

#' Gets values from certain indexes in matrix.
#'
#' Expects for matrix to be symmetrical and removes duplicates.
#'
#' @param matrix: matrix
#' @param indexes: numeric vector
#' @return vector with matrix values where indexes intersect without duplicates
get_distances <- function(matrix, indexes) {
  submatrix <- matrix[indexes, indexes]
  return(submatrix[row(submatrix) > col(submatrix)])
}

#' Gets average standard deviations in clusters if there are 1 to n clusters.
#'
#' @param matrix: matrix
#' @param max_cluster: integer n > 0, maximum number of clusters
#' @return vector with average standard deviations in clusters. In index 1 is average SD with 1 cluster, in index 2 with 2 clusters and so on.
get_sd_in_clusters <- function(matrix, max_clusters = 10) {
  max_clusters <- min(max_clusters, nrow(matrix))
  average_sds <- c()

  clustering_result <- flashClust(dist(matrix))
  for (i in seq(max_clusters)) {
    clusters <- cutree(clustering_result, i)
    sds_in_clusters <- c()
    for (j in seq(i)) {
      cluster <- which(clusters == j)
      if (length(cluster) > 2) {
        sd_in_cluster <- sd(get_distances(matrix, cluster))
        sds_in_clusters <- append(sds_in_clusters, sd_in_cluster)
      }
    }
    average_sd <- sum(sds_in_clusters) / length(sds_in_clusters)
    average_sds <- append(average_sds, average_sd)
  }
  return(average_sds)
}

#' Plots the elbow for matrix.
#'
#' @param matrix: matrix
#' @return plot which shows average standard deviation in clusters if matrix is clustered into 1, 2, 3, ... clusters
plot_elbow <- function(matrix) {
  y_axis <- get_sd_in_clusters(matrix)
  x_axis <- c(seq_along(y_axis))
  plot(x_axis, y_axis, type = "l", ylab = "Standardhälve", xlab = "Klastrite arv")
}

#' Count occurrences of all drugs in patients records.
#'
#' @param data: dataframe
#' @param id: integer
#' @return dataframe which shows count of all different drugs written to patient
count_all_patient_drugs <- function(data, id) {
  patient_data <- subset(data, SUBJECT_ID == id)
  drugs_count <- as.data.frame(table(patient_data$DRUG_SOURCE_VALUE))
  colnames(drugs_count) <- c("drug_code", "count")
  return(drugs_count)
}

#' Get the count of a drug from data.
#'
#' @param data: dataframe of drug codes with their frequencies
#' @param drug: string
#' @return count of the drug in data, 0 if drug is not in data
get_drug_count <- function(data, drug) {
  drug_row <- data[data$drug_code == drug, ]
  ifelse(length(drug_row$count) == 1, return(drug_row$count), return(0))
}

#' Find additional data about clusters.
#'
#' @param data: dataframe
#' @param clusters: vector with lists for every cluster
#' @return dataframe with additional data about every cluster
get_additional_clusters_data <- function(data, clusters, patients_data) {
  clusters_data <- data.frame()
  max_traj_len <- length(patients_data$drugs) / length(patients_data$IDs)

  for (i in seq(along = clusters)) {
    ids_indexes <- match(clusters[[i]], patients_data$IDs)
    traj_len_sum <- 0
    for (index in ids_indexes) {
      traj<- patients_data$drugs[((index - 1) * max_traj_len + 1):(index * max_traj_len)]
      traj_len_sum <- traj_len_sum + length(traj[! traj %in% c(-1)])
    }
    traj_len_avgs <- round(traj_len_sum / length(clusters[[i]]), 2)
    clusters_data <- rbind(clusters_data, c(i, length(clusters[[i]]), get_most_frequent_drug(data, clusters[[i]]),traj_len_avgs))
  }
  colnames(clusters_data) <- c("Klastri number", "Patsientide arv", "Sagedaseim ravim", "Keskmine trajektoori pikkus")
  return(clusters_data)
}

#' Compares all patients drug trajectories with each other
#'
#' @param patient_data: list with patients IDs and trajectories
#' @param punishments: numeric vector
#' @param dtw: boolean if we want to use DTW algorithm
#' @return symmetrical matrix which values are distances between patients trajectories
compare_trajectories <- function(patients_data, punishments, dtw = FALSE) {
  if (dtw) return(compareDrugTrajectoriesWithDTW(patients_data, punishments))
  return(compareDrugTrajectories(patients_data, punishments))
}

#' Collect data about patients.
#'
#' Finds drug trajectory for every patient.
#'
#' @param data: dataframe
#' @param ids: numeric vector
#' @param n: integer, n > 0 - length of drug trajectory
#' @return list with vectors named "IDs" and "drugs". "IDs" hold information about patients IDs and "drugs" contains n-length drug trajectory per patient
collect_patient_data <- function(data, ids, n) {
  return(collectPatientData(data, ids, n))
}

#' Converts binary code to an image file.
#'
#' @param bin: binary code of image
#' @param name: name of the created file
bin_to_img <- function(bin, name) {
  writeBin(bin, paste(name, ".png", sep = ""))
}

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
