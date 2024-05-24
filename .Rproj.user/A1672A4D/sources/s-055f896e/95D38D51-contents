library(pheatmap)
library(flashClust)
library(dplyr)
library(grDevices)
library(grid)

background_color <- "#e9f4fb"

#' Gets the count of most extreme value in a column.
#'
#' Can find the maximal and minimal frequency.
#'
#' @param data: dataframe
#' @param column: string
#' @param max_value: boolean if we want to get the maximal frequency
#' @return maximal or minimal frequency in data column
get_column_extreme_value_count <- function(data, column, max_value = TRUE) {
  ifelse(max_value, return(max(table(data[column]))), return(min(table(data[column]))))
}

#' Gets every patients' records count.
#'
#' @param data: dataframe
#' @return vector with patients' records count
get_patients_records_count <- function(data) {
  patients_by_id <- data %>%
    group_by(SUBJECT_ID) %>%
    dplyr::summarise(records_count = n())
  return(patients_by_id$records_count)
}

#' Gets mean of patients' records count.
#'
#' @param data: dataframe
#' @return patients' records count mean as numeric
get_records_mean <- function(data) {
  return(mean(get_patients_records_count(data)))
}

#' Gets records range for data, that includes at least 50% of patients.
#'
#' @param data: dataframe
#' @return vector with records range
get_records_range <- function(data) {
  patients_records <- sort(get_patients_records_count(data))
  patient_count <- 0.5 * length(patients_records)
  # Gets the patients from the middle
  records_range_min <- patients_records[as.integer((length(patients_records) - patient_count) / 2)]
  records_range_max <- patients_records[as.integer(((length(patients_records) - patient_count) / 2) + patient_count)]
  return(c(records_range_min, records_range_max))
}

#' Rounds given value to first upper or lower tens.
#'
#' Can round up or down.
#'
#' @param value: numeric
#' @param round-up: boolean if we want to round up
#' @return numeric rounded value
get_rounded_value <- function(value, round_up = TRUE) {
  ifelse(round_up, f <- ceiling, f <- floor)
  return(plyr::round_any(value, 10, f = f))
}

server <- function(input, output, session) {
  example_file_path <- "../tmp/datasets/example_data.csv"

  data <- reactive({
    if (!is.null(input$data_file)) {
      return(read.csv(input$data_file$datapath))
    }
    return(read.csv(example_file_path))
  })

  observe({
    if (!is.null(data())) {
      max_records <- get_rounded_value(get_column_extreme_value_count(data(), "SUBJECT_ID", TRUE), TRUE)
      min_records <- get_rounded_value(get_column_extreme_value_count(data(), "SUBJECT_ID", FALSE), FALSE)
      max_trajectory_length <- get_rounded_value(max_records * 0.5, FALSE)
      updateSliderInput(session, "records_range", min = min_records, max = max_records)
      updateSliderInput(session, "trajectory_length", max = max_trajectory_length)
    }
  })


  # Remembers previous user inputs
  previous_inputs <- reactiveValues(data_file = "", popular_drugs_amount = -1, records_range = c(-1, -1), trajectory_length = -1,
                                    level_one = -1, level_two = -1, level_three = -1, level_four = -1, level_five = -1,
                                    clusters_amount = -1, use_dtw_algorithm = FALSE)
  # Remembers previous values of variables
  previous_values <- reactiveValues(first_run = TRUE, data = NULL, most_frequent_drugs = NULL, patients = NULL,
                                    patients_sample = NULL, filtered_data = NULL,
                                    patients_data = NULL, punishments = NULL, compare_matrix = NULL,
                                    clustered_trajectories = NULL, clusters_plots = NULL, clusters = NULL)

  observeEvent(input$submit, {
    output$clustered_trajectories <- renderPlot(NULL)
    output$elbow_plot <- renderPlot(NULL)
    output$clusters_plots <- renderPlot(NULL)
    output$clusters_data <- renderPlot(NULL)
    # Find which inputs have changed
    changed_inputs <- c()
    for (name in names(input)) {

      if (name == "data_file") {
        prev <- previous_inputs[[name]]
        if (is.null(input[[name]])) {
          previous_inputs[[name]] <- example_file_path
        } else {
          previous_inputs[[name]] <- input[[name]]$datapath
        }
        if (!setequal(prev, previous_inputs[[name]])) {
          changed_inputs <- append(changed_inputs, name)
          previous_values[["first_run"]] <- TRUE
        }
      # If the data file is ran through the first time, then the records range and trajectory length are calculated a default value
      } else if (name == "records_range" & previous_values[["first_run"]] == TRUE) {
        # Calculates records range default value
        changed_inputs <- append(changed_inputs, name)
        default_records_range <- get_records_range(data())
        previous_inputs[[name]] <- default_records_range
        updateSliderInput(session, "records_range", value = default_records_range)

      } else if (name == "trajectory_length" & previous_values[["first_run"]] == TRUE) {
        # Calculates trajectory length default value
        changed_inputs <- append(changed_inputs, name)
        default_trajectory_length <- round(get_records_mean(data()))
        previous_inputs[[name]] <- default_trajectory_length
        updateSliderInput(session, "trajectory_length", value = default_trajectory_length)

      } else if (!is.null(previous_inputs[[name]])) {
        if (!setequal(previous_inputs[[name]], input[[name]])){
          changed_inputs <- append(changed_inputs, name)
          previous_inputs[[name]] <- input[[name]]
        }
      }
    }
    previous_values[["first_run"]] <- FALSE

    # Keep in memory which values have changed
    changed_values <- c()
    if ("data_file" %in% changed_inputs) {
      previous_values[["data"]] <- read.csv(previous_inputs[["data_file"]])
      changed_values <- append(changed_values, "data")
    }

    if ("data" %in% changed_values | "records_range" %in% changed_inputs) {
      previous_values[["patients"]] <- get_filtered_patients(previous_values[["data"]], previous_inputs[["records_range"]])
      # Checks if we have enough patients to cluster and run through the work flow.
      # If not, the records range is widened.
      while (length(previous_values[["patients"]]) < 10) {
        new_records_range <- c(previous_inputs[["records_range"]][1] - 1, previous_inputs[["records_range"]][2] + 1)
        previous_inputs[["records_range"]] <- new_records_range
        previous_values[["patients"]] <- get_filtered_patients(previous_values[["data"]], previous_inputs[["records_range"]])
        updateSliderInput(session, "records_range", value = previous_inputs[["records_range"]])
      }
      previous_values[["patients_sample"]] <- sample(previous_values[["patients"]], min(length(previous_values[["patients"]]), 600))
      previous_values[["filtered_data"]] <- previous_values[["data"]] %>% filter(SUBJECT_ID %in% previous_values[["patients_sample"]])
      changed_values <- append(changed_values, "patients")
    }

    if ("data" %in% changed_values | "patients" %in% changed_values | "trajectory_length" %in% changed_inputs) {
      previous_values[["patients_data"]] <- collect_patient_data(previous_values[["filtered_data"]], previous_values[["patients_sample"]], previous_inputs[["trajectory_length"]])
      changed_values <- append(changed_values, "patients_data")
    }

    levels <- c("level_one", "level_two", "level_three", "level_four", "level_five")
    if (TRUE %in% (levels %in% changed_inputs)) {
      previous_values[["punishments"]] <- c()
      # Checks if any of the values are below 0. If such values are found, they are replaced with 0.
      for (i in seq(levels)) {
        current_level <- levels[i]
        current_level_value <- previous_inputs[[current_level]]
        if (current_level_value < 0) {
          current_level_value <- 0
          updateNumericInput(session, current_level, value = 0)
        }
        previous_values[["punishments"]] <- append(previous_values[["punishments"]], current_level_value)
      }
      changed_values <- append(changed_values, "punishments")
    }

    if ("data" %in% changed_values | "patients_data" %in% changed_values | "punishments" %in% changed_values | "use_dtw_algorithm" %in% changed_inputs) {
      previous_values[["compare_matrix"]] <- compare_trajectories(previous_values[["patients_data"]], previous_values[["punishments"]], previous_inputs[["use_dtw_algorithm"]])
      changed_values <- append(changed_values, "compare_matrix")
    }

    if ("data" %in% changed_values | "compare_matrix" %in% changed_values | "clusters_amount" %in% changed_inputs) {
      row_dendrogram <- flashClust(dist(previous_values[["compare_matrix"]]))
      previous_values[["clustered_trajectories"]] <- pheatmap(previous_values[["compare_matrix"]], cutree_rows = previous_inputs[["clusters_amount"]], cluster_rows = row_dendrogram$order, cluster_cols = FALSE)
      changed_values <- append(changed_values, "clustered_trajectories")
    }

    if("data" %in% changed_values | "patients_data" %in% changed_values | "clustered_trajectories" %in% changed_values | "clusters_amount" %in% changed_inputs) {
      previous_values[["clusters"]] <- clusters_with_ids(previous_values[["patients_data"]]$IDs, cutree(previous_values[["clustered_trajectories"]]$tree_row, previous_inputs[["clusters_amount"]]))
      changed_values <- append(changed_values, "clusters")
    }

    if ("data" %in% changed_values | "popular_drugs_amount" %in% changed_inputs) {
      previous_values[["most_frequent_drugs"]] <- get_most_frequent_values(previous_values[["data"]], "DRUG_SOURCE_VALUE", previous_inputs[["popular_drugs_amount"]])
      changed_values <- append(changed_values, "most_frequent_drugs")
    }

    previous_values[["clusters_plots"]] <- reactive(plot_clusters_drugs_count(previous_values[["filtered_data"]], previous_values[["clusters"]], previous_values[["most_frequent_drugs"]], input$normalize))

    output$clustered_trajectories <- renderPlot({previous_values[["clustered_trajectories"]]}, bg = background_color)
    output$elbow_plot <- renderPlot({plot_elbow(previous_values[["compare_matrix"]])}, bg = "#f5f5f5")
    output$clusters_plots <- renderPlot({previous_values[["clusters_plots"]]()}, bg = background_color)
    output$clusters_data <- renderTable({get_additional_clusters_data(previous_values[["filtered_data"]], previous_values[["clusters"]], previous_values[["patients_data"]])}, width = "100%", striped = TRUE, bg = "white")

  }, ignoreInit = FALSE, ignoreNULL = FALSE)

  output$download_results <- downloadHandler(filename = function() {paste("results_", Sys.Date(), ".rds", sep = "")}, content = function(file) {
    # Creates images of graphs and then converts them to binary. Created image file is deleted.
    # This makes it easier to share graphs inside R object.
    graph_1 <- png(filename = "clustering_graph.png")
    print(previous_values[["clustered_trajectories"]])
    dev.off()
    graph_1_bin <- readBin("clustering_graph.png", "raw", file.info("clustering_graph.png")$size)
    file.remove("clustering_graph.png")

    graph_2 <- png(filename = "clusters_graph.png")
    grid.draw(previous_values[["clusters_plots"]]())
    dev.off()
    graph_2_bin <- readBin("clusters_graph.png", "raw", file.info("clusters_graph.png")$size)
    file.remove("clusters_graph.png")

    results <- list(parameters = reactiveValuesToList(previous_inputs),
                    clustering_graph = graph_1_bin,
                    clusters_graph = graph_2_bin,
                    clusters = NULL)
    # If some of the patients were left out in the application, then for final result we cluster all patients.
    if (setequal(previous_values[["patients"]], previous_values[["patients_sample"]])) {
      results$clusters <- previous_values[["clusters"]]
    } else {
      filtered_data <- previous_values[["data"]] %>% filter(SUBJECT_ID %in% previous_values[["patients"]])
      all_patients_data <- collect_patient_data(filtered_data, previous_values[["patients"]], previous_inputs[["trajectory_length"]])
      compare_matrix <- compare_trajectories(all_patients_data, previous_values[["punishments"]], previous_inputs[["use_dtw_algorithm"]])
      all_patients_clusters <- clusters_with_ids(all_patients_data$IDs, cutree(flashClust(dist(compare_matrix)), previous_inputs[["clusters_amount"]]))
      results$clusters <- all_patients_clusters
    }
    saveRDS(results, file)
  })

  output$download_parameters <- downloadHandler(filename = function() {paste("parameters_", Sys.Date(), ".csv", sep = "")}, content = function(file) {
    used_parameters <- reactiveValuesToList(previous_inputs)
    used_parameters$records_range <- paste(used_parameters$records_range, collapse = "-")
    write.csv(as.data.frame(used_parameters), file, row.names = FALSE)
  })
}
