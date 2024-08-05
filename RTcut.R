# Load necessary libraries
library(dplyr)
library(ggplot2)

# Define constants
Chance <- 1/4  # the random guessing chance correct for an item with 4 options
ULimit <- 2000  # the upper limit for threshold; user choice
LLimit <- 0    # the lower limit for threshold; user choice

# Define the function to calculate threshold value for item j
threshold.value <- function(j, x, y) {
  # Create contingency table for response time and response (0/1)
  item <- table(x[, j], y[, j])
  item <- as.matrix(item)
  
  # Extract response times
  S <- as.numeric(row.names(item))
  
  # Calculate P+ at each response time
  Pplus <- item[, 2] / apply(item, 1, sum)
  
  # Calculate cumulative P+
  cumP <- Pplus
  K <- length(item[, 2])
  for (k in 2:K) {
    cumP[k] <- sum(item[1:k, 2]) / sum(item[1:k, ])
  }
  
  # Identify indices where cumulative P+ is less than the Chance
  ind <- (cumP < Chance)
  value.p <- cumP[ind]
  
  # Find the last index where cumP is less than Chance
  LL <- which(ind, arr.ind = TRUE)
  LL <- as.vector(LL)
  
  # Check if LL is empty
  if (length(LL) == 0) {
    value.T <- ULimit  # or any default value within limits
  } else {
    value.T <- LL[length(LL)]
  }
  
  # Adjust threshold value within the specified limits
  if (value.T > ULimit) {
    value.T <- ULimit
  } else if (value.T < LLimit) {
    value.T <- LLimit
  }
  
  # Return the threshold value
  return(drop(value.T))
}

# Define the function to estimate cutoff values for each item
estimate_cutoffs <- function(response_times) {
  cutoff_values_plus_1sd <- apply(response_times, 2, function(rt) {
    mean_rt <- mean(rt, na.rm = TRUE)
    sd_rt <- sd(rt, na.rm = TRUE)
    cutoff <- mean_rt + (1 * sd_rt)
    return(cutoff)
  })
  
  mean_values <- apply(response_times, 2, function(rt) {
    mean_rt <- mean(rt, na.rm = TRUE)
    return(mean_rt)
  })
  
  return(list(plus_1sd = cutoff_values_plus_1sd, mean = mean_values))
}

# Function to read response time data from a .csv file and calculate thresholds
calculate_thresholds_from_csv <- function(response_time_csv, response_data_csv, participant_row) {
  # Read response time data from .csv file
  x <- read.csv(response_time_csv, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  
  # Read response (0/1) data from another .csv file
  y <- read.csv(response_data_csv, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  
  # Convert data frames to matrices
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  # Calculate the threshold values for all items
  num_items <- ncol(x)
  threshold_values <- sapply(1:num_items, function(j) threshold.value(j, x, y))
  print(threshold_values)
  
  # Select the specified participant for plotting
  participant_response_times <- x[participant_row, ]
  
  # Calculate the cutoff values for each item
  cutoff_values <- estimate_cutoffs(x)
  
  # Create a data frame for plotting
  threshold_data <- data.frame(Item = 1:num_items, Threshold = threshold_values)
  participant_data <- data.frame(Item = 1:num_items, ResponseTime = participant_response_times)
  cutoff_data_plus_1sd <- data.frame(Item = 1:num_items, Cutoff = cutoff_values$plus_1sd)
  mean_data <- data.frame(Item = 1:num_items, Mean = cutoff_values$mean)
  
  # Remove rows with NA values to avoid plotting issues
  threshold_data <- threshold_data[complete.cases(threshold_data), ]
  participant_data <- participant_data[complete.cases(participant_data), ]
  cutoff_data_plus_1sd <- cutoff_data_plus_1sd[complete.cases(cutoff_data_plus_1sd), ]
  mean_data <- mean_data[complete.cases(mean_data), ]
  
  # Determine the limits for the y-axis
  y_limit <- max(c(threshold_data$Threshold, participant_data$ResponseTime, 
                   cutoff_data_plus_1sd$Cutoff, mean_data$Mean), na.rm = TRUE)
  
  # Plot the threshold values for all items along with the response times for a single participant
  ggplot() +
    geom_line(data = threshold_data, aes(x = Item, y = Threshold, color = "CUMP RT"), size = 1) +
    geom_point(data = threshold_data, aes(x = Item, y = Threshold, color = "CUMP RT"), size = 3) +
    geom_line(data = participant_data, aes(x = Item, y = ResponseTime, color = "Response Time"), linetype = "dotted", size = 1) +
    geom_point(data = participant_data, aes(x = Item, y = ResponseTime, color = "Response Time"), shape = 17, size = 3) +
    geom_line(data = cutoff_data_plus_1sd, aes(x = Item, y = Cutoff, color = "1SD RT"), linetype = "dashed", size = 1) +
    geom_point(data = cutoff_data_plus_1sd, aes(x = Item, y = Cutoff, color = "1SD RT"), shape = 22, size = 3, fill = "green") +
    geom_line(data = mean_data, aes(x = Item, y = Mean, color = "Mean RT"), size = 1.5) +
    xlab("Item") +
    ylab("Response Times in Seconds") +
    scale_color_manual(name = "Legend", values = c("CUMP RT" = "blue", "Response Time" = "red", "1SD RT" = "green", "Mean RT" = "black")) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.line = element_line(color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey", linetype = "dotted"),
      panel.grid.minor.y = element_blank(),
      axis.ticks.length = unit(.25, "cm"),
      axis.ticks.x = element_line(size = 0.5, color = "black"),
      axis.ticks.y = element_line(size = 0.5, color = "black"),
      axis.ticks.length.x = unit(0.3, "cm"),
      axis.ticks.length.y = unit(0.3, "cm")
    ) +
    scale_x_continuous(breaks = 1:num_items, labels = paste("Item", 1:num_items)) +
    scale_y_continuous(limits = c(0, y_limit))
}
