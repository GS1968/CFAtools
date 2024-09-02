# Load necessary libraries
library(lavaan)
library(ggplot2)
library(gridExtra)

# Function to perform CFA with non-parametric bootstrapping and variability analysis
perform_cfa_with_bootstrap <- function(data, model, estimator = "ML", data_type = "continuous", n_bootstrap = 1000) {
  
  # Check the data type and set the appropriate lavOptions
  if (data_type == "ordered") {
    data <- as.data.frame(lapply(data, as.ordered))  # Convert all variables to ordered factors
  } else if (data_type == "dichotomous") {
    data <- as.data.frame(lapply(data, as.numeric))  # Convert all variables to numeric
  }
  
  # Fit the CFA model on the original data
  fit <- cfa(model, data = data, estimator = estimator)
  
  # Extract fit indices from the original model
  fit_indices <- fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "gfi", "agfi", "rmsea", "srmr", "aic", "bic"))
  
  # Calculate CGFI for the original data
  GFI <- fit_indices["gfi"]
  dfT <- fit_indices["df"]
  k <- ncol(data)  # Number of observed variables
  N <- nrow(data)  # Sample size
  CGFI <- GFI + (2 * (1 - (2 * dfT / (k * (k + 1)))) / N)
  
  # Initialize lists to store bootstrap fit indices
  bootstrap_fit_indices <- list(chisq = numeric(n_bootstrap), df = numeric(n_bootstrap), pvalue = numeric(n_bootstrap),
                                cfi = numeric(n_bootstrap), tli = numeric(n_bootstrap), gfi = numeric(n_bootstrap),
                                agfi = numeric(n_bootstrap), rmsea = numeric(n_bootstrap), srmr = numeric(n_bootstrap),
                                aic = numeric(n_bootstrap), bic = numeric(n_bootstrap), cgfi = numeric(n_bootstrap))
  
  # Perform bootstrapping
  for (i in 1:n_bootstrap) {
    # Resample data with replacement
    bootstrap_sample <- data[sample(1:N, replace = TRUE), ]
    
    # Fit the CFA model on the bootstrap sample
    bootstrap_fit <- cfa(model, data = bootstrap_sample, estimator = estimator)
    
    # Extract fit indices from the bootstrap model
    bootstrap_fit_indices$chisq[i] <- fitMeasures(bootstrap_fit, "chisq")
    bootstrap_fit_indices$df[i] <- fitMeasures(bootstrap_fit, "df")
    bootstrap_fit_indices$pvalue[i] <- fitMeasures(bootstrap_fit, "pvalue")
    bootstrap_fit_indices$cfi[i] <- fitMeasures(bootstrap_fit, "cfi")
    bootstrap_fit_indices$tli[i] <- fitMeasures(bootstrap_fit, "tli")
    bootstrap_fit_indices$gfi[i] <- fitMeasures(bootstrap_fit, "gfi")
    bootstrap_fit_indices$agfi[i] <- fitMeasures(bootstrap_fit, "agfi")
    bootstrap_fit_indices$rmsea[i] <- fitMeasures(bootstrap_fit, "rmsea")
    bootstrap_fit_indices$srmr[i] <- fitMeasures(bootstrap_fit, "srmr")
    bootstrap_fit_indices$aic[i] <- fitMeasures(bootstrap_fit, "aic")
    bootstrap_fit_indices$bic[i] <- fitMeasures(bootstrap_fit, "bic")
    
    # Calculate CGFI for the bootstrap sample
    bootstrap_GFI <- bootstrap_fit_indices$gfi[i]
    bootstrap_dfT <- bootstrap_fit_indices$df[i]
    bootstrap_fit_indices$cgfi[i] <- bootstrap_GFI + (2 * (1 - (2 * bootstrap_dfT / (k * (k + 1)))) / N)
  }
  
  # Calculate the mean, standard deviation, and confidence intervals of the bootstrap fit indices
  bootstrap_means <- sapply(bootstrap_fit_indices, mean, na.rm = TRUE)
  bootstrap_sd <- sapply(bootstrap_fit_indices, sd, na.rm = TRUE)
  bootstrap_ci_2.5 <- sapply(bootstrap_fit_indices, function(x) quantile(x, 0.025, na.rm = TRUE))
  bootstrap_ci_95 <- sapply(bootstrap_fit_indices, function(x) quantile(x, 0.975, na.rm = TRUE))
  
  # Capitalized names for the fit indices
  fit_names <- c("Chi-Square", "DF", "P-Value", "CFI", "TLI", "GFI", "AGFI", "RMSEA", "SRMR", "AIC", "BIC", "CGFI")
  
  # Create a table with capitalized names, original, bootstrap, and confidence intervals
  fit_table <- data.frame(
    Measure = fit_names,
    Original = round(c(fit_indices["chisq"], fit_indices["df"], fit_indices["pvalue"], fit_indices["cfi"],
                       fit_indices["tli"], fit_indices["gfi"], fit_indices["agfi"], fit_indices["rmsea"],
                       fit_indices["srmr"], fit_indices["aic"], fit_indices["bic"], CGFI), 3),
    Bootstrap_Mean = round(bootstrap_means, 3),
    Bootstrap_SD = round(bootstrap_sd, 3),
    CI_2.5 = round(bootstrap_ci_2.5, 3),  # Add 2.5% confidence interval
    CI_95 = round(bootstrap_ci_95, 3)  # Add 95% confidence interval
  )
  
  # Print the table with variability information
  print(fit_table)
  
  # Extract parameter estimates if needed
  parameter_estimates <- parameterEstimates(fit)
  print(parameter_estimates)
  
  # Function to create density plots with customizations
  create_density_plot <- function(data, title, color) {
    ggplot(data.frame(x = data), aes(x = x)) +
      geom_density(color = color, size = 1.2) +  # Color the outline only
      ggtitle(paste("Density Plot of Bootstrap", title)) +  # Add custom title
      geom_vline(aes(xintercept = quantile(data, 0.025, na.rm = TRUE)), linetype = "dashed", color = "black", size = 1) +  # 2.5% CI line
      geom_vline(aes(xintercept = quantile(data, 0.975, na.rm = TRUE)), linetype = "dashed", color = "black", size = 1) +  # 97.5% CI line
      geom_vline(aes(xintercept = mean(data, na.rm = TRUE)), color = "blue", size = 1) +  # Mean line
      theme_minimal() +
      theme(axis.title.x = element_blank(),  # Remove x-axis title
            axis.title.y = element_blank(),  # Remove y-axis title
            plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold title
            axis.text = element_text(face = "bold"),  # Bold axis text
            axis.text.x = element_text(hjust = 0.5, face = "bold"),  # Center and bold x-axis labels
            axis.text.y = element_text(face = "bold"),  # Bold y-axis labels
            panel.grid.major = element_line(size = 0.8),  # Thicker grid lines
            panel.grid.minor = element_line(size = 0.6))  # Thicker minor grid lines
  }
  
  # Define colors for the plots
  colors <- c("red", "blue", "green", "purple", "orange", "cyan", "pink", "brown")
  
  # Create density plots for the specified fit indices
  plot_list <- list(
    create_density_plot(bootstrap_fit_indices$chisq, "Chi-Square", colors[1]),
    create_density_plot(bootstrap_fit_indices$cfi, "CFI", colors[2]),
    create_density_plot(bootstrap_fit_indices$tli, "TLI", colors[3]),
    create_density_plot(bootstrap_fit_indices$gfi, "GFI", colors[4]),
    create_density_plot(bootstrap_fit_indices$agfi, "AGFI", colors[5]),
    create_density_plot(bootstrap_fit_indices$rmsea, "RMSEA", colors[6]),
    create_density_plot(bootstrap_fit_indices$srmr, "SRMR", colors[7]),
    create_density_plot(bootstrap_fit_indices$cgfi, "CGFI", colors[8])
  )
  
  # Arrange the plots in a grid (2 plots per row, 4 rows in total)
  grid.arrange(grobs = plot_list, ncol = 2)
}
