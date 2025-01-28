# Libraries ---------------------------------------------------------------
# Data Manipulation and Wrangling 
library(tidyverse)   # v2.0.0   Comprehensive package for data wrangling, visualization, and manipulation.
# Statistical Analysis 
library(gee)         # v4.13-27 Generalized estimating equations.
library(geepack)     # v1.3.12  Additional tools for GEE modeling.
library(car)         # v3.1-3   Companion to Applied Regression; includes skewness tests.
library(moments)     # v0.14.1  Skewness and kurtosis calculations.
library(contrast)    # v0.24.2  Pairwise and custom contrasts (e.g., `apc()` for contrast matrices).
library(bmet)        # v0.1.0   Tools for Bayesian modeling and metabolic analysis
library(fmsb)        # v0.7.6   For pairwise fisher tests
# Visualization 
library(ggplot2)     # v3.5.1   Flexible data visualization package.
library(ggsignif)    # v0.6.4   Add significance brackets to ggplot2.
library(ggstatsplot) # v0.12.5  Enhanced ggplots with statistical annotations.
library(ggbeeswarm)  # v0.7.2   Creates bee swarm plots.
library(grid)        # v4.3.2   Lower-level graphics system for arranging plots.
library(gridExtra)   # v2.3     Allows arranging multiple plots using grid.
library(scales)      # v1.3.0   Tools for scaling and transformation in plots.
library(svglite)     # v2.1.3   For saving ggplots as SVG files.
#################

# Set Experiment-Specific Variables for Analysis ---------------------------------------------------------------
#Working Directory Containing Results CSVs
data_dir <- "/Volumes/T7_Shield/Imaging (monolayer)/241219/TimeStacks"
expt_name <- "2401219"

#define replicate wells for each experimental condition (aka wellgroups)
wellgroup_names <- c(
  "Veh",
  "SST",
  "YM - Veh",
  "YM - SST"
)

#define wellgroup identities
wellgroup_defs <- list(
  c("D07", "E07", "D08", "E08", "D09", "E09"), # Condition 1 (columns 7 to 9)
  c("D10", "E10", "D11", "E11", "D12", "E12"), # Condition 2 (columns 10 to 12)
  c("D13", "E13", "D14", "E14", "D15", "E15"), # Condition 3 (columns 13 to 15)
  c("D16", "E16", "D17", "E17", "D18", "E18") # Condition 4 (columns 6 to 18)
)

names(wellgroup_defs) <- wellgroup_names

control_wellgroup <- "Veh" #define 'control' wellgroup for showing pairwise comparisons

baseline_period <- c(0, 40) #'baseline' period start and end times

drug_times = c(41, 81) #define drug addition timepoints (leave empty for no lines)

#Define start and end times for each time period to be measured for statistical analysis
time_period_defs <- data.frame(StartTime = c(20, 41, 60), 
                               EndTime = c(40, 46, 80)
)

#trace plot parameters
y_limits = c(-.3, 0.5)
x_limits = c(20, 80) #NA is full range

#stats plot parameters
periods_to_analyze <-c("P3")
set.cex <- 0.9
set.stats.alpha <- 0.5
hide_comparisons <- "FALSE"

#Define parameters for exclusion of cells from analysis 
baseTooHigh <- 0.4
tooLow <- -0.7
tooHigh <- 0.9

#define individual plot sizes for printing
plot_width <- 4
plot_height <- 3
stats_plot_width <- 2
stats_plot_height <- 3
#################
#################
##Functions##########################
normalize_column <- function(column) {
  #baseline <- column[1]
  baseline <- min(column[baseline_period[1]:baseline_period[2]])
  final_value <- column[length(column)]
  normalized_values <- (column - baseline) / (final_value - baseline)
  return(normalized_values)
}
measure_cells_by_period <- function(data, time_period_defs) {
  cell_measurements <- data.frame()
  first_period_volatility <- list()
  first_period_max <- list()
  
  for(i in 1:nrow(time_period_defs)) {
    start_time <- time_period_defs$StartTime[i]
    end_time <- time_period_defs$EndTime[i]
    period_data <- subset(data, Time >= start_time & Time <= end_time)
    
    for(cell_id in unique(period_data$cell)) {
      #subset only data in period
      cell_data <- subset(period_data, cell == cell_id)
      #calculate period max
      max_intensity <- max(cell_data$intensity, na.rm = TRUE)
      #calculate period volatility
      volatility_intensity <- (1/nrow(cell_data))*sum(abs(diff(cell_data$intensity)), na.rm = TRUE)
      
      # Calculate difference in max (P-P1) For first period (assumed baseline) max, recorded in list, then take ratio.
      if (i == 1) {
        first_period_max[[as.character(cell_id)]] <- max_intensity
      }
      if (!is.null(first_period_max[[as.character(cell_id)]])) {
        dif_max <- max_intensity - first_period_max[[as.character(cell_id)]]
      } else {
        dif_max <- NA  # Handle unexpected missing baseline
        warning(paste("Missing baseline max intensity for cell", cell_id, "in period", i))
      }
      
      # Calculate relative volatility (P/P1) For first (assumed baseline) volatility, recorded in list, then take ratio.
      if (i == 1) {
        first_period_volatility[[as.character(cell_id)]] <- volatility_intensity
      }
      if (!is.null(first_period_volatility[[as.character(cell_id)]])) {
        ln_rel_volatility <- log(volatility_intensity / first_period_volatility[[as.character(cell_id)]])
      } else {
        rel_volatility <- NA
        warning(paste("Missing baseline volatility intensity for cell", cell_id, "in period", i))
      }
      
      each_cell_measurement <- data.frame(
        wellgroup = unique(cell_data$wellgroup),
        expt = unique(cell_data$expt),
        well = unique(cell_data$well),
        cell = unique(cell_data$cell),
        tperiod.name = paste0("P", i),
        tperiod.range = paste0(start_time, "-", end_time),
        max = max_intensity,
        dif_max = dif_max,
        volatility = volatility_intensity,
        ln_rel_volatility = ln_rel_volatility
      )
      
      cell_measurements <- rbind(cell_measurements, each_cell_measurement)
    }
  }
  
  return(cell_measurements)
}
find_wellgroup <- Vectorize(function(well) {
  for (group_name in wellgroup_names) {
    if (well %in% wellgroup_defs[[group_name]]) {
      return(group_name)
    }
  }
  return(NA) # Return NA if no group is found
})

createTracePlot <- function(data, title, y_limits = c(NA, NA), x_limits = c(NA, NA), line_x_intercepts = c(), text_scale = 12) {
  average_trace <- aggregate(intensity ~ Time, data = data, FUN = mean)
  x_midpoint <- mean(range(data$Time, na.rm = TRUE))
  nCells <- length(unique(data$cell))

  p <- ggplot() +
    geom_line(data = data,
              aes(x = Time, y = intensity, group = cell),
              color = "black",
              alpha = 0.06) +
    geom_line(data = average_trace,
              aes(x = Time, y = intensity),
              color = "red",
              linewidth = .5) +
    scale_y_continuous(limits = y_limits, breaks = seq(from = 0, to = max(y_limits), by = 0.2)) +  # Define y-axis breaks
    scale_x_continuous(
      name = "Time (min)",
      limits = ifelse(is.na(x_limits), range(data$Time, na.rm = TRUE), x_limits),
      breaks = function(x) seq(from = 0, to = max(x, na.rm = TRUE), by = 20),
      labels = function(x) x / 2  # Convert frames to minutes; adjust as necessary
    ) +
    labs(title = title, y = "Norm. dF/F0") +
    theme_minimal(base_size = text_scale) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "none",
      axis.title.x = element_text(size = rel(.8), face = "bold"),
      axis.title.y = element_text(size = rel(.8), face = "bold"),
      axis.text.x  = element_text(size = rel(1)),
      axis.text.y  = element_text(size = rel(1)),
      plot.title   = element_text(size = rel(1.5), face = "bold"),
      axis.ticks.x = element_line(size = 0.5),  # Add x-axis ticks
      axis.ticks.y = element_line(size = 0.5)   # Add y-axis ticks
    ) +
    annotate("text", x = x_midpoint, y = y_limits[2] - 0.1,
             label = paste("n =", nCells),
             size = text_scale * 0.5,
             colour = "black")
    if (length(line_x_intercepts) > 0) {
    for (x_intercept in line_x_intercepts) {
      p <- p + geom_vline(xintercept = x_intercept, linetype = "dashed", color = "black")
    }
    }

  #report plot name when complete
  print(paste(plotname, "plotted with", nCells, "cells"))

  return(p)
}

run_GEE_and_plot <- function(data,
                             measurement_type,
                             transform_data = FALSE,
                             FDR = FALSE,
                             y_label,
                             period,
                             label_height,
                             hline) {

  # Specify wellgroups to analyze
  if (exists("wellgroups_to_analyze") && length(wellgroups_to_analyze) > 0) {
    ordered_wellgroups <- wellgroups_to_analyze
  } else {
    ordered_wellgroups <- wellgroup_names
  }

  # Filter the measurements data for the specified periods
  data_filtered <- data %>%
    filter(tperiod.name == period) %>%
    filter(wellgroup %in% ordered_wellgroups) %>%
    mutate(wellgroup = factor(wellgroup, levels = ordered_wellgroups))
  data_filtered <- data_filtered[order(data_filtered$well),]

  # skewtest and qqPlot for raw data
  skewtest.raw <- agostino.test(data_filtered[, measurement_type])
  qqPlot(data_filtered[, measurement_type], main = paste("raw", measurement_type,period))
  qplot.raw <- recordPlot()

  if (transform_data == TRUE) {
    if (skewtest.raw$statistic[["skew"]] > 0 & skewtest.raw$p.value < 0.05) { ## if sig. right-skewed, log-transform
      ## note that we can't log-transform negative or zero values so will create a temporary variable that shifts
      ## all values up for the transformation and modeling (here creating a new variable)
      measurement_gee <- paste0("ln_", measurement_type)
      if (min(data_filtered[, measurement_type], na.rm=TRUE) <= 0) {
        data_filtered[, measurement_gee] <- data_filtered[, measurement_type] + abs(min(data_filtered[, measurement_type], na.rm=TRUE)) + 0.001 ## add min plus a tiny value so we don't have 0s
      } else {
        data_filtered[, measurement_gee] <- data_filtered[, measurement_type] ## no need to add tiny value
      }
      data_filtered[, measurement_gee] <- log(data_filtered[, measurement_gee]) #now we can log transform to reduce positive skew

      print("Data is positively skewed. Log transformed.")
    } else if (skewtest.raw$statistic[["skew"]] < 0 & skewtest.raw$p.value < 0.05) { ## if sig. left-skewed, reflect the data and log-transform
      measurement_gee <- paste0("ln_", measurement_type)
      maxraw <- max(data_filtered[, measurement_type], na.rm=TRUE) + 0.001 ## get max plus a tiny value so we don't have 0s
      data_filtered[, measurement_gee] <- maxraw - data_filtered[, measurement_type] ## reflected data by subtracting from max (now pos. skewed)
      ## now we can log transform to reduce positive skew
      data_filtered[, measurement_gee] <- log(data_filtered[, measurement_gee])
      print("Data is negatively skewed. Reflected and log transformed.")
    } else {
      measurement_gee <- paste0("raw_", measurement_type)
      data_filtered[, measurement_gee] <- data_filtered[, measurement_type] ## no skew
      print("Data is not skewed. Not transformed.")
    }
  } else {
    measurement_gee <- paste0("raw_", measurement_type)
    data_filtered[, measurement_gee] <- data_filtered[, measurement_type] ## data already transformed - do nothing
  }

  #look at transformed data
  skewtest.transf <- agostino.test(data_filtered[, measurement_gee])
  qqPlot(data_filtered[, measurement_gee], main = paste("transformed", measurement_type,period))
  qplot.transf <- recordPlot()

  #perform GEE without pairwise comparisons 
  formula_text <- paste(measurement_gee, "~ wellgroup")
  gee_model <- gee(as.formula(formula_text),
                   data = data_filtered,
                   id = as.factor(well),
                   family = gaussian,
                   corstr = "exchangeable")
  gee_sum <- summary(gee_model)
  gee_res <- as.data.frame(gee_sum$coefficients)
  gee_res$robust_p <- 2*pnorm(q=abs(gee_res$`Robust z`), lower.tail=FALSE)
  gee_res <- gee_res[which(rownames(gee_res) != "(Intercept)"), ]
  gee_res <- gee_res %>%
    mutate(Significance = case_when(
      `robust_p` < 0.0001 ~ "****",
      `robust_p` < 0.001  ~ "***",
      `robust_p` < 0.01   ~ "**",
      `robust_p` < 0.05   ~ "*",
      TRUE             ~ "ns"  # Not significant
    ))
  rownames(gee_res) <- gsub("wellgroup", "", rownames(gee_res))

  ##Check Residuals to see if GEE assumptions are met
  res0test <- t.test(gee_model$residuals, mu = 0, alternative = "two.sided")
  print(paste("P-value testing if residual mean=0 (small p-value indicates violation)=",
              round(res0test$p.value,3)))
  qqPlot(gee_model$residuals, main = paste("residuals", measurement_gee,period))
  qplot.residuals <- recordPlot()

  #Execute GEE test
  gee_contrasts <- geese(as.formula(formula_text),
                         data = data_filtered,
                         id = as.factor(well),
                         family = gaussian,
                         corstr = "exchangeable")

  ## perform all contrasts and store raw p-values
  contmat <- apc(length(levels(data_filtered$wellgroup)), labs = paste(levels(data_filtered$wellgroup),"_")) #identity matrix for comparisons 
  contp <- data.frame(group1=character(), group2=character(), p_raw=numeric()) ## empty data frame that will get filled with p-values
  for (w in 1:length(rownames(contmat))) {
    grps <- rownames(contmat)[w]
    grp1 <- str_split_i(grps, "_", 1) ## extract first group name
    grp2 <- str_split_i(grps, "_", 2) ## extract second group name
    grp2 <- sub('.', '', grp2) ## get rid of the "-" at the beginning of 2nd group name
    grp1 <- str_trim(grp1)
    grp2 <- str_trim(grp2)
    contgrp <- contrast(gee_contrasts, list(wellgroup = grp1), list(wellgroup = grp2)) #get P values for contrasts here
    contp[w,] <- c(grp1, grp2, contgrp$Pvalue)
  }

  # ensure P values are numeric
  contp <- contp %>%
    mutate(p_raw = as.numeric(p_raw))

  # adjust p-values for multiple comparisons (Bonferroni)
  if (FDR == TRUE) {
      contp$p_adjusted <- p.adjust(contp$p_raw, method = "bonferroni")
      } else {contp$p_adjusted <- contp$p_raw
      }

  # assign stars for significance levels
  contp <- contp %>%
    mutate(Significance = case_when(
      `p_adjusted` < 0.0001 ~ "****",
      `p_adjusted` < 0.001  ~ "***",
      `p_adjusted` < 0.01   ~ "**",
      `p_adjusted` < 0.05   ~ "*",
      TRUE             ~ "ns"  # Not significant
    ))

  # Create a matrix of P values for easy review
  groups <- levels(data_filtered$wellgroup)
  p_matrix <- matrix(NA, nrow = length(groups), ncol = length(groups),
                     dimnames = list(groups, groups))
  for (i in 1:nrow(contp)) {
    group1 <- contp$group1[i]
    group2 <- contp$group2[i]
    sig <- contp$Significance[i]
    p_matrix[group2, group1] <- sig
  }

  # Calculate quick statistics to use for error bars
  quick_stats <- data_filtered %>%
    group_by(wellgroup) %>%
    summarise(
      meas_mean = mean(!!sym(measurement_type), na.rm = TRUE),
      sd = sd(!!sym(measurement_type), na.rm = TRUE),
      n = n(),
      se = sd / sqrt(n)  # Standard Error
    )

  #make list of significant comparisons for GGsignif to plot
  comparisons <- list()
  significant_indices <- which(contp$Significance != "ns")
  for (i in 1:nrow(contp)) {
    comparisons[[i]] <- c(contp$group1[i], contp$group2[i])
  }
  comparisons_sig <- comparisons[significant_indices]

  #remove top and bottom 10 ourliers from plotted points
  top_bottom_outliers <- c(
    order(data_filtered[, measurement_type], decreasing = TRUE)[1:10],
    order(data_filtered[, measurement_type])[1:10]
  )
  data_to_show <- data_filtered %>%
    filter(!row_number() %in% top_bottom_outliers)


  # Generate the plot with the original data for the violin and jitter
  plot <- ggplot() +
    geom_beeswarm(data = data_to_show, aes(x = wellgroup, y = !!sym(measurement_type), color = wellgroup), shape = 16, size = .9, alpha = set.stats.alpha, cex = set.cex)+
    geom_crossbar(data = quick_stats, aes(x = wellgroup, y = meas_mean, ymin = meas_mean, ymax = meas_mean), width = 0.6) +
    geom_errorbar(data = quick_stats, aes(x = wellgroup, ymin = meas_mean - sd, ymax = meas_mean + sd), width = 0.3) +
    labs(title = paste(measurement_type, period), x = NULL, y = y_label) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x  = element_text(size = 11, color = "black", face = "bold", angle = 20, hjust = 1),  # Set x-axis tick label color to black
      axis.text.y  = element_text(size = 12, color = "black", face = "bold"),  # Set y-axis tick label color to black
      axis.title = element_text(face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(),
      axis.line=element_line(),
    )
  if (length(comparisons_sig) > 0 && (!exists("hide_comparisons") || hide_comparisons != "TRUE")) {
    plot <- plot + geom_signif(data = data_to_show,
                               aes(x = wellgroup, y = !!sym(measurement_type), group = wellgroup),
                               comparisons = comparisons_sig,
                               annotations=contp$Significance[significant_indices],
                               step_increase = 0.1,
                               size = 1,
                               textsize = 6,
                               fontface = "bold",
                               vjust = .4)
  }

  if (!is.na(hline)) {
    plot <- plot + geom_hline(yintercept = hline, linetype = "dashed", color = "black")
  }
  if (!is.na(control_wellgroup)) {
    control_mean <- quick_stats %>%
      filter(wellgroup == control_wellgroup) %>%
      pull(meas_mean)
    plot <- plot + geom_hline(yintercept = control_mean, linetype = "dashed", color = "black")
  }

  # Create data_summary dataframe by grouping by wellgroup
  data_summary <- data_filtered %>%
    group_by(wellgroup) %>%
    summarise(nWells = n_distinct(well), nCells = n_distinct(cell)) %>%
    ungroup()
  
  # Calculate total nWells and nCells
  total_row <- data_summary %>%
    summarise(wellgroup = "Total", nWells = sum(nWells), nCells = sum(nCells))
  
  data_summary <- bind_rows(data_summary, total_row)   # Append the total row to data_summary


  #assemble stats results
  stats_details <- list(skewtest.raw = skewtest.raw,
                        qplot.raw = qplot.raw,
                        skewtest.transf = skewtest.transf,
                        qplot.transf = qplot.transf,
                        qplot.residuals = qplot.residuals,
                        gee_res = gee_res,
                        contmat = contmat,
                        gee_contrasts = gee_contrasts,
                        contp = contp,
                        p_matrix = p_matrix,
                        data_summary = data_summary
  )

  #create and name output
  GEE_and_plot <- list(plot = plot,
                       stats_details = stats_details
  )
  names(GEE_and_plot) <- c(paste0(period,"_",measurement_type,"_plot"),
                           paste0(period,"_",measurement_type,"stats_details")
  )
  return(GEE_and_plot)
}

create_stats_summary <- function(output_dir, expt_name, all_stats_details) {
  
  # Define the output file path
  output_file <- paste0(output_dir, expt_name, "_p_matrices_with_summary.txt")
  
  # Open connection to the file
  con <- file(output_file, open = "wt")
  
  # Iterate over all elements in all_stats_details
  for (test_name in names(all_stats_details)) {
    
    # Write the test name as a header
    writeLines(paste0("Test: ", test_name), con)
    
    # Extract the data_summary from each test's details
    data_summary <- all_stats_details[[test_name]]$data_summary
    
    # Write the data_summary if it exists
    if (!is.null(data_summary)) {
      writeLines("\nData Summary:", con)
      writeLines(capture.output(print(data_summary)), con)
    } else {
      writeLines("\nNo data_summary available for this test.", con)
    }
    
    # Extract the p_matrix from each test's details
    p_matrix <- all_stats_details[[test_name]]$p_matrix
    
    # Write the p_matrix if it exists
    if (!is.null(p_matrix)) {
      writeLines("\nP_matrix:", con)
      writeLines(capture.output(print(p_matrix)), con)
    } else {
      writeLines("\nNo p_matrix available for this test.", con)
    }
    
    # Extract the contp (contrast results) from each test's details
    contp <- all_stats_details[[test_name]]$contp
    
    # Write the contp if it exists
    if (!is.null(contp)) {
      writeLines("\nContrast Results (contp):", con)
      writeLines(capture.output(print(contp)), con)
    } else {
      writeLines("\nNo contrast results (contp) available for this test.", con)
    }
    
    # Add a separator for clarity
    writeLines("\n--------------------------------\n", con)
  }
  
  # Close the connection
  close(con)
  
  print(paste("p_matrices with summaries and contrast results saved to:", output_file))
}
locate_file_in_dir <- function(directories, pattern) {
  file_list <- c()
  for (directory in directories) {
    # Get the full path of files that match the pattern in each directory
    files <- list.files(directory, pattern = pattern, full.names = TRUE)
    file_list <- c(file_list, files)
  }
  return(file_list)
}

#################
# Open files, normalize data, compile, perform period measurements--------------------------------------------------------------------

#set up file manipulation
setwd(data_dir)
timestamp <- format(Sys.time(), "%y%m%d-%H%M%S")
output_dir <- paste0(dirname(data_dir), "/Analysis/e", expt_name, "_t", timestamp, "/")
if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

# Initialize lists for collecting plots, data, and long_data
csv_files <- list.files(pattern = "*.csv")
wellname_list <-list()
plot_list <- list()
data_list <- list()
norm_data_list <- list()
alldata_long <- data.frame()
removed_cells_tally <- 0

#compile data well by well
for(file in csv_files) {
  data <- read.csv(file)
  print(paste0(file, " opened"))
  print(dim(data))
  wellname <- sub("Results.csv$", "", file)
  wellname_list <- c(wellname_list, wellname)
  colnames(data)[1] <- "Time"   #rename columns
  colnames(data)[2] <- "title"
  colnames(data)[3] <- "bottomBG"
  data_list[[wellname]] = data #add raw data to list
  norm_data <- subset(data, select = -c(title))   #Remove labels
  norm_data[, c(3:ncol(norm_data))] <- norm_data[, c(3:ncol(norm_data))] - norm_data$bottomBG   # Background Subtraction: subtract BG from the relevant columns
  norm_data <- norm_data[,c(-2,-3)]   #remove background and slice label columns
  new_names <- c(names(norm_data)[1],
                 gsub("Mean\\.ROI(\\d+)\\.", paste0(wellname, "_\\1"), names(norm_data)[-1]))   #reformat ROI names
  names(norm_data) <- make.names(new_names, unique = TRUE)  # Ensure column names are syntactically valid and unique
  norm_data[,-1] <- data.frame(lapply(norm_data[,-1], normalize_column))   #normalize data to baseline and maximal stimulus

  #Identify bad traces and remove from dataset
  lowCells <- sapply(norm_data[-1], function(column) any(column < tooLow))
  highCells <- sapply(norm_data[-1], function(column) any(column[1:max(time_period_defs$EndTime)] > tooHigh))
  baseHighCells <- sapply(norm_data[-1], function(column) any(column[1:30] > baseTooHigh))
  cells_to_remove <- lowCells | highCells | baseHighCells
  removed_cells_tally <- removed_cells_tally + length(which(cells_to_remove == TRUE))
  norm_data <- data.frame(Time = norm_data$Time, norm_data[-1][!cells_to_remove])

  norm_data_list[[wellname]] = norm_data  #add normalized data to list

  #convert to longform
  long_data <- norm_data %>%
    pivot_longer(cols = -Time, names_to = "cell", values_to = "intensity") %>%
    mutate(well = wellname)
  
  alldata_long <- rbind(alldata_long, long_data)   #compile all data


}

print(paste0(removed_cells_tally, " bad traces removed"))

#add metadata to compiled data and then reorder columns
alldata_long <- alldata_long %>%
  mutate(expt = expt_name) %>%
  mutate(wellgroup = find_wellgroup(well)) %>%
  mutate(well = paste0(expt, "_", well)) %>%
  mutate(cell = paste0(expt, "_", cell))
alldata_long <- alldata_long[, c("Time", "intensity", "cell", "well", "expt","wellgroup")]

write.csv(alldata_long, paste0(output_dir, expt_name, "_normData.csv"), row.names = FALSE) #Save compiled data

all_cell_measurements <- measure_cells_by_period(alldata_long, time_period_defs) # Measure cells by time period
write.csv(all_cell_measurements, paste0(output_dir, expt_name, "_measurements.csv"), row.names = FALSE) #Save measurements

##########################
##Add datasets together###############

# Specify the output folder for saving combined results
top_output_dir <- "/Volumes/T7_Shield/Imaging (monolayer)/Multi_experiment_analysis"
expt_name <- "Fig_s2_SST_YM_base_test"

# Specify the list of input directories
input_dirs <- c(
  '/Volumes/T7_Shield/Imaging (monolayer)/241219/Analysis/e2401219_t241223-145927',
  '/Volumes/T7_Shield/Imaging (monolayer)/241220/Analysis/e2401220_t241223-160120'
)

wellgroup_names <- c(
  "Veh",
  "SST",
  "YM - Veh",
  "YM - SST"
)

wellgroups_to_analyze <- wellgroup_names

# Trace plot parameters
y_limits = c(-.2, .7)
x_limits = c(20, 80) # NA is full range

# Stats plot parameters
periods_to_analyze <- c("P3")
set.stats.alpha <- 0.5
set.cex <- 0.8
hide_comparisons <- "FALSE"
stats_plot_width <- 0.8 * length(wellgroups_to_analyze)
stats_plot_height <- 3

# Locate and compile timeseries data
alldata_long <- data.frame()
timedata_list <- locate_file_in_dir(input_dirs, "_normData.csv$")
for (file in timedata_list) {
  tData <- read.csv(file)
  alldata_long <- rbind(alldata_long, tData)
}

# Locate and compile period measurement data
all_cell_measurements <- data.frame()
measurements_list <- locate_file_in_dir(input_dirs, "_measurements.csv$")
for (file in measurements_list) {
  meas <- read.csv(file)
  all_cell_measurements <- rbind(all_cell_measurements, meas)
}

#Fully name output directory and create if absent
timestamp <- format(Sys.time(), "%y%m%d-%H%M%S")
output_dir <- paste0(top_output_dir, "/", expt_name, "/e", expt_name, "_t", timestamp, "/")
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}

# Save the combined data
write.csv(alldata_long, paste0(output_dir, expt_name, "_combined_normData.csv"), row.names = FALSE)
write.csv(all_cell_measurements, paste0(output_dir, expt_name, "_combined_measurements.csv"), row.names = FALSE)

###########################
#####Create TracePlots############################
# Create Group Trace Grid -------------------------------------------------------------
wellgroup_plot_list <- list()

for (wellgroup in wellgroup_names) {
  plotname <- paste(wellgroup, collapse = ", ")
  data_to_plot <- alldata_long[alldata_long$wellgroup == wellgroup, ]
 
  trace_plot <- createTracePlot(
    data = data_to_plot, 
    title = plotname, 
    y_limits = y_limits, 
    x_limits = x_limits, 
    line_x_intercepts = drug_times, 
    text_scale = 12
  )

  #add plot to list
  wellgroup_plot_list[[wellgroup]] = trace_plot
  
}

#generate grid
wellgroup_plot_grid <- grid.arrange(grobs = wellgroup_plot_list, ncol = length(wellgroup_plot_list))
grid.draw(wellgroup_plot_grid)

#save traceplotgrid
total_width <- plot_width * length(wellgroup_names)
total_height <- plot_height 
output_filename <- paste0(output_dir, expt_name, "_grouptrace.pdf")
ggsave(output_filename, wellgroup_plot_grid, width = total_width, device = "pdf", height = total_height, limitsize = FALSE, dpi = 300)


# Create Well Trace Grid -------------------------------------------------------------
plot_list <- NULL
if (exists("wellname_list") != TRUE) {
  print("wellname_list not defined, Well Trace Grid not generated") 
} else {
  plot_list <- NULL
  for (wellname in unique(alldata_long$well)){
    plotname <- paste(wellname, collapse = ", ")
    data_to_plot <- alldata_long[alldata_long$well %in% wellname, ]
    
    trace_plot <- createTracePlot(
      data = data_to_plot, 
      title = plotname, 
      y_limits = y_limits, 
      x_limits = x_limits, 
      line_x_intercepts = drug_times, 
      text_scale = 12
    )
    
    short_wellname <- sub(".*_(.*)", "\\1", wellname)
    
    #add plot to list
    plot_list[[short_wellname]] = trace_plot
    
  }
  
  #Create GRID
  #Get Rows and Columns
  ordered_letters <- sort(unique(substr(wellname_list, 1, 1)))
  ordered_numbers <- sort(unique(substr(wellname_list, 2, 3)), decreasing = FALSE)
  
  # Create an empty list for the full grid
  ordered_plots <- vector("list", length=length(ordered_letters)*length(ordered_numbers))
  names(ordered_plots) <- expand.grid(ordered_numbers, ordered_letters) %>% 
    rowwise() %>% 
    mutate(comb = paste0(Var2, Var1)) %>% 
    pull(comb)
  
  # Place each plot in its rightful location in the grid
  for (wellname in wellname_list) {
    ordered_plots[[wellname]] <- plot_list[[wellname]]
    print(paste0(wellname, " added to grid"))
  }
  
  # Fill in empty plots with a blank ggplot
  empty_plot <- ggplot() + theme_void() + theme(plot.background = element_blank())
  grid.draw(empty_plot)
  ordered_plots[sapply(ordered_plots, is.null)] <- list(empty_plot)
  combined_plot <- arrangeGrob(grobs = ordered_plots, ncol = length(ordered_numbers))
  
  #save traceplotgrid
  total_width <- plot_width * length(ordered_numbers)
  total_height <- plot_height * length(ordered_letters)
  output_filename <- paste0(output_dir, expt_name, "_welltrace.pdf")
  ggsave(output_filename, combined_plot, width = total_width, height = total_height, limitsize = FALSE, dpi = 300)
}




########################
# Wellgroup GEE and Plot ---------------------------------------------------------

data_to_plot <- all_cell_measurements %>%
  filter(!is.na(all_cell_measurements$wellgroup))

all_stats_plots <- list()
all_stats_details <- list()

#* -- Dif Max Plot -----
metric <- "dif_max"
for (period in periods_to_analyze) {
  results <- run_GEE_and_plot(data = data_to_plot,
                              measurement_type = metric,
                              transform_data = FALSE,
                              FDR = TRUE,
                              y_label = "Dif. Max dF/F0 (treat.-base.)",
                              period = period,
                              label_height = 0.1,
                              hline = NA
  )
  name <- paste0(period,"_",metric)
  all_stats_plots[name] <- results[1]
  all_stats_details[name] <- results[2]
}

# #* -- Max Plot -----
# metric <- "max"
# for (period in periods_to_analyze) {
#   results <- run_GEE_and_plot(data = data_to_plot,
#                               measurement_type = metric,
#                               transform_data = TRUE,
#                               FDR = TRUE,
#                               y_label = "Max Norm. dF/F0",
#                               period = period,
#                               label_height = 0.1,
#                               hline = NA
#   )
#   name <- paste0(period,"_",metric)
#   all_stats_plots[name] <- results[1]
#   all_stats_details[name] <- results[2]
# }
# 
#* -- Rel Volatility Plot -----
metric <- "ln_rel_volatility"
for (period in periods_to_analyze) {
  results <- run_GEE_and_plot(data = data_to_plot,
                              measurement_type = metric,
                              transform_data = FALSE,
                              FDR = TRUE,
                              y_label = "Volatility (ln[treat./base.])",
                              period = period,
                              label_height = 0.1,
                              hline = NA
  )
  name <- paste0(period,"_",metric)
  all_stats_plots[name] <- results[1]
  all_stats_details[name] <- results[2]
}

stats_grid <- grid.arrange(grobs = all_stats_plots, ncol = length(all_stats_plots))
grid.draw(stats_grid)
# #* -- Abs Volatility Plot -----
# metric <- "volatility"
# for (period in periods_to_analyze) {
#   results <- run_GEE_and_plot(data = data_to_plot,
#                               measurement_type = metric,
#                               transform_data = TRUE,
#                               FDR = TRUE,
#                               y_label = "Volatility",
#                               period = period,
#                               label_height = 0.1,
#                               hline = NA
#   )
#   name <- paste0(period,"_",metric)
#   all_stats_plots[name] <- results[1]
#   all_stats_details[name] <- results[2]
# }
# 
# stats_grid <- grid.arrange(grobs = all_stats_plots, ncol = length(all_stats_plots))
# grid.draw(stats_grid)
# 
#* -- Save Stats Results and Plots --------------------------------------------------------------------
#save plots
stats_total_width <- stats_plot_width * length(all_stats_plots)
stats_total_height <- stats_plot_height
output_filename <- paste0(output_dir, expt_name, "_groupstats.pdf")
ggsave(output_filename, stats_grid, device = "pdf", width = stats_total_width, height = stats_total_height, limitsize = FALSE, dpi = 300)

# Save stats results
output_filename <- paste0(output_dir, expt_name, "_all_stats_details.RData")
save(all_stats_details, file = output_filename)
create_stats_summary(output_dir, expt_name, all_stats_details)

########################
# Identify "Responders" ---------------------------------------------------
# Filter measurements for pre, stim, and sustained periods
pre_measurements <- all_cell_measurements %>% 
  filter(tperiod.name == "P1") %>% select(cell, well, max) %>% rename(max_pre = max)
stim_measurements <- all_cell_measurements %>% 
  filter(tperiod.name == "P2") %>% select(cell, well, max) %>% rename(max_stim = max)
sust_measurements <- all_cell_measurements %>%
  filter(tperiod.name == "P3") %>% select(cell, well, max) %>% rename(max_sust = max)

# Merge pre, stim, and sust measurements by cell and well
pre_stim_sust_meas <- pre_measurements %>%
  merge(stim_measurements, by = c("cell", "well")) %>%
  merge(sust_measurements, by = c("cell", "well"))

# Filter to find cells where max_stim is greater than max_pre
responder_meas <- pre_stim_sust_meas %>%
  filter(max_stim / max_pre > 2) %>%
  filter(max_stim - max_pre > 0.09) %>%
  filter(max_stim - max_sust > 0.07) 
  
responder_list <- responder_meas$cell # List of cells with higher max in stim period

# Add 'class' columns to data and meas
alldata_long <- alldata_long %>%
  mutate(class = ifelse(cell %in% responder_list, "responder", "non"))

all_cell_measurements <- all_cell_measurements %>%
  mutate(class = ifelse(cell %in% responder_list, "responder", "non"))

# Create Responder Group Trace Grid -------------------------------------------------------------
resp_wellgroup_plot_list <- list()
total_count_vector <- c()
resp_count_vector <- c()
percent_resp_vector <-c()
for (wellgroup_name in wellgroup_names) {
  data_to_plot <- alldata_long %>%
    filter(wellgroup == wellgroup_name)
  
  nCells <- length(unique(data_to_plot$cell))
  total_count_vector <- c(total_count_vector, nCells)
  data_to_plot <- data_to_plot[data_to_plot$class == "responder", ]
  
  nResponders <- length(unique(data_to_plot$cell))
  resp_count_vector <- c(resp_count_vector, nResponders)
  
  percent_resp <- if (nCells > 0) round(nResponders / nCells * 100, 1) else 0
  percent_resp_vector <- c(percent_resp_vector, percent_resp)
  
  if (any(data_to_plot$class == "responder") == TRUE){
    
    plotname <- paste(wellgroup_name, collapse = ", ")
    trace_plot <- createTracePlot(
      data = data_to_plot, 
      title = plotname, 
      y_limits = c(-0.3, 0.7), 
      x_limits = c(0, 80), 
      line_x_intercepts = drug_times, 
      text_scale = 12
    )
    
    # Remove nCells annotation and replace
    x_midpoint <- mean(range(data_to_plot$Time, na.rm = TRUE))
    if (length(trace_plot$layers) >= 3) {
      trace_plot$layers <- trace_plot$layers[-3] # Assuming the 3rd layer is the nCells annotation
    }
    trace_plot <- trace_plot +
      annotate("text", x = x_midpoint, y = y_limits[2] - 0.05,
               label = paste0("n =", nResponders, "/", nCells, " (", percent_resp, "%)"),
               size = 3,
               colour = "black")
    
    # Add the plot to the list
    resp_wellgroup_plot_list[[wellgroup_name]] <- trace_plot
    
  } else {
    # Print message to indicate no responders for this wellgroup
    print(paste("No responders found for wellgroup:", wellgroup_name, "- Skipping plot"))
  }
}

# Generate grid only if there are plots to display
if (length(resp_wellgroup_plot_list) > 0) {
  resp_wellgroup_plot_grid <- grid.arrange(grobs = resp_wellgroup_plot_list, ncol = length(resp_wellgroup_plot_list))
  grid.draw(resp_wellgroup_plot_grid)
  
  # Save traceplot grid
  total_width <- plot_width * length(wellgroup_names)
  total_height <- plot_height
  output_filename <- paste0(output_dir, expt_name, "_responder_grouptrace.pdf")
  ggsave(output_filename, resp_wellgroup_plot_grid, width = total_width, device = "pdf", height = total_height, limitsize = FALSE, dpi = 300)
} else {
  print("No responder plots to generate.")
}


# Fishers Exact Test of Responders And Bar Plot----------------------------------------
responder_test_results <- pairwise.fisher.test(resp_count_vector, total_count_vector, p.adjust.method = "bonferroni")
colnames(responder_test_results[["p.value"]]) <- wellgroup_names[1:(length(wellgroup_names) - 1)]
rownames(responder_test_results[["p.value"]]) <- wellgroup_names[2:length(wellgroup_names)]

# Create the data frame with wellgroup as a factor
resp_summary_df <- data.frame(
  wellgroup = factor(unlist(wellgroup_names), levels = unlist(wellgroup_names)),
  percent_responders = percent_resp_vector,
  responders = resp_count_vector,
  total_cells = total_count_vector
)


# Create the bar plot
bar_plot <- ggplot(resp_summary_df, aes(x = wellgroup, y = percent_responders, fill = wellgroup)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_text(
    aes(label = paste0(responders, "/", total_cells)),
    vjust = -0.5, size = 3.5
  ) +
  labs(
    title = "Percentage of Responders per Group",
    x = "Wellgroup",
    y = "Percentage of Responders"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 11, color = "black", angle = 20, hjust = 1),  # Set x-axis tick label color to black
    axis.text.y  = element_text(size = 12, color = "black"),  # Set y-axis tick label color to black
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(),
    axis.line=element_line(),
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    expand = expansion(mult = c(0, 0.1))
  )
    
# Display and save the plot
print(bar_plot)
output_filename <- paste0(output_dir, expt_name, "_percent_responders_barplot.pdf")
ggsave(output_filename, bar_plot, width = 3.5, height = 3.5, device = "pdf", limitsize = FALSE, dpi = 300)

output_filename <- paste0(output_dir, expt_name, "_percent_responders_summary.csv")
write.csv(resp_summary_df, file = output_filename, row.names = TRUE)

output_filename <- paste0(output_dir, expt_name, "_percent_responders_test.csv")
write.csv(as.data.frame(responder_test_results$p.value), file = output_filename, row.names = TRUE)


