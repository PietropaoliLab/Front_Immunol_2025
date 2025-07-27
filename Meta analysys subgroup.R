library(here)
library(meta)
library(readxl)
library(ggplot2)
library(dplyr)
library(readr)
library(grid)
library(kableExtra)
library(patchwork)

# Define results folder
ResultFolder <- here("AQ Projects/IL1 Meta-Analysis/Results/")

# https://www.metafor-project.org/doku.php/analyses:morris2008
here()
Database_MetaAnalisys <- read_excel(here("AQ Projects/IL1 Meta-Analysis/Dataset.xlsx"))

# str(Database_MetaAnalisys)

df <- as.data.frame(Database_MetaAnalisys)


# df <- df[c(1:12,22)]

# Calculate effect sizes (mean difference)
df <- df %>%
  mutate(Mean_Diff = POST_mean - PRE_mean,
         SE = sqrt((PRE_SD^2 / Subjects) + (POST_SD^2 / Subjects)))


df$Smoke_binary <- ifelse(df$Smokers > 0, "Sample with smokers", "Sample without smokers")
df$Diabetes_binary <- ifelse(df$With_Diabetes > 0, "Sample with diabetic", "Sample without diabetic")
df$Female_binary <- ifelse(df$FemaleRatio == 100, "Only Female",
                           ifelse(df$FemaleRatio == 0, "Only Male", "Sex mixed"))

df$MonthsAfterTRT <- cut(df$Elapsed_time_days, breaks = c(0, 42, 90, 180, Inf), labels = c("4-6Wks", "Within 3 Mo.", "Within 6 Mo.", "6Mo. and more"))

# Inizialize dataset for storing all the analysis
results <- data.frame()

# Define the function
create_forest_plot <- function(cytokine_name, treatment_type, scale_factor, x_min, x_max, Subgroup, Label, Moderators) {
  
  # Define results folder and PDF file path
  ResultFolder <- here("AQ Projects/IL1 Meta-Analysis/Results/")
  pdf_file <- file.path(ResultFolder, paste0(cytokine_name, "_", treatment_type ,"_results.pdf"))
  
  # Filter data for the specific cytokine and treatment type
  df_filtered <- df %>% filter(Cytokine == cytokine_name, Type_of_Treatment == treatment_type)
  
  # Debugging: Check if df_filtered is properly defined
  if (is.null(df_filtered) || nrow(df_filtered) == 0) {
    stop(paste("df_filtered is not properly defined or has no rows. Cytokine:", cytokine_name, "- Treatment:", treatment_type))
  } else {
    print(paste("df_filtered has", nrow(df_filtered), "rows"))
  }
  
  # Check if there's sufficient data for meta-analysis
  if (nrow(df_filtered) < 2) {
    message(paste("Not enough data for cytokine:", cytokine_name, "- Treatment:", treatment_type))
    return(NULL)
  }
  
  # Rescale SMD and SE values
  df_filtered$Mean_Diff <- df_filtered$Mean_Diff * scale_factor
  df_filtered$SE <- df_filtered$SE * scale_factor
  
  # Open a PDF device to save the plots
  pdf(pdf_file, width = 16, height = 12)  # Adjust width and height as needed
  
  csv_destination <- file.path(ResultFolder, paste0("Results.csv"))
  csv_meta <- file.path(ResultFolder, paste0("MetaAnalysis_results.csv"))
  csv_SubgroupMeta <- file.path(ResultFolder, paste0("Subgroup_MetaAnalysis_results.csv"))
  
  results <- data.frame()
  MetaAnalisysEffect <- data.frame()
  SubgroupEffect <- data.frame()
  
  if (file.exists(csv_destination)) {
    message("Storing CSV for Overall results present")
  } else {
    write_csv(results, file = csv_destination)
  }
  
  
  if (file.exists(csv_meta)) {
    message("Storing CSV for MetaAnalysis effect present")
  } else {
    write_csv(MetaAnalisysEffect, file = csv_meta)
  }
  
  if (file.exists(csv_SubgroupMeta)) {
    message("Storing CSV for Subgroups effect present")
  } else {
    write_csv(SubgroupEffect, file = csv_SubgroupMeta)
  }

  k = 1
  for (i in Subgroup) {
    
    # Debugging: Print subgroup being processed
    print(paste("Processing subgroup:", i))
    
    df_filteredNA <- subset(df_filtered, !is.na(df_filtered[[i]]))
    
    # Check if there's sufficient data for meta-analysis
    if (nrow(df_filteredNA) < 2) {
      message(paste("Not enough data for subgroup analysis:", i))
      next
    }
    
    # Run meta-analysis with/without Hartung-Knapp adjustment to handle large variances
    meta_analysis <- metacont(n.e = Subjects,  # Sample size of experimental group
                              mean.e = POST_mean,  # Post-treatment mean (experimental group)
                              sd.e = POST_SD,  # Post-treatment SD
                              n.c = Subjects,  # Sample size of control group (same as experimental in paired design)
                              mean.c = PRE_mean,  # Pre-treatment mean (control group)
                              sd.c = PRE_SD,  # Pre-treatment SD
                              studlab = Study_ID,  # Study labels
                              data = df_filteredNA,  # Your dataframe
                              sm = "SMD",  # Mean Difference
                              hakn = FALSE,  # Hartung-Knapp adjustment
                              comb.random = TRUE,  # Random effects model
                              byvar = as.factor(df_filteredNA[[i]]))  # Subgroup analysis based on elapsed time
    
    # Adjust forest labels/indications
    meta_analysis[["subgroup.name"]] <- Label[k]
    meta_analysis[["sep.subgroup"]] <- ":"
    meta_analysis[["label.e"]] <- "After TRT"
    meta_analysis[["label.c"]] <- "Before TRT"
    
    # Adjust values for visualization 
    meta_analysis[["mean.e"]] <- meta_analysis[["mean.e"]] * scale_factor
    meta_analysis[["sd.e"]] <- meta_analysis[["sd.e"]] * scale_factor
    meta_analysis[["mean.c"]] <- meta_analysis[["mean.c"]] * scale_factor
    meta_analysis[["sd.c"]] <- meta_analysis[["sd.c"]] * scale_factor
    
    # Create forest plot with adjusted xlim and scaled axis label
    forest(meta_analysis,
           overall = TRUE,           # Show overall effect
           studlab = TRUE,            # Show study labels
           comb.fixed = FALSE,        # Use random-effects model
           xlab = paste0("\nSMD (mg/L) - ", "Scale factor:", scale_factor),
           col.diamond = "#c85200",      # Diamond color for overall effect
           col.square = "#6b8ea4",     # Square color for individual studies
           col.diamond.lines = "black", # Diamond outline color
           col.subgroup = "#c85200",# # Labels colors for subgroups
           fontsize = 7, spacing = 0.5, plotwidth = "80mm", 
           xlim = c(x_min, x_max))  # Adjust the x-axis limits
    
    # Add title using grid.text
    grid.text(paste("Cytokine:", cytokine_name, "- Treatment:", treatment_type),
              x = 0.5, y = 0.95, gp = gpar(cex = 1))
    
    # Extract overall p-value from meta-analysis summary
    meta_summary <- summary(meta_analysis)
    
    # Print overall p-value
    if (!is.null(meta_summary$pval.random)) {
      p_overall <- meta_summary$pval.random
      e_overall <- paste0("SMD: ", round(meta_summary$TE.random, 2), "; 95%C.I.[",round(meta_summary$lower.random, 2), "; ",round(meta_summary$upper.random, 2), "]")
      grid.text(paste(e_overall, "P value for random effect=", round(p_overall, 4)), 
                x = 0.5, y = 0.10, gp = gpar(cex = 0.7))
    }
    
    # Print common p-value
    if (!is.null(meta_summary$pval.common)) {
      p_common <- meta_summary$pval.common
      e_common <- paste0("SMD: ", round(meta_summary$TE.common, 2), "; 95%C.I.[",round(meta_summary$lower.common, 2), "; ",round(meta_summary$upper.common, 2), "]")
      grid.text(paste(e_common, "P value for common effect=", round(p_common, 4)), 
                          x = 0.5, y = 0.12, gp = gpar(cex = 0.7))
    }
    
    # Store Meta Analysys effect for further use
    MetaAnalisysEffect <- data.frame("Cytokine" = cytokine_name,
                                     "Treatment"= treatment_type,
                                     "Subgroup" = Subgroup[k],
                                     "NumberOfStudis" = length(meta_summary$studlab),
                                     "NumberOfSubjects" = meta_analysis$n.e.w,
                                     "SMD_random" = round(meta_summary$TE.random, 2),
                                     "seSMD_random" = round(meta_summary$seTE.random.w, 2),
                                     "SMD_random_lower" = round(meta_summary$lower.random, 2),
                                     "SMD_random_upper" = round(meta_summary$upper.random, 2),
                                     "Pvalue_random" = round(meta_summary$pval.random, 4),
                                     "SMD_common" = round(meta_summary$TE.common, 2),
                                     "seSMD_common" = round(meta_summary$seTE.common.w, 2),
                                     "SMD_common_lower" = round(meta_summary$lower.common, 2),
                                     "SMD_common_upper" = round(meta_summary$upper.common, 2),
                                     "Pvalue_common" = round(meta_summary$pval.common, 4))
    
    # Verify if csv file is present
    if ( nrow(data.frame(read_csv(csv_meta, show_col_types = FALSE))) == 0 ) {
      write_csv(MetaAnalisysEffect, file = csv_meta, col_names = TRUE)
    } else {
      write_csv(MetaAnalisysEffect, file = csv_meta, append = TRUE)
    }
    
    
   
    SubgroupEffect <- data.frame( "Cytokine_name" = cytokine_name,
                                  "Treatment_type" = treatment_type,
                                  "Subgroup" = Subgroup[k],
                                  "Subgroup_levels" = meta_analysis$subgroup.name,
                                  "NumberOfStudies" = meta_analysis$k.w,
                                  "NumberOfSubjects" = meta_analysis$n.e.w,
                                  "TE_subgroup_common" = meta_analysis$TE.common.w,
                                  "Subgroup_lower_common" = meta_analysis$lower.common.w,
                                  "Subgroup_upper_common"= meta_analysis$upper.common.w,
                                  "Subgroup_pvalue_common"= meta_analysis$pval.common.w,
                                  "TE_subgroup_random" = meta_analysis$TE.random.w,
                                  "Subgroup_lower_random" = meta_analysis$lower.random.w,
                                  "Subgroup_upper_random"= meta_analysis$upper.random.w,
                                  "Subgroup_pvalue_random"= meta_analysis$pval.random.w)
    SubgroupEffect$Subgroup_levels <- rownames(SubgroupEffect)
    
    
    # Verify if csv file is present
    if ( nrow(data.frame(read_csv(csv_SubgroupMeta, show_col_types = FALSE))) == 0 ) {
      write_csv(SubgroupEffect, file = csv_SubgroupMeta, col_names = TRUE)
    } else {
      write_csv(SubgroupEffect, file = csv_SubgroupMeta, append = TRUE)
    }

   
     # Save meta-analysis results as HTML table
    html_destination <- file.path(ResultFolder, paste0("HTML/",cytokine_name, "_", treatment_type, "_", Subgroup[k], "_results.html"))
    billy <- as.data.frame(meta_analysis)[1:18]
    
    # Store all data into a single dataset
    tmp <- billy
    tmp$cytokine_name <- cytokine_name
    tmp$treatment_type <- treatment_type
    tmp$subgrouping_names <- Subgroup[k]
    results <- rbind(results, tmp)
    
    # Verify if csv file is present
    if ( nrow(data.frame(read_csv(csv_destination, show_col_types = FALSE))) == 0 ) {
    write_csv(results, file = csv_destination, col_names = TRUE)
    } else {
      write_csv(results, file = csv_destination, append = TRUE)
    }
    
    # Generate HTML results file
    kable(billy, format = "html", 
          caption = paste0('<div style="text-align: center; color: black;">',
                           cytokine_name, " - ", treatment_type, 
                           " treatment - Subgrouping: ", Subgroup[k], 
                           '</div>'),          
          col.names = c("N", "Mean", "SD", "N", "Mean", "SD", "Study", 
                        "TE", "seTE", "Statistic", "P value", "Lower limit", "Upper limit", "W common", 
                        "W random", "Z value", "W fixed", "Subgroup")) %>% 
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 9, full_width = FALSE) %>%
      add_header_above(c("After TRT" = 3, "Before TRT" = 3, "Data" = 12)) %>%
      add_footnote(c(paste0("Scaling factor:", scale_factor)), notation = "number") %>%
      save_kable(file = html_destination, self_contained = TRUE)
    
    k = k + 1
    
    rm(tmp, billy)
  }
  
  
  ##################################
  #
  #  Begin meta-regression section
  #
  ##################################
  # List to store plots
  p <- list()
  k = 1
  
  # Loop through each moderator for meta-regression
  for (m in Moderators) {
    # Debugging: Print the moderator being processed
    print(paste("Processing moderator:", m))
    
    # Subset the data to remove rows with NA in the moderator column
    df_filtered_non_na <- df %>% filter(Cytokine == cytokine_name, 
                                        Type_of_Treatment == treatment_type,
                                        !is.na(.data[[m]]))
    
    # Check if the filtered data is valid (non-empty)
    if (nrow(df_filtered_non_na) < 2) {
      message(paste("Not enough data for meta-regression with moderator:", m))
      next  # Skip to the next moderator if there is insufficient data
    }
    
    # Debugging: Check if the subset is created correctly
    print(paste("df_filtered_non_na has", nrow(df_filtered_non_na), "rows"))
    
    # Perform meta-analysis again for meta-regression on the filtered data
    meta_analysis <- metacont(n.e = Subjects,  # Sample size of experimental group
                              mean.e = POST_mean,  # Post-treatment mean (experimental group)
                              sd.e = POST_SD,  # Post-treatment SD
                              n.c = Subjects,  # Sample size of control group (same as experimental in paired design)
                              mean.c = PRE_mean,  # Pre-treatment mean (control group)
                              sd.c = PRE_SD,  # Pre-treatment SD
                              studlab = Study_ID,  # Study labels
                              data = df_filtered_non_na,  # Filtered dataframe
                              sm = "SMD",  # Mean Difference
                              hakn = TRUE,  # Hartung-Knapp adjustment
                              comb.random = TRUE)
    message("Results of metacont:\n")
    print(meta_analysis)
    
    # Perform meta-regression with the specified moderator
    meta_regression <- metareg(meta_analysis, as.formula(paste("~", m)))
    
    message("Results of metareg:\n")
    print(meta_regression)
    
    # Extract meta-regression results
    summary(meta_regression)
    
    # Prepare data for ggplot
    df_filtered_non_na$SMD <- meta_analysis$TE  # Extract the effect sizes (SMD)
    df_filtered_non_na$SMD_se <- meta_analysis$seTE  # Extract standard errors of SMD
    
    # MyTheme
    MyTheme <- theme_linedraw(10)+ 
      theme(axis.text = element_text(colour = "black"),
            panel.grid = element_blank(),
            legend.key = element_blank(),
            legend.background=element_blank(),
            aspect.ratio = 1,
            strip.background = element_rect(fill = NA, colour = NA),
            text = element_text(family = "Helvetica"))
    
    # Plot using ggplot2
    p[[k]] <- ggplot(df_filtered_non_na, aes(x = .data[[m]], y = SMD)) +
      geom_point(aes(size = 1 / SMD_se), alpha = 0.5) +  # Scatter plot with point size inversely proportional to SE
      scale_size_continuous(range = c(2, 4)) +  # Adjust point size range (lower numbers for smaller points)
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#366785") +  # Regression line
      labs(x = paste0(m),
           y = "Standardized Mean Difference (SMD)",
           caption = paste("Meta-Regression: SMD ~ ", m, "\n",
                           "Slope: ", round(meta_regression$beta[2], 3),
                           "\nP-value: ", round(meta_regression$pval[2], 3),
                           "\nR²: ", round(meta_regression$R2, 3))) +
      MyTheme
    
    k = k + 1
  }
  
  # Print all the plot in a single row
  combined_plot <- wrap_plots(p, nrow = 2) + 
    plot_annotation(title = paste("Meta-Regression for SMD on ", cytokine_name, "- Treatment:", treatment_type))
  
  print(combined_plot)
 
   # Close the PDF device
  dev.off()
  
  # Message to indicate the PDF has been saved
  message(paste0("All plots have been saved to", pdf_file, "\n"))
  
  rm(df_filtered_non_na, df_filtered, df_filteredNA, meta_analysis, meta_summary,
     results, MetaAnalisysEffect, SubgroupEffect)
}



# TNF-a
create_forest_plot("TNF-a", "Standard", scale_factor = 1000000, x_min = -2, x_max = 2, 
                   Subgroup = c("Elapsed_time_days", "MonthsAfterTRT", "Assay_type", "Bias_judgment", "Study_design",  "Additional_devices","Female_binary", "Smoke_binary", "Diabetes_binary"),   # Subgroups for the meta analysis
                   Label = c("Days after treatment", "MonthsAfterTRT","Assay type", "Bias judgment", "Study design",   "Additional_devices","Female_binary", "Smoke_binary", "Diabetes_binary"),   # Labels for the subgroups
                   Moderators = c("Elapsed_time_days", "Mean_age", "Study_Yr", "PRE_mean_N_teeth", "With_Periodontitis", "With_Diabetes", "FemaleRatio","Smokers"))         # Meta-regression moderators

create_forest_plot("TNF-a", "Intensive", scale_factor = 1000000, x_min = -4, x_max = 4, 
                   Subgroup = c("Elapsed_time_days", "MonthsAfterTRT", "Assay_type", "Bias_judgment", "Study_design", "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Subgroups for the meta analysis
                   Label = c("Days after treatment", "MonthsAfterTRT", "Assay type", "Bias judgment", "Study design", "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Labels for the subgroups
                   Moderators = c("Elapsed_time_days","Mean_age", "Study_Yr", "PRE_mean_N_teeth", "With_Periodontitis", "With_Diabetes", "FemaleRatio","Smokers"))                             # Meta-regression moderators


# IL-1beta
create_forest_plot("IL-1beta", "Standard", scale_factor = 1000000, x_min = -.9, x_max = .9, 
                   Subgroup = c("Elapsed_time_days", "MonthsAfterTRT", "Assay_type", "Bias_judgment", "Study_design",  "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Subgroups for the meta analysis
                   Label = c("Days after treatment", "MonthsAfterTRT", "Assay type", "Bias judgment", "Study design",  "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Labels for the subgroups
                   Moderators = c("Elapsed_time_days","Mean_age", "Study_Yr", "PRE_mean_N_teeth","FemaleRatio"))                             # Meta-regression moderators
create_forest_plot("IL-1beta", "Intensive", scale_factor = 1000000, x_min = -10, x_max = 4, 
                   Subgroup = c("Elapsed_time_days", "MonthsAfterTRT", "Assay_type", "Bias_judgment", "Study_design", "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Subgroups for the meta analysis
                   Label = c("Days after treatment", "MonthsAfterTRT", "Assay type", "Bias judgment", "Study design", "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Labels for the subgroups
                   Moderators = c("Elapsed_time_days","Mean_age", "Study_Yr", "PRE_mean_N_teeth", "With_Periodontitis", "FemaleRatio", "With_Diabetes", "Smokers"))                             # Meta-regression moderators

# IL-6
create_forest_plot("IL-6", "Standard", scale_factor = 1000000, x_min = -2, x_max = 2, 
                   Subgroup = c("Elapsed_time_days", "MonthsAfterTRT", "Assay_type", "Bias_judgment", "Study_design", "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Subgroups for the meta analysis
                   Label = c("Days after treatment", "MonthsAfterTRT", "Assay type", "Bias judgment", "Study design", "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Labels for the subgroups
                   Moderators = c("Elapsed_time_days","Mean_age", "Study_Yr", "PRE_mean_N_teeth", "With_Periodontitis", "FemaleRatio","With_Diabetes", "Smokers"))                             # Meta-regression moderators
create_forest_plot("IL-6", "Intensive", scale_factor = 1000000, x_min = -3, x_max = 3, 
                   Subgroup = c("Elapsed_time_days", "MonthsAfterTRT", "Assay_type", "Bias_judgment", "Study_design", "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Subgroups for the meta analysis
                   Label = c("Days after treatment", "MonthsAfterTRT", "Assay type", "Bias judgment", "Study_design", "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Labels for the subgroups
                   Moderators = c("Elapsed_time_days","Mean_age", "Study_Yr", "PRE_mean_N_teeth", "With_Periodontitis", "FemaleRatio","With_Diabetes", "Smokers"))                             # Meta-regression moderators

# CRP
create_forest_plot("CRP", "Standard", scale_factor = 10, x_min = -2.2, x_max = 2.2, 
                   Subgroup = c("Elapsed_time_days", "MonthsAfterTRT", "Assay_type", "Bias_judgment", "Study_design", "Additional_devices","Female_binary", "Smoke_binary", "Diabetes_binary"),   # Subgroups for the meta analysis
                   Label = c("Days after treatment", "MonthsAfterTRT", "Assay type", "Bias judgment", "Study design", "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Labels for the subgroups
                   Moderators = c("Elapsed_time_days","Mean_age", "Study_Yr", "PRE_mean_N_teeth", "With_Periodontitis", "FemaleRatio","With_Diabetes", "Smokers"))                             # Meta-regression moderators
create_forest_plot("CRP", "Intensive", scale_factor = 10, x_min = -3, x_max = 3, 
                   Subgroup = c("Elapsed_time_days", "MonthsAfterTRT", "Assay_type", "Bias_judgment", "Study_design", "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Subgroups for the meta analysis
                   Label = c("Days after treatment", "MonthsAfterTRT", "Assay type", "Bias judgment", "Study design", "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Labels for the subgroups
                   Moderators = c("Elapsed_time_days","Mean_age", "Study_Yr", "PRE_mean_N_teeth", "With_Periodontitis", "FemaleRatio","With_Diabetes", "Smokers"))                             # Meta-regression moderators

# hs-CRP
create_forest_plot("hs-CRP", "Standard", scale_factor = 10, x_min = -5, x_max = 5, 
                   Subgroup = c("Elapsed_time_days", "MonthsAfterTRT", "Assay_type", "Bias_judgment", "Study_design", "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Subgroups for the meta analysis
                   Label = c("Days after treatment", "MonthsAfterTRT", "Assay type", "Bias judgment", "Study design", "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Labels for the subgroups
                   Moderators = c("Elapsed_time_days","Mean_age", "Study_Yr", "PRE_mean_N_teeth", "With_Periodontitis", "FemaleRatio","With_Diabetes", "Smokers"))                             # Meta-regression moderators
create_forest_plot("hs-CRP", "Intensive", scale_factor = 10, x_min = -5, x_max = 5, 
                   Subgroup = c("Elapsed_time_days", "MonthsAfterTRT", "Assay_type", "Bias_judgment", "Study_design", "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Subgroups for the meta analysis
                   Label = c("Days after treatment", "MonthsAfterTRT", "Assay type", "Bias judgment", "Study design", "Additional_devices", "Female_binary","Smoke_binary", "Diabetes_binary"),   # Labels for the subgroups
                   Moderators = c("Elapsed_time_days","Mean_age", "Study_Yr", "PRE_mean_N_teeth", "With_Periodontitis", "With_Diabetes", "Smokers"))                             # Meta-regression moderators


