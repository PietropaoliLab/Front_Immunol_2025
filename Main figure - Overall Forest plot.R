library(readxl)
library(grid)
library(forestploter)
library(readr)
library(here)
MetaAnalysis_results <- data.frame(read_csv(here("AQ Projects/IL1 Meta-Analysis/Results/MetaAnalysis_results.csv")))
MetaAnalysis_results <- subset(MetaAnalysis_results, Subgroup == "Elapsed_time_days")

dim <- nrow(MetaAnalysis_results)
res <- data.frame()

# Iterate through rows
for (d in 1:(dim - 1)) {
  # Compare the first three columns of the current and next row
  if (!all(MetaAnalysis_results[d, 1:3] == MetaAnalysis_results[d + 1, 1:3])) {
    # Bind the row to the results dataframe
    res <- rbind(res, MetaAnalysis_results[d, ])
  }
}

# Add the last row if it hasn't been added already
res <- rbind(res, MetaAnalysis_results[dim, ])

# Remove subgroum column
res <- res[-3]



# Initialize an empty dataframe with the desired structure
mar <- data.frame()

# Iterate over unique Cytokine values
for (cytokine in unique(res$Cytokine)) {
  # Extract rows corresponding to the current Cytokine
  cytokine_data <- res[res$Cytokine == cytokine, ]
  
  # Add a row with the Cytokine name in the Subgroup column and NA in all other columns
  mar <- rbind(mar, data.frame(Subgroup = cytokine, Studies = NA, Subjects = NA, SMD_random = NA, seSMD_random = NA,
                               SMD_random_lower = NA, SMD_random_upper = NA, Pvalue_random = NA,
                               SMD_common = NA, seSMD_common = NA, SMD_common_lower = NA, SMD_common_upper = NA,
                               Pvalue_common = NA, stringsAsFactors = FALSE))
  
  # Iterate over the rows of cytokine_data and add them to mar
  for (i in 1:nrow(cytokine_data)) {
    mar <- rbind(mar, data.frame(
      Subgroup = cytokine_data$Treatment[i],
      Studies = cytokine_data$NumberOfStudis[i],
      Subjects = cytokine_data$NumberOfSubjects[i],
      SMD_random = cytokine_data$SMD_random[i],
      seSMD_random = cytokine_data$seSMD_random[i], # 
      SMD_random_lower = cytokine_data$SMD_random_lower[i],
      SMD_random_upper = cytokine_data$SMD_random_upper[i],
      Pvalue_random = cytokine_data$Pvalue_random[i],
      SMD_common = cytokine_data$SMD_common[i],
      seSMD_common = cytokine_data$seSMD_common[i], # 
      SMD_common_lower = cytokine_data$SMD_common_lower[i],
      SMD_common_upper = cytokine_data$SMD_common_upper[i],
      Pvalue_common = cytokine_data$Pvalue_common[i],
      stringsAsFactors = FALSE
    ))
  }
}




# indent the subgroup if there is a number in the placebo column
mar$Subgroup <- ifelse(is.na(mar$Studies), 
                         mar$Subgroup,
                         paste0("   ", mar$Subgroup))

# NA to blank or NA will be transformed to carachter.
mar$Studies <- ifelse(is.na(mar$Studies), "", mar$Studies)
mar$Subjects <- ifelse(is.na(mar$Subjects), "", mar$Subjects)


# Add a blank column for the forest plot to display CI.
# Adjust the column width with space, and increase the number of spaces below 
# to have a larger area to draw the CI. 
mar$`Common effect` <- paste(rep(" ", 20), collapse = " ")



# Common Effect - Create a confidence interval column to display
mar$`Common effect\nSMD (95% CI)` <- ifelse(is.na(mar$Studies), "",
                           sprintf("%.2f (%.2f to %.2f)",
                                   mar$SMD_common, mar$SMD_common_lower, mar$SMD_common_upper))

mar$`Common effect\nSMD (95% CI)` <- ifelse(mar$`Common effect\nSMD (95% CI)` == "NA (NA to NA)", "", mar$`Common effect\nSMD (95% CI)`)
mar$Pvalue_common <- ifelse(mar$Pvalue_common <0.001, "<0.001", round(mar$Pvalue_common,3))
mar$`Common effect\nSMD (95% CI)` <- paste(mar$`Common effect\nSMD (95% CI)`, "\nP=", mar$Pvalue_common, sep = "")
mar$`Common effect\nSMD (95% CI)` <- ifelse(mar$`Common effect\nSMD (95% CI)` == "\nP=NA", "", mar$`Common effect\nSMD (95% CI)`)



# Add a blank column for the forest plot to display CI.
# Adjust the column width with space, and increase the number of spaces below 
# to have a larger area to draw the CI. 
mar$`Random effect` <- paste(rep(" ", 20), collapse = " ")



# Random Effect - Create a confidence interval column to display
mar$`Random effect\nSMD (95% CI)` <- ifelse(is.na(mar$Studies), "",
                                    sprintf("%.2f (%.2f to %.2f)",
                                            mar$SMD_random, mar$SMD_random_lower, mar$SMD_random_upper))

mar$`Random effect\nSMD (95% CI)` <- ifelse(mar$`Random effect\nSMD (95% CI)` == "NA (NA to NA)", "", mar$`Random effect\nSMD (95% CI)`)
mar$Pvalue_random <- ifelse(mar$Pvalue_random <0.001, "<0.001", round(mar$Pvalue_random,3))
mar$`Random effect\nSMD (95% CI)` <- paste(mar$`Random effect\nSMD (95% CI)`, "\nP=", mar$Pvalue_random, sep = "")
mar$`Random effect\nSMD (95% CI)` <- ifelse(mar$`Random effect\nSMD (95% CI)` == "\nP=NA", "", mar$`Random effect\nSMD (95% CI)`)


rm(cytokine, cytokine_data, d, dim, i, res)

mar$Legend <- mar$Subgroup
MyVars <- c("TNF-a", "IL-1beta", "IL-6", "CRP", "hs-CRP")

mar$Legend <- ifelse(mar$Legend %in% MyVars, "", mar$Legend)
mar$Legend <- gsub("   ", "", mar$Legend)


# Set-up theme
library(extrafont)
loadfonts()
tm <- forest_theme(base_size = 10,
                   base_family = "Helvetica",
                   refline_col = "#aaaaaa",
                   title_gp = gpar(fontface = "bold", col = "black", fontfamily = "Helvetica")
)

p <- forest(mar[,c(1:3, 14:15, 16:17)],
       est = list(mar$SMD_common,
                  mar$SMD_random),
       lower = list(mar$SMD_common_lower,
                    mar$SMD_random_lower), 
       upper = list(mar$SMD_common_upper,
                    mar$SMD_random_upper),
       sizes = list(1/(mar$seSMD_common*15),
                    1/(mar$seSMD_random*13)),
       
       ci_column = c(4, 6),
       ref_line = 0,
       nudge_y = 0,
       xlim = list(c(-1, .5),
                   c(-2.5, .5)),
       xlab = c("SMD", "SMD"),
       
       arrow_lab = c("Decrease after TRT", "Increase after TRT"),
       theme = tm)

g <- edit_plot(p, row = c(3, 6, 9, 12, 15), gp = gpar(col = "#e41a1c")) # Intensive
g <- edit_plot(g, row = c(2, 5, 8, 11, 14), gp = gpar(col = "#377eb8"))     # Standard
g <- edit_plot(g, row = c(1, 4, 7, 10, 13), gp = gpar(fontfamily = "Helvetica"))
g <- edit_plot(g, part = "header", row = 1, gp = gpar(fontfamily = "Helvetica"))

print(g)
here()
pdf(file = here("AQ Projects/IL1 Meta-Analysis/MainFigure_Forest_plot.pdf"), width = 12, height = 8, family = "sans")
print(g)
dev.off()
