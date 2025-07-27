# ================================
# Meta-Analysis and Meta-Regression Workflow
# ================================

# Load libraries
library(here)
library(readxl)
library(meta)
library(metafor)
library(tidyverse)
library(splines)
library(ggplot2)

# -------------------------------
# 1. Load and prepare data
# -------------------------------
data_path <- here("AQ Projects/IL1 Meta-Analysis/Dataset.xlsx")
db <- read_excel(data_path)

# Split by treatment type
standard_df  <- filter(db, Type_of_Treatment == "Standard")
intensive_df <- filter(db, Type_of_Treatment == "Intensive")

# -------------------------------
# 2. Meta-analysis per group
# -------------------------------
run_metacont <- function(df) {
  metacont(
    n.e = Subjects, mean.e = POST_mean, sd.e = POST_SD,
    n.c = Subjects, mean.c = PRE_mean, sd.c = PRE_SD,
    data = df,
    studlab = paste(Study_ID, Cytokine),
    sm = "SMD", method.tau = "DL",
    comb.fixed = TRUE, comb.random = TRUE
  )
}

meta_standard  <- run_metacont(standard_df)
meta_intensive <- run_metacont(intensive_df)

# Combined analysis with subgroup
meta_all <- metacont(
  n.e = Subjects, mean.e = POST_mean, sd.e = POST_SD,
  n.c = Subjects, mean.c = PRE_mean, sd.c = PRE_SD,
  data = db,
  studlab = paste(Study_ID, Cytokine),
  sm = "SMD", byvar = Type_of_Treatment,
  method.tau = "DL", comb.fixed = TRUE, comb.random = TRUE
)

meta_standard
meta_intensive
meta_all

# -------------------------------
# 3. Forest plots
# -------------------------------
forest(meta_standard, sortvar = TE, prediction = TRUE, print.tau2 = TRUE,
       col.random = "blue", col.fixed = "darkgreen", main = "Standard TRT")

forest(meta_intensive, sortvar = TE, prediction = TRUE, print.tau2 = TRUE,
       col.random = "red", col.fixed = "darkgreen", main = "Intensive TRT")

forest(meta_all, bylab = "Treatment Type", prediction = TRUE, print.tau2 = TRUE,
       col.random = "purple", col.fixed = "green4", main = "Standard vs Intensive")

# -------------------------------
# 4. Meta-analysis with metafor
# -------------------------------
meta_data <- data.frame(
  TE = meta_all$TE,
  seTE = meta_all$seTE,
  studlab = meta_all$studlab
)

res_metafor <- rma(yi = TE, sei = seTE, data = meta_data, method = "DL")

forest(res_metafor,
       order = order(meta_data$TE),
       cex = 0.2,
       slab = meta_data$studlab,
       xlab = "Standardized Mean Difference (SMD)")

# -------------------------------
# 5. Meta-regression with splines
# -------------------------------
foo <- data.frame(
  TE = meta_intensive$TE,
  seTE = meta_intensive$seTE,
  Elapsed_time_days = meta_intensive$data$Elapsed_time_days
)

foo$TimeTRT <- as.numeric(cut(foo$Elapsed_time_days,
                              breaks = c(0, 42, 90, 180, Inf),
                              labels = 1:4))

meta_regression_spline <- rma(yi = TE, sei = seTE,
                              mods = ~ ns(Elapsed_time_days, df = 3),
                              data = foo, method = "DL")

summary(meta_regression_spline)

# -------------------------------
# 6. Prediction and plotting
# -------------------------------
new_days <- data.frame(Elapsed_time_days = seq(min(foo$Elapsed_time_days, na.rm = TRUE),
                                               max(foo$Elapsed_time_days, na.rm = TRUE),
                                               length.out = 100))

spline_terms <- model.matrix(~ ns(Elapsed_time_days, df = 3), data = new_days)[, -1]
preds <- predict(meta_regression_spline, newmods = spline_terms)

plot_data <- data.frame(
  Elapsed_time_days = new_days$Elapsed_time_days,
  SMD = preds$pred,
  CI_lb = preds$ci.lb,
  CI_ub = preds$ci.ub
)

# -------------------------------
# 7. Plot with ggplot2
# -------------------------------
MyTheme <- theme_bw(base_size = 11) +
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "none",
        aspect.ratio = 1,
        legend.key = element_blank(),
        legend.background = element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Helvetica"))

ggplot(data = plot_data, aes(x = Elapsed_time_days, y = SMD)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.2) +
  geom_line(color = "#e41a1c", linewidth = 1.2) +
  geom_ribbon(aes(ymin = CI_lb, ymax = CI_ub), fill = "#e41a1c", alpha = 0.05) +
  coord_cartesian(ylim = c(-1.5, 0.7), expand = FALSE) +
  labs(x = "Days after treatment", 
       y = "Effect on TNF-a, IL1b, IL-6, (hs)CRP\nafter intensive TRT (SMD)") +
  MyTheme


foo <- rbind(
data.frame(
"TE" = meta_intensive$TE.random,
"Lower" = meta_intensive$lower.random,
"Upper" = meta_intensive$upper.random,
"seTE" = meta_intensive$seTE.random,
"Pval" = meta_intensive$pval.random,
"Label" = "Intensive"),

data.frame(
  "TE" =  meta_standard$TE.random,
  "Lower" = meta_standard$lower.random,
  "Upper" =meta_standard$upper.random,
  "seTE" = meta_standard$seTE.random,
  "Pval" = meta_standard$pval.random,
  "Label" = "Standard"))

foo <- rbind(foo,
data.frame(
  "TE" =  meta_all$TE.random,
  "Lower" = meta_all$lower.random,
  "Upper" = meta_all$upper.random,
  "seTE" = meta_all$seTE.random,
  "Pval" = meta_all$pval.random,
  "Label" = "Overall"))


MyTheme <- theme_bw(base_size = 11) +
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "none",
        # aspect.ratio = 1,
        legend.key = element_blank(),
        legend.background = element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Helvetica"))


ggplot(foo, aes(x=TE, y=Label, color = Label))+
  geom_vline(xintercept = 0, linetype = 2, linewidth = .2)+
  
   # geom_vline(xintercept = -0.3192194, linetype = 1, linewidth = .1, color = "#377eb8")+
   # geom_vline(xintercept = -0.4090076, linetype = 1, linewidth = .1, color = "gray70")+
   # geom_vline(xintercept = -0.4622603, linetype = 1, linewidth = .1, color = "#e41a1c")+
  
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0, linewidth = .2)+
  geom_point(aes(size = 1/seTE), shape = 15)+
  scale_size_continuous(range = c(5,13))+
  scale_color_manual(values = c("#e41a1c", "gray70","#377eb8"))+
  scale_x_continuous(breaks = c(-0.6, -0.46, -0.41, -0.32, 0))+
  coord_cartesian(xlim = c(-0.6, 0))+
  labs(x="Effect on TNF-a, IL1b, IL-6, (hs)CRP\nafter TRTs (SMD)", y="")+
  MyTheme




library(cowplot)
MyTheme <- theme_bw(base_size = 11) +
  theme(axis.text = element_text(colour = "black", size = 10),
        panel.grid = element_blank(),
        legend.position = "none",
        aspect.ratio = 1,
        legend.key = element_blank(),
        legend.background = element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(family = "Helvetica"))




plot_grid(
  ggplot(foo, aes(x=TE, y=Label, color = Label))+
    geom_vline(xintercept = 0, linetype = 2, linewidth = .2)+
    
    # geom_vline(xintercept = -0.3192194, linetype = 1, linewidth = .1, color = "#377eb8")+
    # geom_vline(xintercept = -0.4090076, linetype = 1, linewidth = .1, color = "gray70")+
    # geom_vline(xintercept = -0.4622603, linetype = 1, linewidth = .1, color = "#e41a1c")+
    
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0, linewidth = .2)+
    geom_point(aes(size = 1/seTE), shape = 15)+
    scale_size_continuous(range = c(5,13))+
    scale_color_manual(values = c("#e41a1c", "gray70","#377eb8"))+
    # scale_x_continuous(breaks = c(-0.6, -0.46, -0.41, -0.32, 0))+
    coord_cartesian(xlim = c(-0.6, 0))+
    labs(x="Effect on TNF-a, IL1b, IL-6, (hs)CRP\nafter TRTs (SMD)", y="")+
    MyTheme,
  ggplot(data = plot_data, aes(x = Elapsed_time_days, y = SMD)) +
            geom_hline(yintercept = 0, linetype = 2, linewidth = 0.2) +
            geom_line(color = "#e41a1c", linewidth = 1.2) +
            geom_ribbon(aes(ymin = CI_lb, ymax = CI_ub), fill = "#e41a1c", alpha = 0.05) +
            coord_cartesian(ylim = c(-1.5, 0.7), expand = FALSE) +
            labs(x = "Days after treatment", 
                 y = "Effect on TNF-a, IL1b, IL-6, (hs)CRP\nafter intensive TRT (SMD)") +
            MyTheme, labels = "AUTO"
  
  )
  