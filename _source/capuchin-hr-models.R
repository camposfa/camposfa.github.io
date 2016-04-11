# ---- prepare_model_workspace --------------------------------------------

Sys.setenv(TZ = 'UTC')
list.of.packages <- list("lme4", "plyr", "lubridate", "scales", "reshape2",  
                         "ggplot2", "RColorBrewer", "MuMIn", "multcomp",
                         "colorspace")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(unlist(new.packages))
lapply(list.of.packages, require, character.only = T)

options(na.action = "na.fail")

# ---- load_hr_data -------------------------------------------------------

# Read hr data from csv file created 2013-10-09
# Skips the entire home range calculation step
# Must redo if new data are added
ob <- read.csv("_data/capuchin-home-ranges/ob.csv")
ob$ints <- interval(start = ob$start, end = ob$end)

### Subset
mon <- subset(ob, block_type == "month")
qua <- subset(ob, block_type == "quarter")
hal <- subset(ob, block_type == "half")
yea <- subset(ob, block_type == "year")


# ---- area_models --------------------------------------------------------

# Loop over scale/zone combinations and apply model selection
# Not a very elegant solution, but whatevs..
df_set <- list(mon, qua, hal)
df_set_names <- c("mon", "qua", "hal")
y_set <- c("area_core", "area_primary", "area_total")

# Blank dfs for storing results
area_model_results <- NULL
area_selection_table <- NULL

for (i in 1:length(df_set)) {
  for (j in 1:length(y_set)) {
    
    # Run models
    m1 <- lmer(log(df_set[[i]][, y_set[j]]) ~ sqrt(nb_reloc) + 
                 (1|id), 
               REML = FALSE, data = df_set[[i]])
    
    m2 <- lmer(log(df_set[[i]][, y_set[j]]) ~ sqrt(nb_reloc) + 
                 mean_tmax_s + (1|id), 
               REML = FALSE, data = df_set[[i]])
    
    m3 <- lmer(log(df_set[[i]][, y_set[j]]) ~ sqrt(nb_reloc) + 
                 mean_fruit_s + (1|id), 
               REML = FALSE, data = df_set[[i]])
    
    m4 <- lmer(log(df_set[[i]][, y_set[j]]) ~ sqrt(nb_reloc) + 
                 adult_mass_s + (1|id), 
               REML = FALSE, data = df_set[[i]])
    
    m5 <- lmer(log(df_set[[i]][, y_set[j]]) ~ sqrt(nb_reloc) + 
                 mean_tmax_s + mean_fruit_s + (1|id), 
               REML = FALSE, data = df_set[[i]])
    
    m6 <- lmer(log(df_set[[i]][, y_set[j]]) ~ sqrt(nb_reloc) + 
                 mean_tmax_s + adult_mass_s + (1|id), 
               REML = FALSE, data = df_set[[i]])
    
    m7 <- lmer(log(df_set[[i]][, y_set[j]]) ~ sqrt(nb_reloc) + 
                 mean_fruit_s + adult_mass_s + (1|id), 
               REML = FALSE, data = df_set[[i]])
    
    m8 <- lmer(log(df_set[[i]][, y_set[j]]) ~ sqrt(nb_reloc) + 
                 mean_tmax_s + mean_fruit_s + adult_mass_s + (1|id), 
               REML = FALSE, data = df_set[[i]])
    
    # Create model selection table dataframe
    temp <- suppressWarnings(model.sel(list(m1, m2, m3, m4, m5, m6, m7, m8)))
    temp$m_name <- row.names(temp)
    temp$data_set <- rep(df_set_names[i], times = nrow(temp))
    temp$zone <- rep(y_set[j], times = nrow(temp))
    
    # Add to area_selection_table
    area_selection_table <- rbind(area_selection_table, data.frame(temp))
    
    # Obtain average model
    fma <- suppressWarnings(model.avg(temp, fit = TRUE))
    avg_model <- data.frame(summary(fma)[["coefmat.subset"]])
    avg_model <- cbind(term = rownames(avg_model), avg_model)
    avg_model$term <- rownames(avg_model)
    rownames(avg_model) <- NULL
    ci <- data.frame(confint(fma))
    ci$term <- rownames(ci)
    rownames(ci) <- NULL
    avg_model <- inner_join(avg_model, ci, by = "term")
    imp <- melt(fma$importance)
    imp$term <- rownames(imp)
    rownames(imp) <- NULL
    avg_model <- inner_join(avg_model, imp, by = "term")
    names(avg_model) <- c("term", "estimate", "std_error", "adjusted_se", "Z Value",
                          "P Value", "lower_ci", "upper_ci", "importance")
    
    # Add to area_model_results
    avg_model$response <- y_set[j]
    avg_model$scale <- df_set_names[i]
    area_model_results <- rbind(area_model_results, avg_model)
  }
}


# ---- rename_area_vars ---------------------------------------------------

area_model_results$response <- 
  factor(revalue(area_model_results$response, 
                 c("area_core" = "Core",
                   "area_primary" = "Primary",
                   "area_total" = "Total")),
         levels = c("Core",
                    "Primary",
                    "Total"))

area_model_results$scale <- 
  factor(revalue(area_model_results$scale, 
                 c("mon" = "Monthly",
                   "qua" = "Quarterly",
                   "hal" = "Half-yearly")),
         levels = c("Monthly",
                    "Quarterly",
                    "Half-yearly"))

area_model_results$term <- 
  factor(revalue(area_model_results$term, 
                 c("adult_mass_s" = "Group mass",
                   "mean_fruit_s" = "Fruit biomass",
                   "mean_tmax_s" = "Mean max temperature",
                   "sqrt(nb_reloc)" = "Sqrt num locations")),
         levels = c("Sqrt num locations",
                    "Mean max temperature", 
                    "Fruit biomass", 
                    "Group mass"))


# ---- plot_area_results --------------------------------------------------

# Exclude half-year
ggplot(subset(area_model_results, scale != "Half-yearly" & 
                term != "Sqrt num locations"), 
       aes(x = term, color = importance)) + 
  geom_hline(aes(yintercept = 0), lty = 2, color = "gray50") + 
  geom_errorbar(aes(ymax = lower_ci, 
                    ymin = upper_ci),
                width = 0.2, size = 0.75) + 
  geom_segment(aes(y = estimate + std_error, 
                   yend = estimate - std_error, 
                   xend = term), 
               size = 1.5) +
  geom_point(aes(y = estimate), 
             size = 1.5, 
             color = "black") +
  labs(y = "Regression estimates", x = "") + 
  scale_color_gradientn(name = "Importance", 
                        colours = rev(heat_hcl(12, 
                                             c = c(80, 30), 
                                             l = c(30, 90), 
                                             power = c(1/5, 2))), 
                        values = seq(0, 1, length.out = 12), 
                        rescaler = function(x, ...) x, 
                        oob = identity, limits = c(0,1)) +
  theme_bw() +
  theme(legend.position = "bottom", 
        panel.border = element_rect(fill = NA, colour = "grey50"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) +
  facet_grid(scale ~ response) +
  coord_flip()



# ---- ndvi_models --------------------------------------------------------

df_set <- list(mon, qua, hal)
df_set_names <- c("mon", "qua", "hal")
y_set <- c("ndvi_high", "ndvi_medium", "ndvi_low")

# Blank dfs for storing results
ndvi_model_results <- NULL
ndvi_selection_table <- NULL

for (i in 1:length(df_set)) {
  for (j in 1:length(y_set)) {
    
    # Run models
    m1 <- lmer(log(df_set[[i]][, y_set[j]]) ~ 1 + (1|id), 
               REML = FALSE, data = df_set[[i]])
    
    m2 <- lmer(log(df_set[[i]][, y_set[j]]) ~ mean_tmax_s + (1|id), 
               REML = FALSE, data = df_set[[i]])
    
    m3 <- lmer(log(df_set[[i]][, y_set[j]]) ~ mean_fruit_s + (1|id), 
               REML = FALSE, data = df_set[[i]])
    
    m4 <- lmer(log(df_set[[i]][, y_set[j]]) ~ adult_mass_s + (1|id), 
               REML = FALSE, data = df_set[[i]])
    
    m5 <- lmer(log(df_set[[i]][, y_set[j]]) ~ mean_tmax_s + mean_fruit_s + 
                 (1|id), REML = FALSE, data = df_set[[i]])
    
    m6 <- lmer(log(df_set[[i]][, y_set[j]]) ~ mean_tmax_s + adult_mass_s + 
                 (1|id), REML = FALSE, data = df_set[[i]])
    
    m7 <- lmer(log(df_set[[i]][, y_set[j]]) ~ mean_fruit_s + adult_mass_s + 
                 (1|id), REML = FALSE, data = df_set[[i]])
    
    m8 <- lmer(log(df_set[[i]][, y_set[j]]) ~ mean_tmax_s + mean_fruit_s + 
                 adult_mass_s + (1|id), REML = FALSE, data = df_set[[i]])
    
    # Create model selection table dataframe
    temp <- suppressWarnings(model.sel(list(m1, m2, m3, m4, m5, m6, m7, m8)))
    temp$m_name <- row.names(temp)
    temp$data_set <- rep(df_set_names[i], times = nrow(temp))
    temp$zone <- rep(y_set[j], times = nrow(temp))
    
    # Add to ndvi_selection_table
    ndvi_selection_table <- rbind(ndvi_selection_table, data.frame(temp))
    
    # Obtain average model
    fma <- suppressWarnings(model.avg(temp, fit = TRUE))
    avg_model <- data.frame(summary(fma)[["coefmat.subset"]])
    avg_model <- cbind(term = rownames(avg_model), avg_model)
    avg_model$term <- rownames(avg_model)
    rownames(avg_model) <- NULL
    ci <- data.frame(confint(fma))
    ci$term <- rownames(ci)
    rownames(ci) <- NULL
    avg_model <- inner_join(avg_model, ci, by = "term")
    imp <- melt(fma$importance)
    imp$term <- rownames(imp)
    rownames(imp) <- NULL
    avg_model <- inner_join(avg_model, imp, by = "term")
    names(avg_model) <- c("term", "estimate", "std_error", "adjusted_se", "Z Value",
                          "P Value", "lower_ci", "upper_ci", "importance")
    
    # Add to area_model_results
    avg_model$response <- y_set[j]
    avg_model$scale <- df_set_names[i]
    ndvi_model_results <- rbind(ndvi_model_results, avg_model)
  }
}


# ---- rename_ndvi_vars ---------------------------------------------------

ndvi_model_results$response <- 
  factor(revalue(ndvi_model_results$response, 
                 c("ndvi_high" = "High Use",
                   "ndvi_medium" = "Medium Use",
                   "ndvi_low" = "Low Use")),
         levels = c("High Use",
                    "Medium Use",
                    "Low Use"))

ndvi_model_results$scale <- 
  factor(revalue(ndvi_model_results$scale, 
                 c("mon" = "Monthly",
                   "qua" = "Quarterly",
                   "hal" = "Half-yearly")),
         levels = c("Monthly",
                    "Quarterly",
                    "Half-yearly"))

ndvi_model_results$term <- 
  factor(revalue(ndvi_model_results$term, 
                 c("adult_mass_s" = "Group mass",
                   "mean_fruit_s" = "Fruit biomass",
                   "mean_tmax_s" = "Mean max temperature")), 
         levels = c("Mean max temperature", 
                    "Fruit biomass", 
                    "Group mass"))


# ---- plot_ndvi_results --------------------------------------------------

# Exclude half-year
ggplot(subset(ndvi_model_results, scale != "Half-yearly"), 
       aes(x = term, color = importance)) + 
  geom_hline(aes(yintercept = 0), lty = 2, color = "gray50") + 
  geom_errorbar(aes(ymax = lower_ci, 
                    ymin = upper_ci),
                width = 0.2, size = 0.75) + 
  geom_segment(aes(y = estimate + std_error, 
                   yend = estimate - std_error, 
                   xend = term), 
               size = 1.5) +
  geom_point(aes(y = estimate), 
             size = 1.5, 
             color = "black") +
  labs(y = "Regression estimates", x = "") + 
  scale_color_gradientn(name = "Importance", 
                        colours = rev(heat_hcl(12, 
                                             h = c(130, 70),
                                             c = c(80, 30), 
                                             l = c(30, 90), 
                                             power = c(1/5, 2))), 
                        values = seq(0, 1, length.out = 12), 
                        rescaler = function(x, ...) x, 
                        oob = identity, limits = c(0,1)) +
  theme_bw() +
  theme(legend.position = "bottom", 
        panel.border = element_rect(fill = NA, colour = "grey50"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) +
  facet_grid(scale ~ response) +
  coord_flip()


# ---- zones_model ---------------------------------------------------------

df_set <- list(mon, qua, hal, yea)
df_set_names <- c("mon", "qua", "hal", "yea")

# Blank dfs for storing results
zones_model_results <- NULL
zones_selection_table <- NULL

for (i in 1:length(df_set)) {
  zones <- melt(df_set[[i]], 
                id.vars = c("id", "block_type", "rowid"), 
                measure.vars = c("ndvi_high", "ndvi_medium", "ndvi_low"))
  
  zones$variable <- factor(zones$variable, levels = c("ndvi_medium",
                                                      "ndvi_low",
                                                      "ndvi_high"))
  
  m1 <- lmer(value ~ variable + (1|rowid) + (1|id), data = zones)
  temp <- dredge(m1, REML = FALSE)
  
  # Multiple Comparisons of Means: Tukey Contrasts
  summary(glht(m1, linfct = mcp(variable = "Tukey")))
  
  # Create model selection table dataframe
  temp$m_name <- row.names(temp)
  temp$data_set <- rep(df_set_names[i], times = nrow(temp))
  
  # Add to zones_selection_table
  zones_selection_table <- rbind(zones_selection_table, data.frame(temp))
  
  # Obtain average model
  fma <- suppressWarnings(model.avg(temp, fit = TRUE))
  avg_model <- data.frame(summary(fma)[["coefmat.subset"]])
  avg_model <- cbind(term = rownames(avg_model), avg_model)
  avg_model$term <- rownames(avg_model)
  rownames(avg_model) <- NULL
  ci <- data.frame(confint(fma))
  ci$term <- rownames(ci)
  rownames(ci) <- NULL
  avg_model <- inner_join(avg_model, ci, by = "term")
  names(avg_model) <- c("term", "estimate", "std_error", "adjusted_se", "Z Value",
                        "P Value", "lower_ci", "upper_ci")
  
  # Add to zones_model_results
  avg_model$scale <- df_set_names[i]
  zones_model_results <- rbind(zones_model_results, avg_model)
}

# Remove intercept
zones_model_results <- subset(zones_model_results, term != "(Intercept)")


# ---- rename_zones_vars ---------------------------------------------------

zones_model_results$scale <- 
  factor(revalue(zones_model_results$scale, 
                 c("mon" = "Monthly",
                   "qua" = "Quarterly",
                   "hal" = "Half-yearly",
                   "yea" = "Yearly")),
         levels = c("Monthly",
                    "Quarterly",
                    "Half-yearly",
                    "Yearly"))

zones_model_results$term <- 
  factor(revalue(zones_model_results$term, 
                 c("variablendvi_low" = "Low Use", 
                   "variablendvi_high" = "High Use")), 
         levels = c("Low Use", 
                    "High Use"))


# ---- plot_zones_results --------------------------------------------------

ggplot(zones_model_results, aes(x = term)) + 
  geom_hline(aes(yintercept = 0), lty = 2, color = "gray50") + 
  geom_errorbar(aes(ymax = lower_ci, ymin = upper_ci),
                width = 0.2, size = 0.75, color = "gray10") + 
  geom_segment(aes(y = estimate + std_error, 
                   yend = estimate - std_error, 
                   xend = term), 
               size = 1.5, color = "gray10") +
  geom_point(aes(y = estimate), 
             size = 1.5, 
             color = "white") +
  scale_y_continuous(limits = c(-0.027, 0.027)) +
  labs(y = "Regression estimates", x = "") + 
  theme_bw() +
  theme(legend.position = "bottom", 
        panel.border = element_rect(fill = NA, colour = "grey50")) +
  facet_grid(scale ~ .) +
  coord_flip()


# ---- model_cleanup ------------------------------------------------------

rm(m1)
rm(m2)
rm(m3)
rm(m4)
rm(m5)
rm(m6)
rm(m7)
rm(m8)
rm(temp)
rm(avg_model)
rm(imp)
rm(fma)
gc()