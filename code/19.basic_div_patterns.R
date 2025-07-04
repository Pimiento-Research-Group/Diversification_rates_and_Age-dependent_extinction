library(divDyn)
library(dplyr)
library(ggplot2)
library(patchwork)

data("stages")

df <- read.table(file = "/Volumes/External_memory/Dropbox/Kristina_PhD_K's_version/Kristina's files/Analyses/PyRate/PyRate_Analysis/inputs/2025/June/species.txt", sep = "\t")

n_reps <- 10

results_list <- vector("list", n_reps)

set.seed(42)

# Function: sample a random age between max_ma and min_ma
sample_random_age <- function(min_ma, max_ma) {
  runif(1, min_ma, max_ma)
}

# Function: assign a stage number based on a single age
assign_stage <- function(age) {
  stage_row <- stages[stages$bottom >= age & stages$top < age, ]
  if (nrow(stage_row) > 0) {
    return(stage_row$stg[1])
  } else {
    return(NA)  # outside known bins
  }
}

results_list <- vector("list", n_reps)

for (i in 1:n_reps) {
  cat("Replicate", i, "\n")
  
  # Step 1: sample random ages
  df$rand_age <- mapply(sample_random_age, df$min_ma, df$max_ma)
  
  # Step 2: assign each random age to a stage
  df$stg <- sapply(df$rand_age, assign_stage)
  
  # Remove occurrences outside known stages (optional)
  df_i <- df[!is.na(df$stg), ]
  
  # Compute diversity dynamics
  out <- divDyn(df_i, tax="taxon_name", bin="stg")
  
  results_list[[i]] <- out
}


combined_results <- bind_rows(
  lapply(seq_along(results_list), function(i) {
    if (!is.null(results_list[[i]]) && is.data.frame(results_list[[i]]) && nrow(results_list[[i]]) > 0) {
      cbind(rep = i, results_list[[i]])
    }
  })
)

summary_rates <- combined_results %>%
  group_by(stg) %>%
  summarise(
    oriPC_mean = mean(oriPC, na.rm=TRUE),
    oriPC_lower = quantile(oriPC, probs=0.025, na.rm=TRUE),
    oriPC_upper = quantile(oriPC, probs=0.975, na.rm=TRUE),
    extPC_mean = mean(extPC, na.rm=TRUE),
    extPC_lower = quantile(extPC, probs=0.025, na.rm=TRUE),
    extPC_upper = quantile(extPC, probs=0.975, na.rm=TRUE),
    .groups = "drop"
  )

summary_rates <- summary_rates %>%
  mutate(
    net_div_mean = oriPC_mean - extPC_mean,
    net_div_lower = oriPC_lower - extPC_upper,
    net_div_upper = oriPC_upper - extPC_lower
  )

summary_rates$mid_ma <- stages$mid[match(summary_rates$stg, stages$stg)]

  # Extinction plot
p_ext <- ggplot(summary_rates, aes(x=mid_ma, y=extPC_mean)) +
  geom_line(color="red") +
  geom_ribbon(aes(ymin=extPC_lower, ymax=extPC_upper), fill="red", alpha=0.3) +
  scale_x_reverse() +
  xlim(145, 0) +
  labs(x = "", y="Per capita extinction rate") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Origination plot
p_ori <- ggplot(summary_rates, aes(x=mid_ma, y=oriPC_mean)) +
  geom_line(color="blue") +
  geom_ribbon(aes(ymin=oriPC_lower, ymax=oriPC_upper), fill="blue", alpha=0.3) +
  scale_x_reverse() +  
  xlim(145, 0) +
  coord_cartesian(xlim = c(145, 0)) +
  labs(x = "", y="Per capita origination rate") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_net <- ggplot(summary_rates, aes(x=mid_ma, y=net_div_mean)) +
  geom_line(color="grey40") +
  geom_ribbon(aes(ymin=net_div_lower, ymax=net_div_upper), fill="grey80", alpha=0.4) +
  xlim(145, 0) +
  labs(
    x = "Time (Ma)",
    y = "Net diversification rate",
  ) +
  theme_minimal()

# Combine into two-panel plot
p_combined <- p_ext / p_ori / p_net

ggsave("/Volumes/External_memory/Dropbox/Kristina_PhD_K's_version/Kristina's files/Analyses/PyRate/PyRate_Analysis/outputs/2025/June/species/extinction_origination_rates_raw.pdf", p_combined, width=6.69291, height=4, units="in", dpi=300)
