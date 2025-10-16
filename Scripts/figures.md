
### Figure 3
```
{r}
setwd("~/Library/CloudStorage/Box-Box/Widgingarri/2025/R_stuff")

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(tools)

metadata <- read.table("/Users/dawnlewis/Library/CloudStorage/Box-Box/Widgingarri/2025/R_stuff/Sample_metadata.txt",
  sep = "\t", 
  header = TRUE, 
  fill = TRUE, 
  stringsAsFactors = FALSE)

geochem <- read.table("~/Library/CloudStorage/Box-Box/Widgingarri/2025/R_stuff/Sediment_Report.txt",
  sep = "\t", 
  header = TRUE, 
  fill = TRUE, 
  stringsAsFactors = FALSE)

sediment_geol <- read.table("/Users/dawnlewis/Library/CloudStorage/Box-Box/Widgingarri/2025/R_stuff/sediment_geol.txt",   sep = "\t", 
  header = TRUE, 
  fill = TRUE, 
  stringsAsFactors = FALSE)

## Long format for line plots
data_variables <- geochem %>%
  pivot_longer(cols = c(pH, Water_pc), names_to = "Variable", values_to = "Value")


# Sort by depth
data_variables <- data_variables %>% arrange(Depth_m)

# --- Water plot ---
p_water <- ggplot(data_variables %>% filter(Variable == "Water_pc" & !is.na(Value)),
                  aes(x = Value, y = Depth_m)) +
  geom_path(color = "blue", linewidth = 0.8) +
  geom_point(color = "blue") +
  scale_y_reverse() +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Moisture", x = "Water %", y = "Depth (m)") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5))

# --- pH plot ---
p_ph <- ggplot(data_variables %>% filter(Variable == "pH" & !is.na(Value)),
               aes(x = Value, y = Depth_m)) +
  geom_path(color = "red", linewidth = 0.8) +
  geom_point(color = "red") +
  scale_y_reverse() +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Acidity", x = "pH", y = NULL) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5))

# --- Detailed sediment (geology) plot ---
detailed_cols <- c(
  "VERY.COARSE.SAND....", "COARSE.SAND....", "MEDIUM.SAND....", "FINE.SAND....", "VERY.FINE.SAND....",
  "VERY.COARSE.SILT....", "COARSE.SILT....", "MEDIUM.SILT....", "FINE.SILT....", "VERY.FINE.SILT....", "CLAY...."
)

detailed_data <- sediment_geol %>%
  select(`Depth..m.`, all_of(detailed_cols)) %>%
  pivot_longer(cols = all_of(detailed_cols), names_to = "GrainType", values_to = "Value") %>%
  group_by(`Depth..m.`) %>%
  mutate(Percent = Value / sum(Value) * 100) %>%
  ungroup()

# --- Custom colors ---
grain_colors <- c(
  "VERY.COARSE.SAND...." = "#8c6d47",
  "COARSE.SAND...." = "#a47f53",
  "MEDIUM.SAND...." = "#bb925e",
  "FINE.SAND...." = "#d3a46a",
  "VERY.FINE.SAND...." = "#eab676",
  "VERY.COARSE.SILT...." = "#d3d0cf",
  "COARSE.SILT...." = "#bcb8b6",
  "MEDIUM.SILT...." = "#908986",
  "FINE.SILT...." = "#7a716e",
  "VERY.FINE.SILT...." = "#645a56",
  "CLAY...." = "#873e23"
)

# --- Clean legend labels ---
detailed_data <- detailed_data %>%
  mutate(GrainType_clean = GrainType %>%
           gsub("\\.+", " ", .) %>%
           toTitleCase())

# --- Match factor levels for consistent fill colors ---
detailed_data$GrainType_clean <- factor(detailed_data$GrainType_clean,
                                        levels = toTitleCase(gsub("\\.+", " ", names(grain_colors))))
detailed_data$GrainType <- factor(detailed_data$GrainType, levels = names(grain_colors))

# --- Sediment plot ---
p_geol <- ggplot(detailed_data, aes(x = Percent, y = factor(`Depth..m.`), fill = GrainType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = grain_colors,
                    labels = levels(detailed_data$GrainType_clean)) +
  labs(title = "Sediment Composition", x = "Percent", y = NULL, fill = "Grain Type") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))

# --- Combine plots with relative widths (water & pH skinnier) ---
geochem <- (p_water | p_ph | p_geol) + plot_layout(widths = c(0.4, 0.4, 1))

# --- Save figure ---
ggsave("/Users/dawnlewis/Library/CloudStorage/Box-Box/Writing/4. Widgingarri/waterpH_sediment.png",
       plot = geochem, width = 14, height = 6, dpi = 300)

```

Figure X
