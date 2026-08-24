
library(rgplates)
library(dplyr)                  
library(ggplot2)               
library(astrochron)
library(RColorBrewer)

## Location of the sites ##
Sites2 <- read.csv("RScripts&Data/Sites.csv", header=TRUE, stringsAsFactors=FALSE)
Sites2 <- Sites2[c("Site", "Longitude", "Latitude")]

target_sites2 <- c("Mandorah 1","Hadrian 1", "Crown 1", "Calliance 2", "Omar 1", "Ermine 1", "Brigadier 1", "Eastbrook 1", "Central Gorgon 1", "Bluebell 1", "York 1")

# Extract paleolatitudes for each site
site_coords2 <- Sites2[Sites2$Site %in% target_sites2, ]
paleolatitudes_list2 <- list()

for (i in 1:nrow(site_coords2)) {
  site_name2 <- site_coords2$Site[i]
  coord <- site_coords2[i, 2:3]
  ages <- seq(2.5, 21.6, 1)
  paleo <- reconstruct(coord, age = ages)
  paleolatitudes <- sapply(paleo, function(x) x[2])
  paleolatitudes_list2[[site_name2]] <- data.frame(Age = ages, Paleolatitude = paleolatitudes, Site = site_name2)
}

final_df2 <- data.frame()


# Get a list of NGR CSV files in the directory
csv_files <- list.files("RScripts&Data/Sites Data_Age-NGR/", pattern = "\\.csv$", full.names = TRUE)

for (file in csv_files) {
  site_name2 <- tools::file_path_sans_ext(basename(file))
  if (site_name2 %in% target_sites2) {
    data <- read.csv(file, stringsAsFactors = FALSE)
    data$Site <- site_name2
    
    # Match paleolatitude for same ages
    colnames(data)[tolower(colnames(data)) == "age"] <- "Age"
    colnames(data)[tolower(colnames(data)) == "gr"] <- "GR"
    
    data$Age <- data$Age / 1000
    
    # Interpolate paleolatitude to match Age in data
    paleo_df2 <- paleolatitudes_list2[[site_name2]]
    
    # Add interpolated Paleolatitude
    data$Paleolatitude <- approx(x = paleo_df2$Age, y = paleo_df2$Paleolatitude, xout = data$Age, rule = 2)$y
    
    # Keep relevant columns
    merged_df2 <- data[, c("Site", "Age", "GR", "Paleolatitude")]
    final_df2 <- rbind(final_df2, merged_df2)
    final_df2 <- final_df2 %>%
      filter(Age >= 2.5, Age <= 21.6)
  }
}
custom_colors <- colorRampPalette(brewer.pal(10,"RdBu"))(21)

ipcc_precip <- c("#543005", "#8C510A", "#BF812D", "#F5F5F5", 
                 "#A6DBD8", "#80CDC1", "#35978F")

label_df2 <- final_df2 %>%
  group_by(Site) %>%
  slice_min(Age, n = 1, with_ties = FALSE)

setwd("RScripts&Data/Time_Slice_Map")

pdf(file = "Fig 2.pdf", width = 11, height = 11, paper = "a4")
png(filename = "Fig 2.png", width = 9000, height = 9000, res = 600)

ggplot(final_df2, aes(x = Paleolatitude, y = Age, color = GR)) +
  scale_color_gradientn(colors = custom_colors, limits = c(0, 50), oob = scales::squish) +
  geom_point(size = 4) +
  scale_y_reverse(breaks = seq(0, 22, 2), expand = c(0, 0)) +  
  scale_x_continuous(breaks = seq(-35, -11, 2), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-31, -11), ylim = c(22, 0), clip = "on") + 
  labs(x = "Paleolatitude (°)", y = "Age (Ma)") +
  theme_minimal(base_size = 18) +
  theme(
    plot.background  = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.6),  
    axis.line.y = element_blank(),  
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length.x = unit(8, "pt"),
    axis.text = element_text(size = 26),
    axis.title = element_text(size = 26),
    plot.margin = margin(7, 5, 17, 20) 
  )
dev.off()



final_df2$Age_rounded <- round(final_df2$Age, 3)

mean_gr <- final_df2 %>%
  group_by(Age_rounded) %>%
  summarise(mean_GR = mean(GR, na.rm = TRUE))

setwd("RScripts&Data/Time_Slice_Map")

pdf(file = "Fig 3.pdf", width = 11, height = 5.5)
png(filename = "Fig 3.png", width = 9000, height = 5000, res = 600)

ggplot(mean_gr, aes(x = Age_rounded, y = 0, color = mean_GR)) +
  scale_color_gradientn(colors = ipcc_precip, limits = c(0, 50), oob = scales::squish) +
  geom_point(size = 8) +
  scale_x_reverse(limits = c(22, 2), breaks = seq(2, 22, 2),
                  expand = c(0, 0), position = "top") + 
  labs(x = "Age (Ma)", y = NULL) +
  guides(color = "none") +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.6),  
    axis.line.y = element_blank(),  
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 0.6),
    axis.title.x.top = element_text(margin = margin(b = 20)),
    axis.ticks.length.x = unit(8, "pt"),
    axis.text = element_text(color = "black", size = 26, ),
    axis.title = element_text(size = 26),
    plot.margin = margin(7, 15, 17, 25) 
  )
dev.off()
