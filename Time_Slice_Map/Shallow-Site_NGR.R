
library(rgplates)
library(dplyr)                  
library(ggplot2)               
library(astrochron)
library(RColorBrewer)

## Location of the sites ##
Sites1 <- read.csv("RScripts&Data/Sites.csv", header=TRUE, stringsAsFactors=FALSE)
Sites1 <- Sites1[c("Site", "Longitude", "Latitude")]

target_sites1 <- c("Great Eastern 1", "Woodbine 1", "Kalyptea 1", "Dinichthys 1", "Walkley 1", "U1464", "Picard 1", "Finucane 1", "Fisher 1", "West Barrow 2")

# Extract paleolatitudes for each site
site_coords1 <- Sites1[Sites1$Site %in% target_sites1, ]
paleolatitudes_list1 <- list()

for (i in 1:nrow(site_coords1)) {
  site_name1 <- site_coords1$Site[i]
  coord <- site_coords1[i, 2:3]
  ages <- seq(2.5, 21.6, 1)
  paleo <- reconstruct(coord, age = ages)
  paleolatitudes <- sapply(paleo, function(x) x[2])
  paleolatitudes_list1[[site_name1]] <- data.frame(Age = ages, Paleolatitude = paleolatitudes, Site = site_name1)
}

final_df1 <- data.frame()


# Get a list of NGR CSV files in the directory
csv_files <- list.files("RScripts&Data/Sites Data_Age-NGR/", pattern = "\\.csv$", full.names = TRUE)

for (file in csv_files) {
  site_name1 <- tools::file_path_sans_ext(basename(file))
  if (site_name1 %in% target_sites1) {
    data <- read.csv(file, stringsAsFactors = FALSE)
    data$Site <- site_name1
    
    # Match paleolatitude for same ages
    colnames(data)[tolower(colnames(data)) == "age"] <- "Age"
    colnames(data)[tolower(colnames(data)) == "gr"] <- "GR"
    
    data$Age <- data$Age / 1000
    
    # Interpolate paleolatitude to match Age in data
    paleo_df1 <- paleolatitudes_list1[[site_name1]]
    
    # Add interpolated Paleolatitude
    data$Paleolatitude <- approx(x = paleo_df1$Age, y = paleo_df1$Paleolatitude, xout = data$Age, rule = 2)$y
    
    # Keep relevant columns
    merged_df1 <- data[, c("Site", "Age", "GR", "Paleolatitude")]
    final_df1 <- rbind(final_df1, merged_df1)
    final_df1 <- final_df1 %>%
      filter(Age >= 2.5, Age <= 21.6)
  }
}
custom_colors <- colorRampPalette(brewer.pal(10,"RdBu"))(21)

ipcc_precip <- c("#543005", "#8C510A", "#BF812D", "#F5F5F5", 
                 "#A6DBD8", "#80CDC1", "#35978F")

label_df1 <- final_df1 %>%
  group_by(Site) %>%
  slice_min(Age, n = 1, with_ties = FALSE)

setwd("RScripts&Data/Time_Slice_Map")

pdf(file = "Supplementary Figure 1.pdf", width = 11, height = 11, paper = "a4")
png(filename = "Supplementary Figure 1.png", width = 9000, height = 9000, res = 600)

ggplot(final_df1, aes(x = Paleolatitude, y = Age, color = GR)) +
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
