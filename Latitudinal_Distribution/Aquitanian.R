#### NGR Data ####
library(rgplates)
library(dplyr)                  
library(ggplot2)               
library(ncdf4)
library(scales)
library(sf)

## Location of the sites ##
Sites2 <- read.csv("RScripts&Data/Sites.csv", header=TRUE, stringsAsFactors=FALSE)
Sites2 <- Sites2[c("Site", "Longitude", "Latitude")]

target_sites2 <- c("Mandorah 1", "Hadrian 1", "Calliance 2", "Omar 1", "Ermine 1", "Brigadier 1", "Eastbrook 1", "Central Gorgon 1", "Bluebell 1", "York 1")

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

#### Model Data ####

Aqui_PE <- nc_open("RScripts&Data/NetCDF/tflzW.djf_precip_evap_data.nc")
names(Aqui_PE$var)

lon_vec <- ncvar_get(Aqui_PE, "longitude")
lat_vec <- ncvar_get(Aqui_PE, "latitude")
Aqui_prec <- ncvar_get(Aqui_PE, "precip_mm_srf")
Aqui_lh <- ncvar_get(Aqui_PE, "lh_mm_srf")

nc_close(Aqui_PE)

Aqui_evap <- Aqui_lh / 2.5e6
Aqui_pe <- (Aqui_prec - Aqui_evap) * 86400

Aqui <- expand.grid(lon = lon_vec, lat = lat_vec)
Aqui$pe <- as.vector(Aqui_pe)
Aqui_sf <- st_as_sf(Aqui, coords = c("lon", "lat"), crs = 4326)

# Masking

coords <- matrix(c(110, -30, 130, -9, 140, -9, 120, -30, 110, -30),
                 ncol = 2,  byrow = TRUE)

poly <- st_polygon(list(coords))
mask_sf <- st_sfc(poly, crs = 4326)
Aqui_mask <- Aqui_sf[st_within(Aqui_sf, mask_sf, sparse = FALSE), ]


#### Data slicing ####

lat_breaks <- seq(-30, -9, by = 1)
lat_mids   <- head(lat_breaks, -1) + 0.5

# NGR Hist
ngr_Aqui_hist <- final_df2 %>%
  filter(Age >= 20.8, Age <= 22.8) %>%
  mutate(lat_bin = cut(Paleolatitude, breaks = lat_breaks, include.lowest = TRUE, labels = lat_mids)) %>%
  group_by(lat_bin) %>%
  summarise(median_GR = median(GR, na.rm = TRUE), Q1 = quantile(GR, 0.25, na.rm = TRUE),
            Q3 = quantile(GR, 0.75, na.rm = TRUE), site_count = n_distinct(Site), .groups = 'drop') %>%
  mutate(lat = as.numeric(as.character(lat_bin)), obs_scaled = (median_GR / 100),
         Q1_scaled = Q1 / 100, Q3_scaled = Q3 / 100)

# Model Hist

Aqui_pe_min <- min(Aqui_mask$pe, na.rm = TRUE)
Aqui_pe_max <- max(Aqui_mask$pe, na.rm = TRUE)

pe_Aqui_hist <- Aqui_mask %>%
  mutate(lat_val = st_coordinates(.)[,2]) %>%
  mutate(lat_bin = cut(lat_val, breaks = lat_breaks, include.lowest = TRUE, labels = lat_mids)) %>%
  group_by(lat_bin) %>%
  summarise(model_raw = median(pe, na.rm = TRUE), .groups = 'drop') %>%
  mutate(lat = as.numeric(as.character(lat_bin)),
         model_scaled = (model_raw - Aqui_pe_min) / (Aqui_pe_max - Aqui_pe_min))

common_Aqui_lats <- seq(min(ngr_Aqui_hist$lat), max(ngr_Aqui_hist$lat), by = 1)
Aqui_interp_obs <- approx(ngr_Aqui_hist$lat, ngr_Aqui_hist$obs_scaled, xout = common_Aqui_lats, rule = 2)$y
Aqui_interp_model <- approx(pe_Aqui_hist$lat, pe_Aqui_hist$model_scaled, xout = common_Aqui_lats, rule = 2)$y


#### Plotting ####

setwd("RScripts&Data/Latitudinal Distribution")

png(filename = "Aquitanian_PE_DJF.png", width = 6000, height = 4500, res = 600)

par(mar = c(4.5, 4.5, 1.5, 1.5), cex.lab = 1.5, cex.axis = 1.5) 

plot(pe_Aqui_hist$lat, pe_Aqui_hist$model_scaled, type = "b", pch = 16, col = "gray60", lty = 1, lwd = 2,
     ylim = c(0, 1), xlim = c(-30,-10), xaxt = "n", yaxt = "n",
     bty = "n", xlab = "Paleolatitude", ylab = "")
axis(1, at = seq(-30, -10, by = 10), labels = paste0(abs(seq(-30, -10, by = 10)), "°S"),
     lwd = 1, lwd.ticks = 1)
axis(2, at = seq(0, 1, by = 0.5), labels = paste0(abs(seq(0, 100, by = 50)), "%"),  lwd = 1, las = 1, lwd.ticks = 1)
polygon(x = c(common_Aqui_lats, rev(common_Aqui_lats)), y = c(Aqui_interp_model, rev(Aqui_interp_obs)), density = 20, angle = 45,
  col = "#F2B32C", border = NA)
lines(pe_Aqui_hist$lat, pe_Aqui_hist$model_scaled, col = "gray60", lty = 1, lwd = 2)
points(ngr_Aqui_hist$lat, ngr_Aqui_hist$obs_scaled, pch = 18, col = "black", cex = 1.5)
lines(ngr_Aqui_hist$lat, ngr_Aqui_hist$obs_scaled, col = "black", lwd = 2, lty = 2)
arrows(x0 = ngr_Aqui_hist$lat, y0 = ngr_Aqui_hist$Q1_scaled, x1 = ngr_Aqui_hist$lat,
  y1 = ngr_Aqui_hist$Q3_scaled, angle = 90, code = 3, length = 0.04, col = "black", lwd = 1.2)

dev.off()
