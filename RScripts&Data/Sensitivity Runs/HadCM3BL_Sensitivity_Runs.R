library(ncdf4)
library(ggplot2)
library(scales)
library(dplyr)
library(purrr)

timeslice_ids    <- c("tflzW", "tflzX", "tflzY")   
timeslice_labels <- c("Aquitanian", "Burdigalian", "Langhian")

co2_suffix <- c("1", "2")                
co2_labels <- c("1" = "560 ppm", "2" = "280 ppm")

precip_dir <- "RScripts&Data/Sensitivity Runs/"
evap_dir   <- "RScripts&Data/Sensitivity Runs/"
lsm_dir    <- "RScripts&Data/Sensitivity Runs/"

lon_range <- c(105, 155)
lat_range <- c(-55, -5)

get_pe_slice <- function(base_id, suffix) {
  
  id <- paste0(base_id, suffix)
  
  precip_file <- file.path(precip_dir, paste0(id, ".djf_precip_temp_data.nc"))
  evap_file   <- file.path(evap_dir,   paste0(id, ".djf_total_evap.nc"))
  lsm_file    <- file.path(lsm_dir,    paste0(id, ".qrparm.mask.nc"))
  
  nc_p   <- nc_open(precip_file)
  nc_e   <- nc_open(evap_file)
  nc_lsm <- nc_open(lsm_file)
  
  lon  <- ncvar_get(nc_p, "longitude")
  lat  <- ncvar_get(nc_p, "latitude")
  prec <- ncvar_get(nc_p, "precip_mm_srf")
  lh   <- ncvar_get(nc_e, "lh_mm_srf")
  lsm  <- ncvar_get(nc_lsm, "lsm")
  
  nc_close(nc_p)
  nc_close(nc_e)
  nc_close(nc_lsm)
  
  lat_order <- order(lat)
  lon_order <- order(lon)
  
  lat  <- lat[lat_order]
  lon  <- lon[lon_order]
  prec <- prec[lon_order, lat_order]
  lh   <- lh[lon_order, lat_order]
  lsm  <- lsm[lon_order, lat_order]
  
  lon <- seq(lon[1], lon[length(lon)], length.out = length(lon))
  lat <- seq(lat[1], lat[length(lat)], length.out = length(lat))
  
  evap    <- lh / 2.5e6
  prec_mm <- prec * 86400
  evap_mm <- evap * 86400
  pe      <- prec_mm - evap_mm
  
  lon_idx <- which(lon >= lon_range[1] & lon <= lon_range[2])
  lat_idx <- which(lat >= lat_range[1] & lat <= lat_range[2])
  
  list(lon = lon[lon_idx],
       lat = lat[lat_idx],
       pe  = pe[lon_idx, lat_idx],
       lsm = lsm[lon_idx, lat_idx])
}

combo_grid <- expand.grid(base_id = timeslice_ids, suffix = co2_suffix,
                          stringsAsFactors = FALSE)

all_data <- pmap_dfr(combo_grid, function(base_id, suffix) {
  d <- get_pe_slice(base_id, suffix)
  
  df <- expand.grid(lon = d$lon, lat = d$lat)
  df$pe        <- as.vector(d$pe)
  df$elevation <- as.vector(d$lsm)
  df$timeslice <- timeslice_labels[match(base_id, timeslice_ids)]
  df$co2       <- co2_labels[suffix]
  df
})

all_data$timeslice <- factor(all_data$timeslice, levels = timeslice_labels)
all_data$co2       <- factor(all_data$co2, levels = c("560 ppm", "280 ppm"))

ipcc_precip <- c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5",
                 "#E8F6E8", "#CDECE6", "#A6DBD8", "#80CDC1", "#35978F")

deg_lon <- function(x) paste0(abs(x), "°", ifelse(x < 0, "W", "E"))
deg_lat <- function(x) paste0(abs(x), "°", ifelse(x < 0, "S", "N"))


p <- ggplot(all_data) +
  geom_raster(aes(x = lon, y = lat, fill = pe), interpolate = TRUE) +
  geom_contour(aes(x = lon, y = lat, z = elevation),
               breaks = 0.5, color = "black", linewidth = 0.5) +
  scale_fill_gradientn(name = "P-E (mm/day)", colors = ipcc_precip,
                       limits = c(-8, 8),
                       breaks = c(-8, -4, 0, 4, 8),
                       values = scales::rescale(c(-8, -4, -2.5, -1.0, -0.5, 0, 0.5, 1.0, 2.5, 4, 8)),
                       oob = scales::squish) +
  scale_x_continuous(breaks = seq(110, 150, 40), labels = deg_lon, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-50, 0, 40),  labels = deg_lat, expand = c(0, 0)) +
  coord_sf(expand = FALSE) +
  facet_grid(co2 ~ timeslice) +
  guides(fill = guide_colorbar(title.position = "top",
                               barwidth = unit(0.5, "cm"),
                               barheight = unit(2.4, "cm"))) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16, colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        axis.text.x = element_text(margin = margin(t = 6)),
        axis.text.y = element_text(margin = margin(r = 6)),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        strip.text = element_text(size = 15, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 13, hjust = 0.5),
        legend.text  = element_text(size = 12),
        panel.spacing = unit(0.8, "lines"),
        plot.margin = margin(5, 12, 5, 5)
  )

setwd("RScripts&Data/Sensitivity Runs/")

ggsave("Extended Figure 4.png", plot = p,
       width = 14, height = 10, dpi = 600)
