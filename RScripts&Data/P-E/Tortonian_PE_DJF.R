# Load required libraries
library(ncdf4)
library(ggplot2)
library(scales)

Tor_PE <- nc_open("RScripts&Data/NetCDF/tfgSa.djf_precip_evap_data.nc")
Tor_LSM <- nc_open("RScripts&Data/NetCDF/tfgSa.qrparm.mask.nc")
names(Tor_PE$var)
names(Tor_LSM$var)

lon_vec <- ncvar_get(Tor_PE, "longitude")
lat_vec <- ncvar_get(Tor_PE, "latitude")
Tor_prec <- ncvar_get(Tor_PE, "precip_mm_srf")
Tor_lh <- ncvar_get(Tor_PE, "lh_mm_srf")
Tor_lsm <- ncvar_get(Tor_LSM, "lsm")

nc_close(Tor_PE)
nc_close(Tor_LSM)

Tor_evap <- Tor_lh / 2.5e6
Tor_pe <- (Tor_prec - Tor_evap) * 86400
Tor_prec_mm <- Tor_prec * 86400

Tor <- expand.grid(lon = lon_vec, lat = lat_vec)
Tor$elevation <- as.vector(Tor_lsm)
Tor$pe <- as.vector(Tor_pe)
Tor$prec <- as.vector(Tor_prec_mm)

Tor_aus <- subset(Tor, lon >= 105 & lon <= 155 & lat >= -55 & lat <= -5)
deg_lon <- function(x) paste0(abs(x), "°", ifelse(x < 0, "W", "E"))
deg_lat <- function(x) paste0(abs(x), "°", ifelse(x < 0, "S", "N"))


ipcc_p <- c("#F5F5F5", "#E8F6E8", "#CDECE6", "#A6DBD8", "#80CDC1", "#35978F")

ipcc_precip <- c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", 
                 "#E8F6E8", "#CDECE6", "#A6DBD8", "#80CDC1", "#35978F")

poly_coords <- data.frame(lon = c(110, 130, 140, 120, 110),
                          lat = c(-30, -9, -9, -30, -30))

# Precipitation 

setwd("RScripts&Data/P-E/")

png(filename = "Tortonian_P_DJF.png", width = 5000, height = 5000, res = 600)

ggplot(Tor_aus) +
  geom_raster(aes(x = lon, y = lat, fill = prec), interpolate = TRUE) +
  geom_contour(aes(x = lon, y = lat, z = elevation), 
               breaks = 0.5, color = "black", linewidth = 0.5) +
  scale_fill_gradientn(name = "Precipitation\n(mm/day)", colors = ipcc_p,
                       limits = c(0, 8), values = scales::rescale(c(0, 0.5, 1.5, 3, 5, 8)), oob = scales::squish)+
  scale_x_continuous(breaks = seq(110,150,40), labels = deg_lon, expand = c(0,0))+
  scale_y_continuous(breaks = seq(-50,0,40), labels = deg_lat, expand = c(0,0))+
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_blank(),axis.text = element_text(size = 16, colour = "black"),
        panel.grid = element_blank(),panel.background = element_rect(fill = "white",colour = NA), 
        axis.text.x = element_text(margin = margin(t = 6)),
        axis.text.y = element_text(margin = margin(r = 6)),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 11)
  )
dev.off()


# Precipitation minus Evaporation

setwd("RScripts&Data/P-E/")

png(filename = "Tortonian_PE_DJF.png", width = 5000, height = 5000, res = 600)

ggplot(Tor_aus) +
  geom_raster(aes(x = lon, y = lat, fill = pe), interpolate = TRUE) +
  geom_polygon(data = poly_coords, aes(x = lon, y = lat), 
               fill = NA, color = "darkred", linewidth = 1, linetype = 2) +
  geom_contour(aes(x = lon, y = lat, z = elevation), 
               breaks = 0.5, color = "black", linewidth = 0.5) +
  scale_fill_gradientn(name = "P-E (mm/day)", colors = ipcc_precip,
                       limits = c(-8, 8), values = scales::rescale(c(-8, -4, -2.5, -1.0, -0.5, 0, 0.5, 1.0, 2.5, 4, 8)), oob = scales::squish)+
  scale_x_continuous(breaks = seq(110,150,40), labels = deg_lon, expand = c(0,0))+
  scale_y_continuous(breaks = seq(-50,0,40), labels = deg_lat, expand = c(0,0))+
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_blank(),axis.text = element_text(size = 22, colour = "black"),
        panel.grid = element_blank(),panel.background = element_rect(fill = "white",colour = NA), 
        axis.text.x = element_text(margin = margin(t = 6)),
        axis.text.y = element_text(margin = margin(r = 6)),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(-0.25, "cm"),
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 11)
  )
dev.off()