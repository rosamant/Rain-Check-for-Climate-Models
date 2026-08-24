# Load required libraries
library(ncdf4)
library(ggplot2)
library(scales)

Ser_PE <- nc_open("RScripts&Data/NetCDF/tfgsZ.djf_precip_evap_data.nc")
Ser_LSM <- nc_open("RScripts&Data/NetCDF/tfgsZ.qrparm.mask.nc")
names(Ser_PE$var)
names(Ser_LSM$var)

lon_vec <- ncvar_get(Ser_PE, "longitude")
lat_vec <- ncvar_get(Ser_PE, "latitude")
Ser_prec <- ncvar_get(Ser_PE, "precip_mm_srf")
Ser_lh <- ncvar_get(Ser_PE, "lh_mm_srf")
Ser_lsm <- ncvar_get(Ser_LSM, "lsm")

nc_close(Ser_PE)
nc_close(Ser_LSM)

Ser_evap <- Ser_lh / 2.5e6
Ser_pe <- (Ser_prec - Ser_evap) * 86400
Ser_prec_mm <- Ser_prec * 86400

Ser <- expand.grid(lon = lon_vec, lat = lat_vec)
Ser$elevation <- as.vector(Ser_lsm)
Ser$pe <- as.vector(Ser_pe)
Ser$prec <- as.vector(Ser_prec_mm)

Ser_aus <- subset(Ser, lon >= 105 & lon <= 155 & lat >= -55 & lat <= -5)
deg_lon <- function(x) paste0(abs(x), "°", ifelse(x < 0, "W", "E"))
deg_lat <- function(x) paste0(abs(x), "°", ifelse(x < 0, "S", "N"))


ipcc_p <- c("#F5F5F5", "#E8F6E8", "#CDECE6", "#A6DBD8", "#80CDC1", "#35978F")

ipcc_precip <- c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", 
                 "#E8F6E8", "#CDECE6", "#A6DBD8", "#80CDC1", "#35978F")

poly_coords <- data.frame(lon = c(110, 130, 140, 120, 110),
                          lat = c(-30, -9, -9, -30, -30))

# Precipitation 

setwd("RScripts&Data/NetCDF/")

png(filename = "Serravallian_P_DJF.png", width = 5000, height = 5000, res = 600)

ggplot(Ser_aus) +
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

setwd("RScripts&Data/NetCDF/")

png(filename = "Serravallian_PE_DJF.png", width = 5000, height = 5000, res = 600)

ggplot(Ser_aus) +
  geom_raster(aes(x = lon, y = lat, fill = pe), interpolate = TRUE) +
  geom_polygon(data = poly_coords, aes(x = lon, y = lat), 
               fill = NA, color = "darkred", linewidth = 1, linetype = 2) +
  geom_contour(aes(x = lon, y = lat, z = elevation), 
               breaks = 0.5, color = "black", linewidth = 0.5) +
  #geom_contour(aes(lon, lat, z = pe), breaks = 4, color = "black", linewidth = 0.6, linetype = "dashed") +
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