# Load required libraries
library(ncdf4)
library(ggplot2)
library(scales)

Lan_TR <- nc_open("RScripts&Data/NetCDF/tfgsY.djf_runoff_data.nc")
Lan_LSM <- nc_open("RScripts&Data/NetCDF/tfgsY.qrparm.mask.nc")
names(Lan_TR$var)
names(Lan_LSM$var)

lon_vec <- ncvar_get(Lan_TR, "longitude")
lat_vec <- ncvar_get(Lan_TR, "latitude")
Lan_SR <- ncvar_get(Lan_TR, "srfRunoff_mm_srf")
Lan_SSR <- ncvar_get(Lan_TR, "subsrfRunoff_mm_srf")
Lan_R <- Lan_SSR + Lan_SR
Lan_lsm <- ncvar_get(Lan_LSM, "lsm")

nc_close(Lan_TR)
nc_close(Lan_LSM)

Lan_R_mm <- Lan_R * 86400

Lan <- expand.grid(lon = lon_vec, lat = lat_vec)
Lan$elevation <- as.vector(Lan_lsm)
Lan$R <- as.vector(Lan_R_mm)

Lan_aus <- subset(Lan, lon >= 105 & lon <= 155 & lat >= -55 & lat <= -5)
deg_lon <- function(x) paste0(abs(x), "°", ifelse(x < 0, "W", "E"))
deg_lat <- function(x) paste0(abs(x), "°", ifelse(x < 0, "S", "N"))

ipcc_r <- c("#EDF7F2",  "#D0EADF",  "#A9DDD2",  "#79C9C0",  "#49B0AE",  "#217F88",  "#0B5F68")

ipcc_precip <- c("#BF812D", "#DFC27D", "#F5F5F5", "#80CDC1", "#35978F", "#003C30")

poly_coords <- data.frame(lon = c(110, 139, 146, 116, 110),
                          lat = c(-30, -9, -9, -30, -30))

vals <- c(0.01, 0.1, 0.2, 0.35, 0.5, 1, 2, 4)

# Total Runoff 

setwd("RScripts&Data/NetCDF/")

png(filename = "Langhian_R_DJF.png", width = 5000, height = 5000, res = 600)

ggplot(Lan_aus) +
  geom_raster(aes(x = lon, y = lat, fill = R), interpolate = TRUE) +
  geom_contour(aes(x = lon, y = lat, z = elevation), 
               breaks = 0.5, color = "black", linewidth = 0.5) +
  scale_fill_stepsn(name = "Total Runoff\n(mm/day)", colors = ipcc_precip,
                    limits = c(0.01, 4), trans = "log10", breaks = c(0.1, 0.2, 0.35, 0.5, 1, 2), 
                    labels = c(0.1, 0.2, 0.35, 0.5, 1, 2), oob = scales::squish,
                    guide = guide_colorsteps(barheight = unit(6, "cm"))) +
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

