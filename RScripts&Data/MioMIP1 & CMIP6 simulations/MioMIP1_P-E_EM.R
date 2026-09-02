library(ncdf4)
library(ggplot2)
library(scales)

E2MMIO400 <- "E2MMIO400-1deg-AnnualMean-PE.nc"
E2MMIO560 <- "E2MMIO560-1deg-AnnualMean-PE.nc"
E2MMIO850 <- "E2MMIO850-1deg-AnnualMean-PE.nc"


Poblete_pg <- "RScripts&Data/MioMIP1 & CMIP6 simulations/Paleogeography/EM_topo.nc"     # IPSLCM paleogeography
Farnsworth_pg <- "RScripts&Data/MioMIP1 & CMIP6 simulations/Paleogeography/tfgsY.qrparm.mask.nc"     # Langhian HadCM3 paleogeography
Burls_pg   <- "RScripts&Data/MioMIP1 & CMIP6 simulations/Paleogeography/miocene_topo_pollard_antscape_dolan_0.5x0.5.nc"     # Herold/Burls paleogeography
Zhang_pg   <- "RScripts&Data/MioMIP1 & CMIP6 simulations/Paleogeography/texpd.qrparm.mask.nc"     # NorESM paleogeography
Bradshaw_pg   <- "RScripts&Data/MioMIP1 & CMIP6 simulations/Paleogeography/teudm.qrparm.mask.nc"     # mid-Miocene HadCM3 paleogeography


inspect_exp_names <- function(file) {
  nc <- nc_open(file)
  exp_vals <- nc$dim$exp$vals
  nc_close(nc)
  print(exp_vals)
  invisible(exp_vals)
}

inspect_exp_names(E2MMIO400)
inspect_exp_names(E2MMIO560)
inspect_exp_names(E2MMIO850)


read_pe_file <- function(file) {
  nc <- nc_open(file)
  lon      <- ncvar_get(nc, "lon")
  lat      <- ncvar_get(nc, "lat")
  pe       <- ncvar_get(nc, "pe")        
  exp_vals <- as.character(nc$dim$exp$vals)
  nc_close(nc)
  
  out <- do.call(rbind, lapply(seq_along(exp_vals), function(i) {
    pe_djf_mmyr  <- (pe[,,12,i] + pe[,,1,i] + pe[,,2,i]) / 3
    pe_djf_mmday <- pe_djf_mmyr / 365
    
    df <- expand.grid(lon = lon, lat = lat)
    df$PE     <- as.vector(pe_djf_mmday)
    df$file   <- basename(file)
    df$exp_id <- exp_vals[i]
    df
  }))
  out
}

pe_400 <- read_pe_file(E2MMIO400)
pe_560 <- read_pe_file(E2MMIO560)
pe_850 <- read_pe_file(E2MMIO850)

pe_all <- rbind(pe_400, pe_560, pe_850)


exp_meta <- data.frame(
  file    = c(basename(E2MMIO400), basename(E2MMIO400), basename(E2MMIO400), basename(E2MMIO400), basename(E2MMIO400), basename(E2MMIO400),
              basename(E2MMIO400), basename(E2MMIO400), basename(E2MMIO400), basename(E2MMIO400), basename(E2MMIO400), basename(E2MMIO400),
              basename(E2MMIO560), basename(E2MMIO560), basename(E2MMIO560), basename(E2MMIO560), basename(E2MMIO560), basename(E2MMIO560),
              basename(E2MMIO560), basename(E2MMIO560), basename(E2MMIO560),
              basename(E2MMIO850), basename(E2MMIO850), basename(E2MMIO850), basename(E2MMIO850), basename(E2MMIO850)),
  exp_id  = c("CCSM-NH3 355ppm",  "CCSM3 T42 MMCO 400ppm", "CCSM3 T42 MMG 400ppm",  "CCSM4 Mid Miocene 400ppm", "CESM1 Mid Miocene 400 ppm", "COSMOS Mid Miocene 450ppm",
              "HadCM3L Mid Miocene 90SLE 400ppm",  "HadCM3L Mid Miocene 55SLE 400ppm", "HadCM3L Mid Miocene NoICE 400ppm",  "HadCM3L Langhian 400ppm", "IPSLCM 20Ma 420ppm", "NorESM-L 20Ma 350ppm",
              "CCSM-NH3 560ppm",  "CESM1 Mid Miocene 560 ppm", "HadCM3L Mid Miocene 90SLE 560ppm",  "HadCM3L Mid Miocene 55SLE 560ppm", "HadCM3L Mid Miocene NoICE 560ppm",
              "HadCM3L Langhian 560ppm", "IPSLCM 20Ma 560ppm",  "IPSLCM 20Ma NoGIS 560ppm", "NorESM-L 20Ma 560ppm",
              "CESM1 Mid Miocene 840 ppm",  "HadCM3L Mid Miocene 90SLE 850ppm", "HadCM3L Mid Miocene 55SLE 850ppm",  "HadCM3L Mid Miocene NoICE 850ppm", "IPSLCM 20Ma 840ppm"),                      
  model       = c("CCSM-NH3",  "CCSM3", "CCSM3",  "CCSM4", "CESM1", "COSMOS",
                  "HadCM3L",  "HadCM3L", "HadCM3L",  "HadCM3L", "IPSLCM", "NorESM-L",
                  "CCSM-NH3",  "CESM1", "HadCM3L",  "HadCM3L", "HadCM3L",
                  "HadCM3L", "IPSLCM",  "IPSLCM", "NorESM-L",
                  "CESM1",  "HadCM3L", "HadCM3L",  "HadCM3L", "IPSLCM"), 
  co2_ppm     = c(355, 400, 400, 400, 400, 450, 400, 400, 400, 400, 420, 350,
                   560, 560, 560, 560, 560, 560, 560, 560, 560,
                   840, 850, 850, 850, 840),
  age_Ma      = c(16, 16, 14, 16, 16, 16, 16, 16, 16, 16, 20, 20,
                  16, 16, 16, 16, 16, 16, 20, 20, 20,
                  16, 16, 16, 16, 20),
  paleo_source = c("Burls_pg", "Burls_pg", "Burls_pg", "Burls_pg", "Burls_pg", "Burls_pg", "Bradshaw_pg", "Bradshaw_pg", "Bradshaw_pg", "Farnsworth_pg", "Poblete_pg", "Zhang_pg",
                   "Burls_pg", "Burls_pg", "Bradshaw_pg", "Bradshaw_pg", "Bradshaw_pg", "Farnsworth_pg", "Poblete_pg", "Poblete_pg", "Zhang_pg",
                   "Burls_pg", "Bradshaw_pg", "Bradshaw_pg", "Bradshaw_pg", "Poblete_pg"),
  stringsAsFactors = FALSE)
  
exp_meta$facet_label <- exp_meta$exp_id
exp_meta$facet_label <- factor(exp_meta$facet_label, levels = unique(exp_meta$exp_id))

exp_meta$co2_case <- ifelse(exp_meta$file == basename(E2MMIO400), "400 ppm case",
                            ifelse(exp_meta$file == basename(E2MMIO560), "560 ppm case",
                                   ifelse(exp_meta$file == basename(E2MMIO850), "850 ppm case", NA)))

pe_sel <- merge(pe_all, exp_meta, by = c("file", "exp_id"))
facet_order <- exp_meta$facet_label[order(exp_meta$model, exp_meta$co2_ppm)]
pe_sel$facet_label <- factor(pe_sel$facet_label, levels = unique(facet_order))
pe_sel_aus <- subset(pe_sel, lon >= 105 & lon <= 155 & lat >= -55 & lat <= -5)


group_20Ma <- subset(pe_sel_aus, age_Ma == 20)
group_400  <- subset(pe_sel_aus, age_Ma %in% c(14, 16) & co2_case == "400 ppm case")
group_560_850  <- subset(pe_sel_aus, age_Ma %in% c(14, 16) & co2_case %in% c("560 ppm case", "850 ppm case"))

group_20Ma$facet_label <- droplevels(group_20Ma$facet_label)
group_400$facet_label  <- droplevels(group_400$facet_label)
group_560_850$facet_label  <- droplevels(group_560_850$facet_label)

paleo_registry <- list(
  Burls_pg = list(file = Burls_pg, lon_var = "lon", lat_var = "lat", var = "topo", contour_break = 0),
  Bradshaw_pg = list(file = Bradshaw_pg, lon_var = "longitude", lat_var = "latitude", var = "lsm", contour_break = 0.5),
  Farnsworth_pg = list(file = Farnsworth_pg, lon_var = "longitude", lat_var = "latitude", var = "lsm", contour_break = 0.5),
  Poblete_pg = list(file = Poblete_pg, lon_var = "lon", lat_var = "lat", var = "topo", contour_break = 0),
  Zhang_pg = list(file = Zhang_pg, lon_var = "longitude", lat_var = "latitude", var = "lsm", contour_break = 0.5))

load_paleo <- function(name, registry) {
  info <- registry[[name]]
  nc <- nc_open(info$file)
  lon <- ncvar_get(nc, info$lon_var)
  lat <- ncvar_get(nc, info$lat_var)
  z   <- ncvar_get(nc, info$var)
  nc_close(nc)
  
  df <- expand.grid(lon = lon, lat = lat)
  df$elevation     <- as.vector(z)
  df$paleo_source  <- name
  df$contour_break <- info$contour_break
  df}

all_coasts <- do.call(rbind, lapply(names(paleo_registry), load_paleo,
                                    registry = paleo_registry))

all_coasts_aus <- subset(all_coasts, lon >= 105 & lon <= 155 &
                           lat >= -55 & lat <= -5)

exp_meta$show_coastline <- TRUE
exp_meta$show_coastline[exp_meta$model == "COSMOS"] <- FALSE

facet_lookup <- unique(exp_meta[, c("facet_label", "paleo_source", "show_coastline")])
facet_lookup_plotted <- subset(facet_lookup, show_coastline)

no_coast_labels <- subset(facet_lookup, !show_coastline)
no_coast_labels$lon <- 110   
no_coast_labels$lat <- -10

no_coast_labels$facet_label <- factor(no_coast_labels$facet_label,
                                      levels = levels(pe_sel_aus$facet_label))

coast_expanded <- do.call(rbind, lapply(seq_len(nrow(facet_lookup_plotted)), function(i) {
  src   <- facet_lookup_plotted$paleo_source[i]
  label <- facet_lookup_plotted$facet_label[i]
  df <- subset(all_coasts_aus, paleo_source == src)
  df$facet_label <- label
  df
}))

coast_expanded$facet_label <- factor(coast_expanded$facet_label,
                                     levels = levels(pe_sel_aus$facet_label))

coast_expanded_topo <- subset(coast_expanded, contour_break == 0)
coast_expanded_lsm  <- subset(coast_expanded, contour_break == 0.5)

wrap_two_lines <- function(labels, threshold = 20) {
  sapply(labels, function(x) {
    if (nchar(x) <= threshold) return(x)   # short enough, leave as-is
    
    spaces <- gregexpr(" ", x)[[1]]
    if (spaces[1] == -1) return(x)          # no spaces to break on, leave as-is
    
    mid <- nchar(x) / 2
    break_point <- spaces[which.min(abs(spaces - mid))]  # space closest to the middle
    
    paste0(substr(x, 1, break_point - 1), "\n", substr(x, break_point + 1, nchar(x)))
  }, USE.NAMES = FALSE)
}

ipcc_precip <- c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5",
                 "#E8F6E8", "#CDECE6", "#A6DBD8", "#80CDC1", "#35978F")

deg_lon <- function(x) paste0(abs(x), "°", ifelse(x < 0, "W", "E"))
deg_lat <- function(x) paste0(abs(x), "°", ifelse(x < 0, "S", "N"))

make_pe_plot <- function(pe_data, coast_topo, coast_lsm, nrow = 1, ncol = NULL, wrap_threshold = 20) {
  
  facets_present <- levels(pe_data$facet_label)
  
  coast_topo_sub <- subset(coast_topo, facet_label %in% facets_present)
  coast_topo_sub$facet_label <- factor(coast_topo_sub$facet_label, levels = facets_present)
  
  coast_lsm_sub  <- subset(coast_lsm, facet_label %in% facets_present)
  coast_lsm_sub$facet_label <- factor(coast_lsm_sub$facet_label, levels = facets_present)
  
  ggplot(pe_data) +
    geom_raster(aes(x = lon, y = lat, fill = PE), interpolate = TRUE) +
    geom_contour(data = coast_topo_sub,
                 aes(x = lon, y = lat, z = elevation, group = facet_label),
                 breaks = 0, color = "black", linewidth = 0.5) +
    geom_contour(data = coast_lsm_sub,
                 aes(x = lon, y = lat, z = elevation, group = facet_label),
                 breaks = 0.5, color = "black", linewidth = 0.5) +
    facet_wrap(~ facet_label, nrow = nrow, ncol = ncol,
               labeller = labeller(facet_label = function(x) wrap_two_lines(x, threshold = wrap_threshold))) +
    scale_fill_gradientn(
      name = "P-E (mm/day)", colors = ipcc_precip,
      limits = c(-8, 8),
      values = scales::rescale(c(-8, -4, -2.5, -1.0, -0.5, 0, 0.5, 1.0, 2.5, 4, 8)),
      oob = scales::squish) +
    scale_x_continuous(breaks = seq(110, 150, 40), labels = deg_lon, expand = c(0,0)) +
    scale_y_continuous(breaks = seq(-50, 0, 40), labels = deg_lat, expand = c(0,0)) +
    coord_sf(xlim = c(105, 155), ylim = c(-55, -5), expand = FALSE) +
    theme_minimal(base_size = 14) +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 11, colour = "black"),
          axis.text.x = element_text(size = 10),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          panel.spacing.x = unit(1, "lines"),
          panel.spacing.y = unit(1.2, "lines"),
          strip.text = element_text(size = 10, face = "bold", margin = margin(b = 4, t = 4)),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.length = unit(-0.25, "cm"),
          legend.title = element_text(size = 12),
          legend.text  = element_text(size = 11))
}

p_20Ma <- make_pe_plot(group_20Ma, coast_expanded_topo, coast_expanded_lsm, nrow = 2, ncol = 3, wrap_threshold = 20)
p_400  <- make_pe_plot(group_400,  coast_expanded_topo, coast_expanded_lsm, nrow = 3, ncol = 4, wrap_threshold = 16)
p_560_850  <- make_pe_plot(group_560_850,  coast_expanded_topo, coast_expanded_lsm, nrow = 3, ncol = 4, wrap_threshold = 16)

add_no_coast_label <- function(plot, pe_data, no_coast_labels) {
  facets_present <- levels(pe_data$facet_label)
  labels_sub <- subset(no_coast_labels, facet_label %in% facets_present)
  labels_sub$facet_label <- factor(labels_sub$facet_label, levels = facets_present)
  
  plot +
    geom_text(data = labels_sub,
              aes(x = lon, y = lat, label = "coastline\nunavailable"),
              inherit.aes = FALSE, size = 2.5, fontface = "italic",
              color = "grey30", hjust = 0)
}

p_20Ma <- add_no_coast_label(p_20Ma, group_20Ma, no_coast_labels)
p_400  <- add_no_coast_label(p_400,  group_400,  no_coast_labels)
p_560_850  <- add_no_coast_label(p_560_850,  group_560_850,  no_coast_labels)

save_pe_plot <- function(plot, filename, nrow, ncol, panel_width = 1700, panel_height = 1900, res = 600) {
  width  <- ncol * panel_width
  height <- nrow * panel_height
  png(filename = filename, width = width, height = height, res = res)
  print(plot)
  dev.off()
}

setwd("RScripts&Data/MioMIP1 & CMIP6 simulations/")
save_pe_plot(p_20Ma, "Extended Figure 3.png", nrow = 2, ncol = 3)
save_pe_plot(p_400,  "Extended Figure 4.png", nrow = 3, ncol = 4)
save_pe_plot(p_560_850,  "Extended Figure 5.png", nrow = 3, ncol = 4)


