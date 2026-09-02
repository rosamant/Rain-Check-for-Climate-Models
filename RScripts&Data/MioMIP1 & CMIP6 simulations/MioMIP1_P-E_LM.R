library(ncdf4)
library(ggplot2)
library(scales)

M2LMIO280 <- "M2LMIO280-1deg-AnnualMean-PE.nc"
M2LMIO400 <- "M2LMIO400-1deg-AnnualMean-PE.nc"
M2LMIO560 <- "M2LMIO560-1deg-AnnualMean-PE.nc"


Poblete_pg <- "RScripts&Data/MioMIP1 & CMIP6 simulations/Paleogeography/LM_topo.nc"     # IPSLCM paleogeography
Farnsworth_pg <- "RScripts&Data/MioMIP1 & CMIP6 simulations/Paleogeography/tfgSa.qrparm.mask.nc"     # Langhian HadCM3 paleogeography
Farnsworth_pg_mes <- "RScripts&Data/MioMIP1 & CMIP6 simulations/Paleogeography/tfgSb.qrparm.mask.nc"     # Langhian HadCM3 paleogeography
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

inspect_exp_names(M2LMIO280)
inspect_exp_names(M2LMIO400)
inspect_exp_names(M2LMIO560)


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

pe_280 <- read_pe_file(M2LMIO280)
pe_400 <- read_pe_file(M2LMIO400)
pe_560 <- read_pe_file(M2LMIO560)

pe_all <- rbind(pe_280, pe_400, pe_560)


exp_meta <- data.frame(
  file    = c(basename(M2LMIO280), basename(M2LMIO280), basename(M2LMIO280), basename(M2LMIO280),
              basename(M2LMIO400), basename(M2LMIO400), basename(M2LMIO400), basename(M2LMIO400), basename(M2LMIO400), basename(M2LMIO400), basename(M2LMIO400),
              basename(M2LMIO560), basename(M2LMIO560)),
  exp_id  = c("COSMOS Late Miocene 278ppm",  "COSMOS Late Miocene ICEQ 278ppm", "HadCM3L Late Miocene 280ppm",  "IPSLCM Late Miocene 300ppm",
              "COSMOS Late Miocene 450ppm",  "COSMOS Late Miocene ICEQ 450ppm", "HadCM3L Late Miocene 400ppm",  "HadCM3L Tortonian 400ppm", "HadCM3L Messinian 400ppm", "NorESM-L 10Ma 350ppm", "IPSLCM Late Miocene 420ppm",  
              "NorESM-L 10Ma 560ppm", "IPSLCM Late Miocene 560ppm"),                      
  model       = c("COSMOS",  "COSMOS", "HadCM3L", "IPSLCM",
                  "COSMOS",  "COSMOS", "HadCM3L",  "HadCM3L", "HadCM3L", "NorESM-L",  "IPSLCM",
                  "NorESM-L",  "IPSLCM"), 
  co2_ppm     = c(278, 278, 280, 300, 
                  450, 450, 400, 400, 400, 350, 420,
                  560, 560),
  age_Ma      = c(10, 10, 10, 10, 
                  10, 10, 10, 10, 6, 10, 10,
                  10, 10),
  paleo_source = c("Burls_pg", "Burls_pg", "Bradshaw_pg", "Poblete_pg",
                   "Burls_pg", "Burls_pg", "Bradshaw_pg", "Farnsworth_pg", "Farnsworth_pg_mes", "Zhang_pg", "Poblete_pg",
                   "Zhang_pg", "Poblete_pg"),
  stringsAsFactors = FALSE)

exp_meta$facet_label <- exp_meta$exp_id
exp_meta$facet_label <- factor(exp_meta$facet_label, levels = unique(exp_meta$exp_id))

exp_meta$co2_case <- ifelse(exp_meta$file == basename(M2LMIO280), "280 ppm case",
                            ifelse(exp_meta$file == basename(M2LMIO400), "400 ppm case",
                                   ifelse(exp_meta$file == basename(M2LMIO560), "560 ppm case", NA)))

pe_sel <- merge(pe_all, exp_meta, by = c("file", "exp_id"))
facet_order <- exp_meta$facet_label[order(exp_meta$model, exp_meta$co2_ppm)]
pe_sel$facet_label <- factor(pe_sel$facet_label, levels = unique(facet_order))
pe_sel_aus <- subset(pe_sel, lon >= 105 & lon <= 155 & lat >= -55 & lat <= -5)

group_280_560 <- subset(pe_sel_aus, age_Ma %in% c(6, 10) & co2_case %in% c("280 ppm case", "560 ppm case"))
group_400  <- subset(pe_sel_aus, age_Ma %in% c(6, 10) & co2_case == "400 ppm case")

group_400$facet_label  <- droplevels(group_400$facet_label)
group_280_560$facet_label <- droplevels(group_280_560$facet_label)

paleo_registry <- list(
  Burls_pg = list(file = Burls_pg, lon_var = "lon", lat_var = "lat", var = "topo", contour_break = 0),
  Bradshaw_pg = list(file = Bradshaw_pg, lon_var = "longitude", lat_var = "latitude", var = "lsm", contour_break = 0.5),
  Farnsworth_pg = list(file = Farnsworth_pg, lon_var = "longitude", lat_var = "latitude", var = "lsm", contour_break = 0.5),
  Farnsworth_pg_mes = list(file = Farnsworth_pg_mes, lon_var = "longitude", lat_var = "latitude", var = "lsm", contour_break = 0.5),
  Poblete_pg = list(file = Poblete_pg, lon_var = "longitude", lat_var = "latitude", var = "topo", contour_break = 0),
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
    if (nchar(x) <= threshold) return(x)   
    
    spaces <- gregexpr(" ", x)[[1]]
    if (spaces[1] == -1) return(x)          
    
    mid <- nchar(x) / 2
    break_point <- spaces[which.min(abs(spaces - mid))]  
    
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

p_280_560  <- make_pe_plot(group_280_560,  coast_expanded_topo, coast_expanded_lsm, nrow = 2, ncol = 3, wrap_threshold = 20)
p_400  <- make_pe_plot(group_400,  coast_expanded_topo, coast_expanded_lsm, nrow = 3, ncol = 4, wrap_threshold = 16)

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

p_280_560  <- add_no_coast_label(p_280_560,  group_280_560,  no_coast_labels)
p_400  <- add_no_coast_label(p_400,  group_400,  no_coast_labels)

save_pe_plot <- function(plot, filename, nrow, ncol, panel_width = 1700, panel_height = 1900, res = 600) {
  width  <- ncol * panel_width
  height <- nrow * panel_height
  png(filename = filename, width = width, height = height, res = res)
  print(plot)
  dev.off()
}

setwd("RScripts&Data/MioMIP1 & CMIP6 simulations/")
save_pe_plot(p_280_560,  "Extended Figure 6.png", nrow = 2, ncol = 3)
save_pe_plot(p_400,  "Extended Figure 7.png", nrow = 2, ncol = 4)

