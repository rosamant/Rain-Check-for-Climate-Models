#### Import libraries ####

library(ggplot2)                
library(sf)                     
library(sp)                     
library(dplyr)                  
library(RColorBrewer)           
library(raster)                 
library(terra)
library(scales)
library(readxl)
library(stringr)
library(viridis)
library(png)
library(magick)
library(grid)
library(patchwork)
library(rgplates)
theme_set(theme_bw())           

#### Bathymetry Data ####

gebco_raster <- rast("RScripts&Data/Time_Slice_Map/gebco_australia_bathymetry.tif")
gebco_crop <- crop(gebco_raster, extent(112, 128, -23, -10))
gebco_crop <- aggregate(gebco_crop, fact = 10)  
gebco_df <- as.data.frame(gebco_crop, xy = TRUE)
colnames(gebco_df) <- c("x", "y", "elevation")
gebco_df$elevation <- pmin(pmax(gebco_df$elevation, -6000), 2000)

#### NGR Data from all Sites ####

# Get a list of NGR CSV files in the directory
csv_files <- list.files("RScripts&Data/Sites Data_Age-NGR/", pattern = "\\.csv$", full.names = TRUE)

# Create an empty list to store data frames
data_list <- list()

# Iterate over each CSV file
for (file in csv_files) {
  # Read data from CSV file
  data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  
  # Store data frame in the list
  data_list[[file]] <- data
  
  # Retain the 'Site' variable
  data$Site <- sub("\\.csv$", "", basename(file))  # Remove .csv suffix from file name
  
  # Create time slices of 100 years
  time_slices <- lapply(seq(2500, max(data$AGE), by = 100), function(year) {
    subset_data <- subset(data, AGE >= year & AGE < year + 100)
    if (nrow(subset_data) > 0) {
      subset_data$TimeSlice <- paste(year + 100, "-", year, sep = "")
      return(subset_data)
    } else {
      return(NULL)
    }
  })
  
  # Store time slices in the list
  data_list[[file]] <- do.call(rbind, time_slices)
}

# Combine all time slices from all datasets into a single data frame
combined_data <- do.call(rbind, data_list)

# Calculate mean GR values for every site for each time slice
mean_gr <- aggregate(GR ~ Site + TimeSlice, data = combined_data, FUN = mean)

# Extract the numerical part from the 'TimeSlice' column
mean_gr <- mean_gr %>%  mutate(TimeSlice_numeric = as.numeric(sub("([0-9]+)-[0-9]+", "\\1", TimeSlice)))

# Arrange the dataframe by the extracted numerical part of the 'TimeSlice' column in increasing order
mean_gr <- mean_gr %>% arrange(TimeSlice_numeric)

# Ensure mean_gr is a data frame
mean_gr <- as.data.frame(mean_gr)

# Remove the temporary 'TimeSlice_numeric' column
mean_gr <- subset(mean_gr, select = -TimeSlice_numeric)

# Location of the sites
Sites <- read.csv("RScripts&Data/Sites.csv", header=TRUE, stringsAsFactors=FALSE)

# Create a spatial points data frame
sitepoints_sp <- SpatialPointsDataFrame(coords = Sites[, c("Longitude", "Latitude")],
                                        data = Sites,
                                        proj4string = CRS("+proj=longlat +datum=WGS84"))

# Convert to sf object
sitepoints <- st_as_sf(sitepoints_sp)

#### Calculate paleolatitude and paleolatitude for time slice map ####

lat_refs <- data.frame(Site = c("lat_7S", "lat_10S", "lat_13S", "lat_16S", "lat_19S", "lat_22S", "lat_25S"),
                       Longitude = rep(120, 7), Latitude  = seq(-7,-25,-3))

paleolatitudes_list <- list()
for (i in 1:nrow(lat_refs)) {
  site_name <- lat_refs$Site[i]
  coord <- lat_refs[i, 2:3]
  ages <- seq(0, 22, 1)
  paleo <- reconstruct(coord, age = ages)
  paleolatitudes <- sapply(paleo, function(x) x[2])
  paleolatitudes_list[[site_name]] <- data.frame(
    Age           = ages,
    Paleolatitude = paleolatitudes,
    Site          = site_name
  )
}

paleo_df <- do.call(rbind, paleolatitudes_list)


lon_refs <- data.frame(Site = c("lon_108E","lon_112E", "lon_116E", "lon_120E", "lon_124E", "lon_128E", "lon_132E"),
                       Longitude = c(108, 112, 116, 120, 124, 128, 132), Latitude  = rep(-16, 7))

paleolongitudes_list <- list()
for (i in 1:nrow(lon_refs)) {
  site_name <- lon_refs$Site[i]
  coord     <- lon_refs[i, 2:3]
  ages      <- seq(0, 22, 1)
  paleo     <- reconstruct(coord, age = ages)
  paleolons <- sapply(paleo, function(x) x[1])
  paleolongitudes_list[[site_name]] <- data.frame(
    Age          = ages,
    PaleoLon     = paleolons,
    Site         = site_name
  )
}

paleo_lon_df <- do.call(rbind, paleolongitudes_list)

#### Extracting Reef Areal Extent ####

# Reefs Shapefile
Reefs_shapefiles <- list(
  list(file = "RScripts&Data/Reef Shp Files/Late Burdigallian.shp", from = 18600, to = 15700),
  list(file = "RScripts&Data/Reef Shp Files/Langhian.shp", from = 15600, to = 12100),
  list(file = "RScripts&Data/Reef Shp Files/Serravallian.shp", from = 12000, to = 10800),
  list(file = "RScripts&Data/Reef Shp Files/Tortonian.shp", from = 10700, to = 7300),
  list(file = "RScripts&Data/Reef Shp Files/Messinian_Tortoniant.shp", from = 7200, to = 6400),
  list(file = "RScripts&Data/Reef Shp Files/Messinian.shp", from = 6300, to = 5900),
  list(file = "RScripts&Data/Reef Shp Files/Pliocene.shp", from = 5800, to = 2500)
)

# Converting Shapefile to Spatial Feature
Reefs_shapefiles_age <- lapply(Reefs_shapefiles, function(info) {
  sf_obj <- st_read(info$file, quiet = TRUE)
  sf_obj$from_ka <- info$from
  sf_obj$to_ka <- info$to
  return(sf_obj)
})

Reefs_shapefiles_age_all <- bind_rows(Reefs_shapefiles_age)
Reefs_shapefiles_age_all <- st_transform(Reefs_shapefiles_age_all, crs = 4326)


# Reef TIF files

Reef_rasters <- list(
  list(file = "RScripts&Data/Reef Shp Files/Early-Middle_Miocene_Reefs_modified.tif", from = 16300, to = 10500),
  list(file = "RScripts&Data/Reef Shp Files/Late_Miocene_Reefs_modified.tif", from = 10400, to = 5430),
  list(file = "RScripts&Data/Reef Shp Files/Pliocene_Reefs_modified.tif", from = 5330, to = 2500)
)

Reef_raster_sf <- lapply(Reef_rasters, function(info) {
  r <- rast(info$file)
  alpha <- r[[4]]
  alpha[alpha == 0] <- NA 
  poly <- as.polygons(alpha, dissolve = TRUE) |> st_as_sf() |> st_make_valid()
  poly$from_ka <- info$from
  poly$to_ka <- info$to
  return(poly)
})

Reef_raster_sf_all <- do.call(bind_rows, Reef_raster_sf)
Reef_raster_sf_all <- st_transform(Reef_raster_sf_all, crs = 4326)
All_reef_extents <- bind_rows(Reefs_shapefiles_age_all, Reef_raster_sf_all)
reef_legend_df <- data.frame(x = 0, y = 0, label = "Reef extent")

#### Mean NGR Image ####

setwd("RScripts&Data/Time_Slice_Map")
img <- image_read("Time_Slice_Map_Mean_GR_CENO.png")

ngr_img <- image_resize(img, "2500x")
ngr_img <- image_trim(ngr_img)
ngr_img <- as.raster(ngr_img)
img_h <- nrow(ngr_img)
img_w <- ncol(ngr_img)
ngr_aspect <- img_h / img_w

make_ngr_panel <- function(current_age){
  
  # calibration values
  age_top <- 2
  age_bottom <- 22
  
  column_top <- 0.042
  column_bottom <- 0.92
  
  # normalized y position
  ypos <- column_top + (current_age - age_top) *
    (column_bottom - column_top) / (age_bottom - age_top)
  ypos <- 1 - ypos
  
  ggplot() + annotation_custom(rasterGrob(ngr_img, interpolate = TRUE,
                                          width  = unit(1, "npc"), height = unit(1, "npc")),
                               xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    geom_rect(aes(xmin = 0.07, xmax = 0.97, ymin = ypos - 0.005,
                  ymax = ypos + 0.005), fill = "#FF6B6B", alpha = 0.75 ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()+
    theme(
      plot.margin  = margin(0, 0, 0, 0),
      panel.spacing = unit(0, "null"),
      plot.background  = element_blank(),
      panel.background = element_blank(),
      aspect.ratio     = ngr_aspect
    )
}

#### Color Bar ####

topo_colors <- c("#303130","#444645","#606361","#707371","#999c9a","#aaacaa","#babcbb", "#d3d4d3",
                 "#e5dcc9","#f3e3c2","#edd096","#edc980","#eabf6a","#deaa41","#ac8332","#825e17")

elev_breaks <- c(-7000, -6000, -4000, -2000, -1000, -400, -150, 0,
                 20, 100, 200, 300, 400, 550, 750, 1000)

scaled_breaks <- (elev_breaks - min(elev_breaks)) / (max(elev_breaks) - min(elev_breaks))


#### Plotting ####

# Create a ggplot object for each time slice
time_slice_plots <- lapply(unique(mean_gr$TimeSlice)[length(unique(mean_gr$TimeSlice)):1], function(ts) {
  # Subset the data for the current time slice
  ts_data <- subset(mean_gr, TimeSlice == ts)
  
  # Merge latitude and longitude data with GR values
  ts_data <- merge(ts_data, sitepoints, by = "Site")
  
  # Create an sf object for the merged data
  ts_data_sf <- st_as_sf(ts_data, coords = c("Longitude", "Latitude"))
  
  # Assign CRS to ts_data_sf
  ts_data_sf <- st_set_crs(ts_data_sf, 4326)
  ts_vals <- as.numeric(unlist(strsplit(ts, "-")))
  ts_start <- max(ts_vals)  
  ts_end <- min(ts_vals)    
  ts_start_ma <- ts_start / 1000
  ts_end_ma   <- ts_end / 1000
  ts_title <- sprintf("%.2f–%.2f Ma", ts_start_ma, ts_end_ma)
  
  current_age <- mean(c(ts_start_ma, ts_end_ma))
  
  modern_breaks <- seq(-7,-25,-3)
  modern_lon_breaks <- seq(108,132,4)
  
  paleo_breaks <- sapply(modern_breaks, function(lat) {
    site_name <- paste0("lat_", abs(lat), "S")
    site_data <- paleo_df[paleo_df$Site == site_name, ]
    paleo_val <- approx(x = site_data$Age,y = site_data$Paleolatitude,
                        xout = current_age, rule = 2)$y
  })
  
  paleo_lon_breaks  <- sapply(modern_lon_breaks, function(lon) {
    site_name <- paste0("lon_", lon, "E")
    site_data <- paleo_lon_df[paleo_lon_df$Site == site_name, ]
    approx(x = site_data$Age, y = site_data$PaleoLon,
           xout = current_age, rule = 2)$y
  })
  
  paleo_grid <- seq(-10, -35, by = -3)
  paleo_lon_grid <- seq(108, 132, by = 4)
  
  tick_positions <- approx(x = paleo_breaks, y = modern_breaks,   
                           xout = paleo_grid, rule = 1)$y
  
  lon_tick_positions <- approx(x = paleo_lon_breaks,y = modern_lon_breaks,
                               xout = paleo_lon_grid, rule = 1)$y
  
  visible     <- !is.na(tick_positions) & tick_positions >= -23 & tick_positions <= -10
  tick_y      <- tick_positions[visible]
  tick_labels <- paste0(abs(paleo_grid[visible]), "°S")
  
  lon_visible     <- !is.na(lon_tick_positions) & lon_tick_positions >= 112 & lon_tick_positions <= 128
  lon_tick_x      <- lon_tick_positions[lon_visible]
  lon_tick_labels <- paste0(paleo_lon_grid[lon_visible], "°E")
  
  # Filter shapefiles whose range covers this time slice
  reefs_ts <- subset(All_reef_extents, from_ka >= ts_end & to_ka <= ts_start)
  
  # Define a custom color palette with multiple shades between red and blue
  custom_colors <- colorRampPalette(brewer.pal(10,"RdBu"))(21)
  
  # Create the ggplot object for the current time slice
  bm <- ggplot() +
    geom_raster(data = gebco_df, aes(x = x, y = y, fill = elevation)) +
    geom_contour(data = gebco_df, aes(x = x, y = y, z = elevation),
                 breaks = 0, color = "black", size = 0.3) +
    geom_sf(data = reefs_ts, aes(alpha = "Reef extent"), fill = "#FF6B6B",
            color = "black", show.legend = TRUE) +
    geom_point(data = reef_legend_df, aes(x = x, y = y, alpha = label),
               size = 0.001, inherit.aes = FALSE, show.legend = FALSE) +
    scale_alpha_manual(values = c("Reef extent" = 0.8), name = "") +
    geom_sf(data = ts_data_sf, aes(color = GR), size = 4.5) +
    scale_color_gradientn(colors = ipcc_precip,
                          limits = c(0, 50),
                          oob = scales::squish, name = "NGR") +
    scale_fill_gradientn(colors = topo_colors, values = scaled_breaks,
                         breaks = c(-7000,-4000,-1000,0,900),
                         labels = c("-7000","-4000","-1000","Sea-Level","900"),
                         limits = c(-7112, 994), name = "Elevation (m)") +   
    guides(
      color = guide_colorbar(order = 1, barwidth = 1, barheight = 13, title.position = "top", title.hjust = 0.5),
      fill = guide_colorbar(order = 2, barwidth = 1, barheight = 10, title.position = "top", title.hjust = 0.5),
      alpha = guide_legend(order = 3, override.aes = list(fill = "#FF6B6B", color = "black")),
    ) +
    coord_sf(xlim = c(112, 128), ylim = c(-23, -10), expand = FALSE, clip = "off") +
    annotate("text", x = 111.8, y = tick_y, label = tick_labels,
             hjust = 1, size  = 18 / .pt, color = "#4D4D4DFF", fontface = "plain") +
    annotate("text", x = lon_tick_x, y = -23.4, label = lon_tick_labels,
             vjust    = 1, size = 18 / .pt, color    = "#4D4D4DFF", fontface = "plain") +
    labs(x = "Paleolongitude", y = "Paleolatitude") +
    theme_minimal(base_size = 18) +
    theme(
      legend.box = "vertical",
      legend.box.just = "left",        
      legend.position = "right",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      plot.margin = margin(0, 0, 0, 60),
      axis.title = element_text(size = 20),
      axis.title.x = element_text(size = 20, vjust = -7.4),
      axis.title.y = element_text(size = 20, vjust = 14),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text = element_text(size = 18),
      axis.ticks = element_line(color = "black"),
      panel.background = element_rect(fill = "aliceblue"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )  +
    annotate("segment",x = 112, xend = 111.9, y = tick_y, yend = tick_y,
             linewidth = 0.75) +
    annotate("segment", x = lon_tick_x, xend = lon_tick_x, y = -23, yend = -23.15,
             linewidth = 0.75) #+
  
  ngr_panel <- make_ngr_panel(current_age)
  
  combined_plot <- wrap_elements(
    (ngr_panel | bm) + plot_layout(widths = c(1.1, 3)), clip = FALSE
  ) + theme(plot.margin = margin(90, 0, 70, 0))
  
  # Return the ggplot object
  return(combined_plot)
})

# Set the directory where you want to save the PNG files
output_dir <- "RScripts&Data/Time Slice Maps_NW_Australia/"

# Print or plot the time slice maps and save them as PNG files
for (i in seq_along(time_slice_plots)) {
  # Define the file name for the PNG file
  file_name <- paste0("time_slice_plot_", i, ".png")
  
  # Save the plot as a PNG file
  ggsave(file.path(output_dir, file_name), plot = time_slice_plots[[i]], width = 10, height = 8, units = "in", dpi = 600, bg = "white")
}

png_files <- list.files(output_dir, pattern = 'time_slice_plot_.*\\.png$', full.names = TRUE)
sorted_png_files <- png_files[order(as.numeric(gsub(".*time_slice_plot_(\\d+)\\.png$", "\\1", png_files)))]

av::av_encode_video(sorted_png_files, framerate = 4, output = 'RScripts&Data/Time Slice Maps_NW_Australia/Ts_Gif.mp4', vfilter = "scale=1920:-1")

# Import PNG files and create a GIF

dir <- "RScripts&Data/Time Slice Maps_NW_Australia/Time_Slices/"
png_files <- list.files(dir, pattern = 'time_slice_plot_.*\\.png$', full.names = TRUE)
sorted_png_files <- png_files[order(as.numeric(gsub(".*time_slice_plot_(\\d+)\\.png$", "\\1", png_files)))]


annotated_indices <- c(1, 2, 3, 4, 5, 6, 11, 30, 53, 88, 114, 164, 191)  
freeze_indices <- c(6, 11, 30, 53, 88, 114, 164, 191)  
freeze_duration <- 28  

r_frame <- png::readPNG(sorted_png_files[[190]])
target_w <- dim(r_frame)[2]
target_h <- dim(r_frame)[1]

for (i in annotated_indices) {
  img <- magick::image_read(sorted_png_files[[i]])
  img <- magick::image_resize(img, paste0(target_w, "x", target_h, "!"))
  magick::image_write(img, path = sorted_png_files[[i]])
}

expanded_files <- c()

for (f in sorted_png_files) {
  expanded_files <- c(expanded_files, f)
  idx <- as.numeric(gsub(".*time_slice_plot_(\\d+)\\.png$", "\\1", f))
  if (idx %in% freeze_indices) {
    expanded_files <- c(expanded_files, rep(f, freeze_duration))
  }
}
av::av_encode_video(expanded_files, framerate = 4, output = 'RScripts&Data/Time Slice Maps_NW_Australia/Time_Slices/Ts_Gif.mp4', vfilter = "scale=1920:-1")
