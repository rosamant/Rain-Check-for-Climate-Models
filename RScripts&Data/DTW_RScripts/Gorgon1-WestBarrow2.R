install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Gorgon1 and WestBarrow2 datasets

Gorgon1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Gorgon1.csv", header=TRUE, stringsAsFactors=FALSE)
Gorgon1=Gorgon1[c(1:9280),] # Oligocene-Miocene
head(Gorgon1)
plot(Gorgon1, type="l", xlim = c(300, 1700), ylim = c(0, 50))

WestBarrow2 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/WestBarrow-2.csv", header=TRUE, stringsAsFactors=FALSE)
WestBarrow2=WestBarrow2[c(1:5539),] # Oligocene-Miocene
head(WestBarrow2)
plot(WestBarrow2, type="l", xlim = c(100, 1000), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Gorgon1_interpolated <- linterp(Gorgon1, dt = 0.2, genplot = F)
WestBarrow2_interpolated <- linterp(WestBarrow2, dt = 0.2, genplot = F)

# Scaling the data
Gmean = Gmean(Gorgon1_interpolated$GR)
Gstd = Gsd(Gorgon1_interpolated$GR)
Gorgon1_scaled = (Gorgon1_interpolated$GR - Gmean)/Gstd
Gorgon1_rescaled = data.frame(Gorgon1_interpolated$DEPT, Gorgon1_scaled)

Wmean = Gmean(WestBarrow2_interpolated$GR)
Wstd = Gsd(WestBarrow2_interpolated$GR)
WestBarrow2_scaled = (WestBarrow2_interpolated$GR - Wmean)/Wstd
WestBarrow2_rescaled = data.frame(WestBarrow2_interpolated$DEPT, WestBarrow2_scaled)

# Resampling the data using moving window statistics
Gorgon1_scaled = mwStats(Gorgon1_rescaled, cols = 2, win=3, ends = T)
Gorgon1_standardized = data.frame(Gorgon1_scaled$Center_win, Gorgon1_scaled$Average)

WestBarrow2_scaled = mwStats(WestBarrow2_rescaled, cols = 2, win=3, ends = T)
WestBarrow2_standardized = data.frame(WestBarrow2_scaled$Center_win, WestBarrow2_scaled$Average)

# Plotting the rescaled and resampled data
plot(Gorgon1_standardized, type="l", xlim = c(300, 1700), ylim = c(-20, 20), xlab = "Gorgon1 Resampled Depth", ylab = "Normalized GR")
plot(WestBarrow2_standardized, type="l", xlim = c(100, 1000), ylim = c(-20, 20), xlab = "WestBarrow2 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_wb1_g1_ap1 <- dtw(WestBarrow2_standardized$WestBarrow2_scaled.Average, Gorgon1_standardized$Gorgon1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_wb1_g1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
WestBarrow2_on_Gorgon1_depth = tune(WestBarrow2_standardized, cbind(WestBarrow2_standardized$WestBarrow2_scaled.Center_win[al_wb1_g1_ap1$index1s], Gorgon1_standardized$Gorgon1_scaled.Center_win[al_wb1_g1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Gorgon1_standardized, type = "l", ylim = c(-20, 20), xlim = c(300, 1700), xlab = "Gorgon1 Resampled Depth", ylab = "Normalized GR")
lines(WestBarrow2_on_Gorgon1_depth, col = "red")

# DTW Distance
al_wb1_g1_ap1$normalizedDistance
al_wb1_g1_ap1$distance

#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(WestBarrow2_standardized),ncol=nrow(Gorgon1_standardized))
image(x=Gorgon1_standardized[,1],y=WestBarrow2_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Gorgon1_standardized[,1] - 380))
base_1_y <- which.min(abs(WestBarrow2_standardized[,1] - 150))

# Depth values for second datum
base_2_x <- which.min(abs(Gorgon1_standardized[,1] - 525))
base_2_y <- which.min(abs(WestBarrow2_standardized[,1] - 300))

# Depth values for third datum
base_3_x <- which.min(abs(Gorgon1_standardized[,1] - 625))
base_3_y <- which.min(abs(WestBarrow2_standardized[,1] - 380))

# Depth values for fourth datum
base_4_x <- which.min(abs(Gorgon1_standardized[,1] - 870))
base_4_y <- which.min(abs(WestBarrow2_standardized[,1] - 550))

# Depth values for fifth datum
base_5_x <- which.min(abs(Gorgon1_standardized[,1] - 1550))
base_5_y <- which.min(abs(WestBarrow2_standardized[,1] - 950))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(WestBarrow2_standardized), ncol = nrow(Gorgon1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+300):nrow(WestBarrow2_standardized),1:(base_1_x-300)] <- 0
compare.window[1:(base_1_y-100),(base_1_x+100):ncol(compare.window)] <- 0

compare.window[(base_2_y+300):nrow(WestBarrow2_standardized),1:(base_2_x-300)] <- 0
compare.window[1:(base_2_y-200),(base_2_x+200):ncol(compare.window)] <- 0

compare.window[(base_3_y+300):nrow(WestBarrow2_standardized),1:(base_3_x-300)] <- 0
compare.window[1:(base_3_y-300),(base_3_x+300):ncol(compare.window)] <- 0

compare.window[(base_4_y+400):nrow(WestBarrow2_standardized),1:(base_4_x-400)] <- 0
compare.window[1:(base_4_y-300),(base_4_x+300):ncol(compare.window)] <- 0

compare.window[(base_5_y+200):nrow(WestBarrow2_standardized),1:(base_5_x-250)] <- 0
compare.window[1:(base_5_y-300),(base_5_x+300):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Gorgon1_standardized[,1],y=WestBarrow2_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Gorgon1_standardized[,1],y=WestBarrow2_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_wb1_g1_ap1 <- dtw(WestBarrow2_standardized$WestBarrow2_scaled.Average, Gorgon1_standardized$Gorgon1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = F, open.begin = F))
plot(al_wb1_g1_ap1, type = "threeway")

# DTW Distance measure
al_wb1_g1_ap1$normalizedDistance
al_wb1_g1_ap1$distance

image(y = Gorgon1_standardized[,1], x = WestBarrow2_standardized[,1], z = compare.window, useRaster = T)
lines(WestBarrow2_standardized$WestBarrow2_scaled.Center_win[al_wb1_g1_ap1$index1], Gorgon1_standardized$Gorgon1_scaled.Center_win[al_wb1_g1_ap1$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
WestBarrow2_on_Gorgon1_depth = tune(WestBarrow2_standardized, cbind(WestBarrow2_standardized$WestBarrow2_scaled.Center_win[al_wb1_g1_ap1$index1s], Gorgon1_standardized$Gorgon1_scaled.Center_win[al_wb1_g1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Gorgon1_standardized, type = "l", ylim = c(-20, 20), xlim = c(300, 1700), xlab = "Gorgon1 Resampled Depth", ylab = "Normalized GR (West Barrow-2)")
lines(WestBarrow2_on_Gorgon1_depth, col = "red")

# Tuning WestBarrow2 data on Bluebell1 depth
WestBarrow2_on_Bluebell1_depth = tune(WestBarrow2_on_Gorgon1_depth, cbind(Gorgon1_standardized$Gorgon1_scaled.Center_win[al_g1_b1_ap1$index1s], Bluebell1_standardized$Bluebell1_scaled.Center_win[al_g1_b1_ap1$index2s]), extrapolate = F)

# Tuning WestBarrow2 data on Fisher1 depth
WestBarrow2_on_Fisher1_depth = tune(WestBarrow2_on_Bluebell1_depth, cbind(Bluebell1_standardized$Bluebell1_scaled.Center_win[al_bl1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_bl1_fi1_ap1$index2s]), extrapolate = F)

# Tuning WestBarrow2 data on Finucane1 depth
WestBarrow2_on_Finucane1_depth = tune(WestBarrow2_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning WestBarrow2 data on Picard1 depth
WestBarrow2_on_Picard1_depth = tune(WestBarrow2_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (West Barrow-2)")
lines(WestBarrow2_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
WestBarrow2_originalGR_on_Picard1_depth = data.frame(WestBarrow2_on_Picard1_depth$X1, WestBarrow2_interpolated[1:4455,2])

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (West Barrow-2)")
lines(WestBarrow2_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to WestBarrow2 

U1463Age_on_WestBarrow2_depth = tune(WestBarrow2_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_WestBarrow2_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "WestBarrow2")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_WestBarrow2_depth) <- new_column_names
write.csv(U1463Age_on_WestBarrow2_depth, file = "RScripts&Data/Sites Data_Age-NGR/West Barrow 2.csv", row.names = FALSE)


WestBarrow2_age_depth <- approx(x = WestBarrow2_interpolated$GR,
                            y = WestBarrow2_interpolated$DEPT,
                            xout = U1463Age_on_WestBarrow2_depth$GR,
                            rule = 1)$y

WestBarrow2_agemodel = data.frame(WestBarrow2_age_depth, U1463Age_on_WestBarrow2_depth$AGE)
WestBarrow2_agemodel <- na.omit(WestBarrow2_agemodel)
new_column_names1 <- c("Depth", "Age")
colnames(WestBarrow2_agemodel) <- new_column_names1
plot(WestBarrow2_agemodel)
write.csv(WestBarrow2_agemodel, file = "RScripts&Data/Sites Data_Age-NGR/WestBarrow2_DepthAge.csv", row.names = FALSE)
