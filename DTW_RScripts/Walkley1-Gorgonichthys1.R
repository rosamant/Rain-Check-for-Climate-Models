install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Walkley1 and Gorgonichthys1 datasets

Walkley1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Walkley 1.csv", header=TRUE, stringsAsFactors=FALSE)
Walkley1=Walkley1[c(6:6312),] # Oligocene-Miocene
head(Walkley1)
plot(Walkley1, type="l", xlim = c(350, 1600), ylim = c(0, 60))

Gorgonichthys1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Gorgonichthys1.csv", header=TRUE, stringsAsFactors=FALSE)
Gorgonichthys1=Gorgonichthys1[c(168:6001),] # Oligocene-Miocene
head(Gorgonichthys1)
plot(Gorgonichthys1, type="l", xlim = c(250, 1200), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Walkley1_interpolated <- linterp(Walkley1, dt = 0.2, genplot = F)
Gorgonichthys1_interpolated <- linterp(Gorgonichthys1, dt = 0.2, genplot = F)

# Scaling the data
Wmean = Gmean(Walkley1_interpolated$GR)
Wstd = Gsd(Walkley1_interpolated$GR)
Walkley1_scaled = (Walkley1_interpolated$GR - Wmean)/Wstd
Walkley1_rescaled = data.frame(Walkley1_interpolated$DEPT, Walkley1_scaled)

Gmean = Gmean(Gorgonichthys1_interpolated$GR)
Gstd = Gsd(Gorgonichthys1_interpolated$GR)
Gorgonichthys1_scaled = (Gorgonichthys1_interpolated$GR - Gmean)/Gstd
Gorgonichthys1_rescaled = data.frame(Gorgonichthys1_interpolated$DEPT, Gorgonichthys1_scaled)

# Resampling the data using moving window statistics
Walkley1_scaled = mwStats(Walkley1_rescaled, cols = 2, win=3, ends = T)
Walkley1_standardized = data.frame(Walkley1_scaled$Center_win, Walkley1_scaled$Average)

Gorgonichthys1_scaled = mwStats(Gorgonichthys1_rescaled, cols = 2, win=3, ends = T)
Gorgonichthys1_standardized = data.frame(Gorgonichthys1_scaled$Center_win, Gorgonichthys1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Walkley1_standardized, type="l", xlim = c(350, 1600), ylim = c(-20, 20), xlab = "Walkley1 Resampled Depth", ylab = "Normalized GR")
plot(Gorgonichthys1_standardized, type="l", xlim = c(250, 1200), ylim = c(-20, 20), xlab = "Gorgonichthys1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_g1_w1_ap1 <- dtw(Gorgonichthys1_standardized$Gorgonichthys1_scaled.Average, Walkley1_standardized$Walkley1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_g1_w1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Gorgonichthys1_on_Walkley1_depth = tune(Gorgonichthys1_standardized, cbind(Gorgonichthys1_standardized$Gorgonichthys1_scaled.Center_win[al_g1_w1_ap1$index1s], Walkley1_standardized$Walkley1_scaled.Center_win[al_g1_w1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Walkley1_standardized, type = "l", ylim = c(-20, 20), xlim = c(350, 1600), xlab = "Gorgonichthys1 Resampled Depth", ylab = "Normalized GR")
lines(Gorgonichthys1_on_Walkley1_depth, col = "red")

# DTW Distance
al_g1_w1_ap1$normalizedDistance
al_g1_w1_ap1$distance


#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(Gorgonichthys1_standardized),ncol=nrow(Walkley1_standardized))
image(x=Walkley1_standardized[,1],y=Gorgonichthys1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Walkley1_standardized[,1] - 670))
base_1_y <- which.min(abs(Gorgonichthys1_standardized[,1] - 583))

# Depth values for second datum
base_2_x <- which.min(abs(Walkley1_standardized[,1] - 843))
base_2_y <- which.min(abs(Gorgonichthys1_standardized[,1] - 737))

# Depth values for third datum
base_3_x <- which.min(abs(Walkley1_standardized[,1] - 1054))
base_3_y <- which.min(abs(Gorgonichthys1_standardized[,1] - 868))

# Depth values for fourth datum
base_4_x <- which.min(abs(Walkley1_standardized[,1] - 1192))
base_4_y <- which.min(abs(Gorgonichthys1_standardized[,1] - 899))

# Depth values for fifth datum
base_5_x <- which.min(abs(Walkley1_standardized[,1] - 1306))
base_5_y <- which.min(abs(Gorgonichthys1_standardized[,1] - 981))

# Depth values for sixth datum
base_6_x <- which.min(abs(Walkley1_standardized[,1] - 1458))
base_6_y <- which.min(abs(Gorgonichthys1_standardized[,1] - 1108))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(Gorgonichthys1_standardized), ncol = nrow(Walkley1_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+200):nrow(Gorgonichthys1_standardized),1:(base_1_x-200)] <- 0
compare.window[1:(base_1_y-200),(base_1_x+200):ncol(compare.window)] <- 0

compare.window[(base_2_y+200):nrow(Gorgonichthys1_standardized),1:(base_2_x-200)] <- 0
compare.window[1:(base_2_y-200),(base_2_x+200):ncol(compare.window)] <- 0

compare.window[(base_3_y+80):nrow(Gorgonichthys1_standardized),1:(base_3_x-80)] <- 0
compare.window[1:(base_3_y-80),(base_3_x+50):ncol(compare.window)] <- 0

compare.window[(base_4_y+80):nrow(Gorgonichthys1_standardized),1:(base_4_x-80)] <- 0
compare.window[1:(base_4_y-80),(base_4_x+80):ncol(compare.window)] <- 0

compare.window[(base_5_y+80):nrow(Gorgonichthys1_standardized),1:(base_5_x-80)] <- 0
compare.window[1:(base_5_y-80),(base_5_x+80):ncol(compare.window)] <- 0

compare.window[(base_6_y+80):nrow(Gorgonichthys1_standardized),1:(base_6_x-80)] <- 0
compare.window[1:(base_6_y-80),(base_6_x+80):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Walkley1_standardized[,1],y=Gorgonichthys1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Walkley1_standardized[,1],y=Gorgonichthys1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_g1_w1_ap2 <- dtw(Gorgonichthys1_standardized$Gorgonichthys1_scaled.Average, Walkley1_standardized$Walkley1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = F, open.begin = F))
plot(al_g1_w1_ap2, type = "threeway")

# DTW Distance measure
al_g1_w1_ap2$normalizedDistance
al_g1_w1_ap2$distance

image(y = Walkley1_standardized[,1], x = Gorgonichthys1_standardized[,1], z = compare.window, useRaster = T)
lines(Gorgonichthys1_standardized$Gorgonichthys1_scaled.Center_win[al_g1_w1_ap2$index1], Walkley1_standardized$Walkley1_scaled.Center_win[al_g1_w1_ap2$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
Gorgonichthys1_on_Walkley1_depth1 = tune(Gorgonichthys1_standardized, cbind(Gorgonichthys1_standardized$Gorgonichthys1_scaled.Center_win[al_g1_w1_ap2$index1s], Walkley1_standardized$Walkley1_scaled.Center_win[al_g1_w1_ap2$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Walkley1_standardized, type = "l", ylim = c(-20, 20), xlim = c(350, 1600), xlab = "Walkley-1 Resampled Depth", ylab = "Normalized GR (Gorgonichthys-1)")
lines(Gorgonichthys1_on_Walkley1_depth1, col = "red")

# Tuning Gorgonichthys1 data on Caswell1 depth
Gorgonichthys1_on_Caswell1_depth = tune(Gorgonichthys1_on_Walkley1_depth1, cbind(Walkley1_standardized$Walkley1_scaled.Center_win[al_w1_c1_ap1$index1s], Caswell1_standardized$Caswell1_scaled.Center_win[al_w1_c1_ap1$index2s]), extrapolate = F)

# Tuning Gorgonichthys1 data on Calliance2 depth
Gorgonichthys1_on_Calliance2_depth = tune(Gorgonichthys1_on_Caswell1_depth, cbind(Caswell1_standardized$Caswell1_scaled.Center_win[al_c2_c1_ap1$index2s], Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_c1_ap1$index1s]), extrapolate = F)

# Tuning Gorgonichthys1 data on Omar1 depth
Gorgonichthys1_on_Omar1_depth = tune(Gorgonichthys1_on_Calliance2_depth, cbind(Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_o1_ap2$index1s], Omar1_standardized$Omar1_scaled.Center_win[al_c2_o1_ap2$index2s]), extrapolate = F)

# Tuning Gorgonichthys1 data on SG1 depth
Gorgonichthys1_on_SG1_depth = tune(Gorgonichthys1_on_Omar1_depth, cbind(Omar1_standardized$Omar1_scaled.Center_win[al_o1_sg1_ap2$index1s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_o1_sg1_ap2$index2s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 20), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (Gorgonichthys-1)")
lines(Gorgonichthys1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Gorgonichthys1_originalGR_on_SG1_depth = data.frame(Gorgonichthys1_on_SG1_depth$X1, Gorgonichthys1_interpolated[966:4385,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (Gorgonichthys-1)")
lines(Gorgonichthys1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Gorgonichthys1 

SG1Age_on_Gorgonichthys1_depth = tune(Gorgonichthys1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Gorgonichthys1_depth, type = "l", ylim = c(0, 60), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Gorgonichthys1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Gorgonichthys1_depth) <- new_column_names
SG1Age_on_Gorgonichthys1_depth[,1] = SG1Age_on_Gorgonichthys1_depth[,1] * 1000
write.csv(SG1Age_on_Gorgonichthys1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Gorgonichthys 1.csv", row.names = FALSE)
