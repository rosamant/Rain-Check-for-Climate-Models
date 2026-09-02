install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Gorgonichthys 1 and Kaleidoscope 1 datasets

Gorgonichthys1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Gorgonichthys1.csv", header=TRUE, stringsAsFactors=FALSE)
Gorgonichthys1=Gorgonichthys1[c(168:6001),] # Oligocene-Miocene
head(Gorgonichthys1)
plot(Gorgonichthys1, type="l", xlim = c(250, 1200), ylim = c(0, 30))

Kaleidoscope1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Kaleidoscope1.csv", header=TRUE, stringsAsFactors=FALSE)
Kaleidoscope1=Kaleidoscope1[c(1:4112),] # Oligocene-Miocene
head(Kaleidoscope1)
plot(Kaleidoscope1, type="l", xlim = c(300, 1050), ylim = c(0, 30))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Gorgonichthys1_interpolated <- linterp(Gorgonichthys1, dt = 0.2, genplot = F)
Kaleidoscope1_interpolated <- linterp(Kaleidoscope1, dt = 0.2, genplot = F)

# Scaling the data
Gmean = Gmean(Gorgonichthys1_interpolated$GR)
Gstd = Gsd(Gorgonichthys1_interpolated$GR)
Gorgonichthys1_scaled = (Gorgonichthys1_interpolated$GR - Gmean)/Gstd
Gorgonichthys1_rescaled = data.frame(Gorgonichthys1_interpolated$DEPT, Gorgonichthys1_scaled)

Kmean = Gmean(Kaleidoscope1_interpolated$GR)
Kstd = Gsd(Kaleidoscope1_interpolated$GR)
Kaleidoscope1_scaled = (Kaleidoscope1_interpolated$GR - Kmean)/Kstd
Kaleidoscope1_rescaled = data.frame(Kaleidoscope1_interpolated$DEPT, Kaleidoscope1_scaled)

# Resampling the data using moving window statistics
Gorgonichthys1_scaled = mwStats(Gorgonichthys1_rescaled, cols = 2, win=3, ends = T)
Gorgonichthys1_standardized = data.frame(Gorgonichthys1_scaled$Center_win, Gorgonichthys1_scaled$Average)

Kaleidoscope1_scaled = mwStats(Kaleidoscope1_rescaled, cols = 2, win=3, ends = T)
Kaleidoscope1_standardized = data.frame(Kaleidoscope1_scaled$Center_win, Kaleidoscope1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Gorgonichthys1_standardized, type="l", xlim = c(250, 1200), ylim = c(-20, 20), xlab = "Gorgonichthys1 Resampled Depth", ylab = "Normalized GR")
plot(Kaleidoscope1_standardized, type="l", xlim = c(300, 1050), ylim = c(-10, 10), xlab = "Kaleidoscope 1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_k1_g1_ap1 <- dtw(Kaleidoscope1_standardized$Kaleidoscope1_scaled.Average, Gorgonichthys1_standardized$Gorgonichthys1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = F))
plot(al_k1_g1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Kaleidoscope1_on_Gorgonichthys1_depth = tune(Kaleidoscope1_standardized, cbind(Kaleidoscope1_standardized$Kaleidoscope1_scaled.Center_win[al_k1_g1_ap1$index1s], Gorgonichthys1_standardized$Gorgonichthys1_scaled.Center_win[al_k1_g1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Gorgonichthys1_standardized, type = "l", ylim = c(-20, 20), xlim = c(250, 1200), xlab = "Gorgonichthys1 Resampled Depth", ylab = "Normalized GR")
lines(Kaleidoscope1_on_Gorgonichthys1_depth, col = "red")

# DTW Distance
al_k1_g1_ap1$normalizedDistance
al_k1_g1_ap1$distance

# Tuning Kaleidoscope1 data on Walkley1 depth
Kaleidoscope1_on_Walkley1_depth = tune(Kaleidoscope1_on_Gorgonichthys1_depth, cbind(Gorgonichthys1_standardized$Gorgonichthys1_scaled.Center_win[al_g1_w1_ap2$index1s], Walkley1_standardized$Walkley1_scaled.Center_win[al_g1_w1_ap2$index2s]), extrapolate = F)

# Tuning Kaleidoscope1 data on Caswell1 depth
Kaleidoscope1_on_Caswell1_depth = tune(Kaleidoscope1_on_Walkley1_depth, cbind(Walkley1_standardized$Walkley1_scaled.Center_win[al_w1_c1_ap1$index1s], Caswell1_standardized$Caswell1_scaled.Center_win[al_w1_c1_ap1$index2s]), extrapolate = F)

# Tuning Kaleidoscope1 data on Calliance2 depth
Kaleidoscope1_on_Calliance2_depth = tune(Kaleidoscope1_on_Caswell1_depth, cbind(Caswell1_standardized$Caswell1_scaled.Center_win[al_c2_c1_ap1$index2s], Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_c1_ap1$index1s]), extrapolate = F)

# Tuning Kaleidoscope1 data on Omar1 depth
Kaleidoscope1_on_Omar1_depth = tune(Kaleidoscope1_on_Calliance2_depth, cbind(Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_o1_ap2$index1s], Omar1_standardized$Omar1_scaled.Center_win[al_c2_o1_ap2$index2s]), extrapolate = F)

# Tuning Kaleidoscope1 data on SG1 depth
Kaleidoscope1_on_SG1_depth = tune(Kaleidoscope1_on_Omar1_depth, cbind(Omar1_standardized$Omar1_scaled.Center_win[al_o1_sg1_ap2$index1s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_o1_sg1_ap2$index2s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 20), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (Kaleidoscope-1)")
lines(Kaleidoscope1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Kaleidoscope1_originalGR_on_SG1_depth = data.frame(Kaleidoscope1_on_SG1_depth$X1, Kaleidoscope1_interpolated[1088:3817,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (Kaleidoscope-1)")
lines(Kaleidoscope1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Kaleidoscope1 

SG1Age_on_Kaleidoscope1_depth = tune(Kaleidoscope1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Kaleidoscope1_depth, type = "l", ylim = c(0, 60), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Kaleidoscope1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Kaleidoscope1_depth) <- new_column_names
SG1Age_on_Kaleidoscope1_depth[,1] = SG1Age_on_Kaleidoscope1_depth[,1] * 1000
write.csv(SG1Age_on_Kaleidoscope1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Kaleidoscope 1.csv", row.names = FALSE)
