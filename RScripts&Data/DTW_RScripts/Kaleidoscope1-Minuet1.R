install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Kaleidoscope 1 and Minuet 1 datasets

Kaleidoscope1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Kaleidoscope1.csv", header=TRUE, stringsAsFactors=FALSE)
Kaleidoscope1=Kaleidoscope1[c(1:4112),] # Oligocene-Miocene
head(Kaleidoscope1)
plot(Kaleidoscope1, type="l", xlim = c(300, 1050), ylim = c(0, 30))

Minuet1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Minuet 1.csv", header=TRUE, stringsAsFactors=FALSE)
Minuet1=Minuet1[c(531:6441),] # Oligocene-Miocene
head(Minuet1)
plot(Minuet1, type="l", xlim = c(300, 900), ylim = c(0, 30))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Kaleidoscope1_interpolated <- linterp(Kaleidoscope1, dt = 0.2, genplot = F)
Minuet1_interpolated <- linterp(Minuet1, dt = 0.2, genplot = F)

# Scaling the data
Kmean = Gmean(Kaleidoscope1_interpolated$GR)
Kstd = Gsd(Kaleidoscope1_interpolated$GR)
Kaleidoscope1_scaled = (Kaleidoscope1_interpolated$GR - Kmean)/Kstd
Kaleidoscope1_rescaled = data.frame(Kaleidoscope1_interpolated$DEPT, Kaleidoscope1_scaled)

Mimean = Gmean(Minuet1_interpolated$GR)
Mistd = Gsd(Minuet1_interpolated$GR)
Minuet1_scaled = (Minuet1_interpolated$GR - Mimean)/Mistd
Minuet1_rescaled = data.frame(Minuet1_interpolated$DEPT, Minuet1_scaled)

# Resampling the data using moving window statistics
Kaleidoscope1_scaled = mwStats(Kaleidoscope1_rescaled, cols = 2, win=3, ends = T)
Kaleidoscope1_standardized = data.frame(Kaleidoscope1_scaled$Center_win, Kaleidoscope1_scaled$Average)

Minuet1_scaled = mwStats(Minuet1_rescaled, cols = 2, win=3, ends = T)
Minuet1_standardized = data.frame(Minuet1_scaled$Center_win, Minuet1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Kaleidoscope1_standardized, type="l", xlim = c(300, 1050), ylim = c(-20, 20), xlab = "Kaleidoscope1 Resampled Depth", ylab = "Normalized GR")
plot(Minuet1_standardized, type="l", xlim = c(300, 900), ylim = c(-10, 10), xlab = "Minuet 1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_mi1_k1_ap1 <- dtw(Minuet1_standardized$Minuet1_scaled.Average, Kaleidoscope1_standardized$Kaleidoscope1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = F, open.end = F))
plot(al_mi1_k1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Minuet1_on_Kaleidoscope1_depth = tune(Minuet1_standardized, cbind(Minuet1_standardized$Minuet1_scaled.Center_win[al_mi1_k1_ap1$index1s], Kaleidoscope1_standardized$Kaleidoscope1_scaled.Center_win[al_mi1_k1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Kaleidoscope1_standardized, type = "l", ylim = c(-10, 10), xlim = c(300, 1050), xlab = "Kaleidoscope1 Resampled Depth", ylab = "Normalized GR")
lines(Minuet1_on_Kaleidoscope1_depth, col = "red")

# DTW Distance
al_mi1_k1_ap1$normalizedDistance
al_mi1_k1_ap1$distance

# Tuning Minuet1 data on Gorgonichthys1 depth
Minuet1_on_Gorgonichthys1_depth = tune(Minuet1_on_Kaleidoscope1_depth, cbind(Kaleidoscope1_standardized$Kaleidoscope1_scaled.Center_win[al_k1_g1_ap1$index1s], Gorgonichthys1_standardized$Gorgonichthys1_scaled.Center_win[al_k1_g1_ap1$index2s]), extrapolate = F)

# Tuning Minuet1 data on Walkley1 depth
Minuet1_on_Walkley1_depth = tune(Minuet1_on_Gorgonichthys1_depth, cbind(Gorgonichthys1_standardized$Gorgonichthys1_scaled.Center_win[al_g1_w1_ap2$index1s], Walkley1_standardized$Walkley1_scaled.Center_win[al_g1_w1_ap2$index2s]), extrapolate = F)

# Tuning Minuet1 data on Caswell1 depth
Minuet1_on_Caswell1_depth = tune(Minuet1_on_Walkley1_depth, cbind(Walkley1_standardized$Walkley1_scaled.Center_win[al_w1_c1_ap1$index1s], Caswell1_standardized$Caswell1_scaled.Center_win[al_w1_c1_ap1$index2s]), extrapolate = F)

# Tuning Minuet1 data on Calliance2 depth
Minuet1_on_Calliance2_depth = tune(Minuet1_on_Caswell1_depth, cbind(Caswell1_standardized$Caswell1_scaled.Center_win[al_c2_c1_ap1$index2s], Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_c1_ap1$index1s]), extrapolate = F)

# Tuning Minuet1 data on Omar1 depth
Minuet1_on_Omar1_depth = tune(Minuet1_on_Calliance2_depth, cbind(Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_o1_ap2$index1s], Omar1_standardized$Omar1_scaled.Center_win[al_c2_o1_ap2$index2s]), extrapolate = F)

# Tuning Minuet1 data on SG1 depth
Minuet1_on_SG1_depth = tune(Minuet1_on_Omar1_depth, cbind(Omar1_standardized$Omar1_scaled.Center_win[al_o1_sg1_ap2$index1s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_o1_sg1_ap2$index2s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 20), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (Minuet-1)")
lines(Minuet1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Minuet1_originalGR_on_SG1_depth = data.frame(Minuet1_on_SG1_depth$X1, Minuet1_interpolated[893:2885,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (Minuet-1)")
lines(Minuet1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Minuet1 

SG1Age_on_Minuet1_depth = tune(Minuet1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Minuet1_depth, type = "l", ylim = c(0, 60), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Minuet1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Minuet1_depth) <- new_column_names
SG1Age_on_Minuet1_depth[,1] = SG1Age_on_Minuet1_depth[,1] * 1000
write.csv(SG1Age_on_Minuet1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Minuet 1.csv", row.names = FALSE)
