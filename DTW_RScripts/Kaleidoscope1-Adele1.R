install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Kaleidoscope 1 and Adele 1 datasets

Kaleidoscope1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Kaleidoscope1.csv", header=TRUE, stringsAsFactors=FALSE)
Kaleidoscope1=Kaleidoscope1[c(1:4112),] # Oligocene-Miocene
head(Kaleidoscope1)
plot(Kaleidoscope1, type="l", xlim = c(300, 1050), ylim = c(0, 30))

Adele1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Adele 1.csv", header=TRUE, stringsAsFactors=FALSE)
Adele1=Adele1[c(100:4370),] # Oligocene-Miocene
head(Adele1)
Adele1$DEPT <- as.numeric(Adele1$DEPT)
plot(Adele1, type="l", xlim = c(250, 950), ylim = c(0, 50))


#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Kaleidoscope1_interpolated <- linterp(Kaleidoscope1, dt = 0.2, genplot = F)
Adele1_interpolated <- linterp(Adele1, dt = 0.2, genplot = F)

# Scaling the data
Kmean = Gmean(Kaleidoscope1_interpolated$GR)
Kstd = Gsd(Kaleidoscope1_interpolated$GR)
Kaleidoscope1_scaled = (Kaleidoscope1_interpolated$GR - Kmean)/Kstd
Kaleidoscope1_rescaled = data.frame(Kaleidoscope1_interpolated$DEPT, Kaleidoscope1_scaled)

Amean = Gmean(Adele1_interpolated$GR)
Astd = Gsd(Adele1_interpolated$GR)
Adele1_scaled = (Adele1_interpolated$GR - Amean)/Astd
Adele1_rescaled = data.frame(Adele1_interpolated$DEPT, Adele1_scaled)

# Resampling the data using moving window statistics
Kaleidoscope1_scaled = mwStats(Kaleidoscope1_rescaled, cols = 2, win=3, ends = T)
Kaleidoscope1_standardized = data.frame(Kaleidoscope1_scaled$Center_win, Kaleidoscope1_scaled$Average)

Adele1_scaled = mwStats(Adele1_rescaled, cols = 2, win=3, ends = T)
Adele1_standardized = data.frame(Adele1_scaled$Center_win, Adele1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Kaleidoscope1_standardized, type="l", xlim = c(300, 1050), ylim = c(-20, 20), xlab = "Kaleidoscope1 Resampled Depth", ylab = "Normalized GR")
plot(Adele1_standardized, type="l", xlim = c(250, 900), ylim = c(-10, 10), xlab = "Adele 1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_a1_k1_ap1 <- dtw(Adele1_standardized$Adele1_scaled.Average, Kaleidoscope1_standardized$Kaleidoscope1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = F))
plot(al_a1_k1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Adele1_on_Kaleidoscope1_depth = tune(Adele1_standardized, cbind(Adele1_standardized$Adele1_scaled.Center_win[al_a1_k1_ap1$index1s], Kaleidoscope1_standardized$Kaleidoscope1_scaled.Center_win[al_a1_k1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Kaleidoscope1_standardized, type = "l", ylim = c(-10, 10), xlim = c(300, 1050), xlab = "Kaleidoscope1 Resampled Depth", ylab = "Normalized GR")
lines(Adele1_on_Kaleidoscope1_depth, col = "red")

# DTW Distance
al_a1_k1_ap1$normalizedDistance
al_a1_k1_ap1$distance

# Tuning Adele1 data on Gorgonichthys1 depth
Adele1_on_Gorgonichthys1_depth = tune(Adele1_on_Kaleidoscope1_depth, cbind(Kaleidoscope1_standardized$Kaleidoscope1_scaled.Center_win[al_k1_g1_ap1$index1s], Gorgonichthys1_standardized$Gorgonichthys1_scaled.Center_win[al_k1_g1_ap1$index2s]), extrapolate = F)

# Tuning Adele1 data on Walkley1 depth
Adele1_on_Walkley1_depth = tune(Adele1_on_Gorgonichthys1_depth, cbind(Gorgonichthys1_standardized$Gorgonichthys1_scaled.Center_win[al_g1_w1_ap2$index1s], Walkley1_standardized$Walkley1_scaled.Center_win[al_g1_w1_ap2$index2s]), extrapolate = F)

# Tuning Adele1 data on Caswell1 depth
Adele1_on_Caswell1_depth = tune(Adele1_on_Walkley1_depth, cbind(Walkley1_standardized$Walkley1_scaled.Center_win[al_w1_c1_ap1$index1s], Caswell1_standardized$Caswell1_scaled.Center_win[al_w1_c1_ap1$index2s]), extrapolate = F)

# Tuning Adele1 data on Calliance2 depth
Adele1_on_Calliance2_depth = tune(Adele1_on_Caswell1_depth, cbind(Caswell1_standardized$Caswell1_scaled.Center_win[al_c2_c1_ap1$index2s], Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_c1_ap1$index1s]), extrapolate = F)

# Tuning Adele1 data on Omar1 depth
Adele1_on_Omar1_depth = tune(Adele1_on_Calliance2_depth, cbind(Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_o1_ap2$index1s], Omar1_standardized$Omar1_scaled.Center_win[al_c2_o1_ap2$index2s]), extrapolate = F)

# Tuning Adele1 data on SG1 depth
Adele1_on_SG1_depth = tune(Adele1_on_Omar1_depth, cbind(Omar1_standardized$Omar1_scaled.Center_win[al_o1_sg1_ap2$index1s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_o1_sg1_ap2$index2s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 20), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (Adele-1)")
lines(Adele1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Adele1_originalGR_on_SG1_depth = data.frame(Adele1_on_SG1_depth$X1, Adele1_interpolated[1151:3184,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (Adele-1)")
lines(Adele1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Adele1 

SG1Age_on_Adele1_depth = tune(Adele1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Adele1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Adele1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Adele1_depth) <- new_column_names
SG1Age_on_Adele1_depth[,1] = SG1Age_on_Adele1_depth[,1] * 1000
write.csv(SG1Age_on_Adele1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Adele 1.csv", row.names = FALSE)
