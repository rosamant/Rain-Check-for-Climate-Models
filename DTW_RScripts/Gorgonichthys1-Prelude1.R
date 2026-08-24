install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import  Gorgonichthys 1 and Prelude 1 datasets

Gorgonichthys1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Gorgonichthys1.csv", header=TRUE, stringsAsFactors=FALSE)
Gorgonichthys1=Gorgonichthys1[c(168:6001),] # Oligocene-Miocene
head(Gorgonichthys1)
plot(Gorgonichthys1, type="l", xlim = c(250, 1200), ylim = c(0, 30))

Prelude1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Prelude1.csv", header=TRUE, stringsAsFactors=FALSE)
Prelude1=Prelude1[c(200:5314),] # Oligocene-Miocene
head(Prelude1)
plot(Prelude1, type="l", xlim = c(220, 1050), ylim = c(0, 40))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Gorgonichthys1_interpolated <- linterp(Gorgonichthys1, dt = 0.2, genplot = F)
Prelude1_interpolated <- linterp(Prelude1, dt = 0.2, genplot = F)

# Scaling the data
Gmean = Gmean(Gorgonichthys1_interpolated$GR)
Gstd = Gsd(Gorgonichthys1_interpolated$GR)
Gorgonichthys1_scaled = (Gorgonichthys1_interpolated$GR - Gmean)/Gstd
Gorgonichthys1_rescaled = data.frame(Gorgonichthys1_interpolated$DEPT, Gorgonichthys1_scaled)

Pmean = Gmean(Prelude1_interpolated$GR)
Pstd = Gsd(Prelude1_interpolated$GR)
Prelude1_scaled = (Prelude1_interpolated$GR - Pmean)/Pstd
Prelude1_rescaled = data.frame(Prelude1_interpolated$DEPT, Prelude1_scaled)

# Resampling the data using moving window statistics
Gorgonichthys1_scaled = mwStats(Gorgonichthys1_rescaled, cols = 2, win=3, ends = T)
Gorgonichthys1_standardized = data.frame(Gorgonichthys1_scaled$Center_win, Gorgonichthys1_scaled$Average)

Prelude1_scaled = mwStats(Prelude1_rescaled, cols = 2, win=3, ends = T)
Prelude1_standardized = data.frame(Prelude1_scaled$Center_win, Prelude1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Gorgonichthys1_standardized, type="l", xlim = c(250, 1200), ylim = c(-20, 20), xlab = "Gorgonichthys1 Resampled Depth", ylab = "Normalized GR")
plot(Prelude1_standardized, type="l", xlim = c(220, 1050), ylim = c(-10, 10), xlab = "Prelude 1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_p1_g1_ap1 <- dtw(Prelude1_standardized$Prelude1_scaled.Average, Gorgonichthys1_standardized$Gorgonichthys1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_p1_g1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Prelude1_on_Gorgonichthys1_depth = tune(Prelude1_standardized, cbind(Prelude1_standardized$Prelude1_scaled.Center_win[al_p1_g1_ap1$index1s], Gorgonichthys1_standardized$Gorgonichthys1_scaled.Center_win[al_p1_g1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Gorgonichthys1_standardized, type = "l", ylim = c(-15, 15), xlim = c(250, 1200), xlab = "Gorgonichthys1 Resampled Depth", ylab = "Normalized GR")
lines(Prelude1_on_Gorgonichthys1_depth, col = "red")

# DTW Distance
al_p1_g1_ap1$normalizedDistance
al_p1_g1_ap1$distance

# Tuning Prelude1 data on Walkley1 depth
Prelude1_on_Walkley1_depth = tune(Prelude1_on_Gorgonichthys1_depth, cbind(Gorgonichthys1_standardized$Gorgonichthys1_scaled.Center_win[al_g1_w1_ap2$index1s], Walkley1_standardized$Walkley1_scaled.Center_win[al_g1_w1_ap2$index2s]), extrapolate = F)

# Tuning Prelude1 data on Caswell1 depth
Prelude1_on_Caswell1_depth = tune(Prelude1_on_Walkley1_depth, cbind(Walkley1_standardized$Walkley1_scaled.Center_win[al_w1_c1_ap1$index1s], Caswell1_standardized$Caswell1_scaled.Center_win[al_w1_c1_ap1$index2s]), extrapolate = F)

# Tuning Prelude1 data on Calliance2 depth
Prelude1_on_Calliance2_depth = tune(Prelude1_on_Caswell1_depth, cbind(Caswell1_standardized$Caswell1_scaled.Center_win[al_c2_c1_ap1$index2s], Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_c1_ap1$index1s]), extrapolate = F)

# Tuning Prelude1 data on Omar1 depth
Prelude1_on_Omar1_depth = tune(Prelude1_on_Calliance2_depth, cbind(Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_o1_ap2$index1s], Omar1_standardized$Omar1_scaled.Center_win[al_c2_o1_ap2$index2s]), extrapolate = F)

# Tuning Prelude1 data on SG1 depth
Prelude1_on_SG1_depth = tune(Prelude1_on_Omar1_depth, cbind(Omar1_standardized$Omar1_scaled.Center_win[al_o1_sg1_ap2$index1s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_o1_sg1_ap2$index2s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 20), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (Prelude-1)")
lines(Prelude1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Prelude1_originalGR_on_SG1_depth = data.frame(Prelude1_on_SG1_depth$X1, Prelude1_interpolated[1131:3892,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (Prelude-1)")
lines(Prelude1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Prelude1 

SG1Age_on_Prelude1_depth = tune(Prelude1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Prelude1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Prelude1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Prelude1_depth) <- new_column_names
SG1Age_on_Prelude1_depth[,1] = SG1Age_on_Prelude1_depth[,1] * 1000
write.csv(SG1Age_on_Prelude1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Prelude 1.csv", row.names = FALSE)
