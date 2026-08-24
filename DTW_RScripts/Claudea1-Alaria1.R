install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Claudea 1 and Alaria 1 datasets

Claudea1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Claudea 1.csv", header=TRUE, stringsAsFactors=FALSE)
Claudea1=Claudea1[c(1:5764),] # Oligocene-Miocene
head(Claudea1)
plot(Claudea1, type="l", xlim = c(500, 1700), ylim = c(0, 60))

Alaria1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Alaria 1.csv", header=TRUE, stringsAsFactors=FALSE)
Alaria1=Alaria1[c(1:5995),] # Oligocene-Miocene
head(Alaria1)
plot(Alaria1, type="l", xlim = c(450, 1650), ylim = c(0, 60))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Claudea1_interpolated <- linterp(Claudea1, dt = 0.2, genplot = F)
Alaria1_interpolated <- linterp(Alaria1, dt = 0.2, genplot = F)

# Scaling the data
Clmean = Gmean(Claudea1_interpolated$GR)
Clstd = Gsd(Claudea1_interpolated$GR)
Claudea1_scaled = (Claudea1_interpolated$GR - Clmean)/Clstd
Claudea1_rescaled = data.frame(Claudea1_interpolated$DEPT, Claudea1_scaled)

Amean = Gmean(Alaria1_interpolated$GR)
Astd = Gsd(Alaria1_interpolated$GR)
Alaria1_scaled = (Alaria1_interpolated$GR - Amean)/Astd
Alaria1_rescaled = data.frame(Alaria1_interpolated$DEPT, Alaria1_scaled)

# Resampling the data using moving window statistics
Claudea1_scaled = mwStats(Claudea1_rescaled, cols = 2, win=3, ends = T)
Claudea1_standardized = data.frame(Claudea1_scaled$Center_win, Claudea1_scaled$Average)

Alaria1_scaled = mwStats(Alaria1_rescaled, cols = 2, win=3, ends = T)
Alaria1_standardized = data.frame(Alaria1_scaled$Center_win, Alaria1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Claudea1_standardized, type="l", xlim = c(500, 1700), ylim = c(-20, 20), xlab = "Claudea1 Resampled Depth", ylab = "Normalized GR")
plot(Alaria1_standardized, type="l", xlim = c(450, 1650), ylim = c(-10, 20), xlab = "Alaria 1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_a1_cl1_ap1 <- dtw(Alaria1_standardized$Alaria1_scaled.Average, Claudea1_standardized$Claudea1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_a1_cl1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Alaria1_on_Claudea1_depth = tune(Alaria1_standardized, cbind(Alaria1_standardized$Alaria1_scaled.Center_win[al_a1_cl1_ap1$index1s], Claudea1_standardized$Claudea1_scaled.Center_win[al_a1_cl1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Claudea1_standardized, type = "l", ylim = c(-20, 20), xlim = c(500, 1700), xlab = "Claudea1 Resampled Depth", ylab = "Normalized GR")
lines(Alaria1_on_Claudea1_depth, col = "red")

# DTW Distance
al_a1_cl1_ap1$normalizedDistance
al_a1_cl1_ap1$distance

# Tuning Alaria1 data on SG1 depth
Alaria1_on_Brontosaurus1_depth = tune(Alaria1_on_Claudea1_depth, cbind(Claudea1_standardized$Claudea1_scaled.Center_win[al_c1_b1_ap1$index1s], Brontosaurus1_standardized$Brontosaurus1_scaled.Center_win[al_c1_b1_ap1$index2s]), extrapolate = F)

# Tuning Alaria1 data on SG1 depth
Alaria1_on_SG1_depth = tune(Alaria1_on_Brontosaurus1_depth, cbind(Brontosaurus11_standardized$Brontosaurus11_scaled.Center_win[al_sg1_b1_ap2$index2s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_sg1_b1_ap2$index1s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (Alaria-1)")
lines(Alaria1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Alaria1_originalGR_on_SG1_depth = data.frame(Alaria1_on_SG1_depth$X1, Alaria1_interpolated[890:6001,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (Alaria-1)")
lines(Alaria1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Alaria1 

SG1Age_on_Alaria1_depth = tune(Alaria1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Alaria1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Alaria-1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Alaria1_depth) <- new_column_names
SG1Age_on_Alaria1_depth[,1] = SG1Age_on_Alaria1_depth[,1] * 1000
write.csv(SG1Age_on_Alaria1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Alaria 1.csv", row.names = FALSE)

Alaria1_age_depth <- approx(x = Alaria1_interpolated$GR,
                              y = Alaria1_interpolated$DEPT,
                              xout = SG1Age_on_Alaria1_depth$GR,
                              rule = 1)$y

Alaria1_agemodel = data.frame(Alaria1_age_depth, SG1Age_on_Alaria1_depth$AGE)
Alaria1_agemodel <- na.omit(Alaria1_agemodel)
new_column_names1 <- c("Depth", "Age")
colnames(Alaria1_agemodel) <- new_column_names1
plot(Alaria1_agemodel)
write.csv(Alaria1_agemodel, file = "RScripts&Data/Sites Data_Age-NGR/Alaria1_DepthAge.csv", row.names = FALSE)
