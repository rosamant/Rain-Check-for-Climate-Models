install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Caswell1 and Calliance2 datasets

Caswell1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Caswell1.csv", header=TRUE, stringsAsFactors=FALSE)
Caswell1=Caswell1[c(1:6434),] # Oligocene-Miocene
head(Caswell1)
plot(Caswell1, type="l", xlim = c(350, 1650), ylim = c(0, 50))

Calliance2 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Calliance 2.csv", header=TRUE, stringsAsFactors=FALSE)
Calliance2=Calliance2[c(1:14202),] # Oligocene-Miocene
head(Calliance2)
plot(Calliance2, type="l", xlim = c(600, 2000), ylim = c(0, 80))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Caswell1_interpolated <- linterp(Caswell1, dt = 0.2, genplot = F)
Calliance2_interpolated <- linterp(Calliance2, dt = 0.2, genplot = F)

# Scaling the data
Cmean = Gmean(Caswell1_interpolated$GR)
Cstd = Gsd(Caswell1_interpolated$GR)
Caswell1_scaled = (Caswell1_interpolated$GR - Cmean)/Cstd
Caswell1_rescaled = data.frame(Caswell1_interpolated$DEPT, Caswell1_scaled)

C2mean = Gmean(Calliance2_interpolated$GR)
C2std = Gsd(Calliance2_interpolated$GR)
Calliance2_scaled = (Calliance2_interpolated$GR - C2mean)/C2std
Calliance2_rescaled = data.frame(Calliance2_interpolated$DEPT, Calliance2_scaled)

# Resampling the data using moving window statistics
Caswell1_scaled = mwStats(Caswell1_rescaled, cols = 2, win=3, ends = T)
Caswell1_standardized = data.frame(Caswell1_scaled$Center_win, Caswell1_scaled$Average)

Calliance2_scaled = mwStats(Calliance2_rescaled, cols = 2, win=3, ends = T)
Calliance2_standardized = data.frame(Calliance2_scaled$Center_win, Calliance2_scaled$Average)

# Plotting the rescaled and resampled data
plot(Caswell1_standardized, type="l", xlim = c(350, 1700), ylim = c(-20, 20), xlab = "Caswell1 Resampled Depth", ylab = "Normalized GR")
plot(Calliance2_standardized, type="l", xlim = c(600, 2000), ylim = c(-20, 30), xlab = "Calliance2 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_c2_c1_ap1 <- dtw(Calliance2_standardized$Calliance2_scaled.Average, Caswell1_standardized$Caswell1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_c2_c1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Calliance2_on_Caswell1_depth = tune(Calliance2_standardized, cbind(Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_c1_ap1$index1s], Caswell1_standardized$Caswell1_scaled.Center_win[al_c2_c1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Caswell1_standardized, type = "l", ylim = c(-20, 30), xlim = c(350, 1600), xlab = "Calliance 2 Resampled Depth", ylab = "Normalized GR")
lines(Calliance2_on_Caswell1_depth, col = "red")

# DTW Distance

al_c2_c1_ap1$normalizedDistance
al_c2_c1_ap1$distance

# Tuning Caswell1 data on Calliance2 depth
Caswell1_on_Calliance2_depth = tune(Caswell1_standardized, cbind(Caswell1_standardized$Caswell1_scaled.Center_win[al_c2_c1_ap1$index2s], Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_c1_ap1$index1s]), extrapolate = F)

# Tuning Caswell1 data on Omar1 depth
Caswell1_on_Omar1_depth = tune(Caswell1_on_Calliance2_depth, cbind(Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_o1_ap2$index1s], Omar1_standardized$Omar1_scaled.Center_win[al_c2_o1_ap2$index2s]), extrapolate = F)

# Tuning Caswell1 data on SG1 depth
Caswell1_on_SG1_depth = tune(Caswell1_on_Omar1_depth, cbind(Omar1_standardized$Omar1_scaled.Center_win[al_o1_sg1_ap2$index1s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_o1_sg1_ap2$index2s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (Caswell-1)")
lines(Caswell1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Caswell1_originalGR_on_SG1_depth = data.frame(Caswell1_on_SG1_depth$X1, Caswell1_interpolated[1552:6372,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (Caswell-1)")
lines(Caswell1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Caswell1 

SG1Age_on_Caswell1_depth = tune(Caswell1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Caswell1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Caswell1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Caswell1_depth) <- new_column_names
SG1Age_on_Caswell1_depth[,1] = SG1Age_on_Caswell1_depth[,1] * 1000
write.csv(SG1Age_on_Caswell1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Caswell 1.csv", row.names = FALSE)
