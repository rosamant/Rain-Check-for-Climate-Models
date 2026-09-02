install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Crown 1 and Pharos1 datasets

Crown1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Crown 1.csv", header=TRUE, stringsAsFactors=FALSE)
Crown1=Crown1[c(1:8673),] # Oligocene-Miocene
head(Crown1)
plot(Crown1, type="l", xlim = c(450, 2250), ylim = c(0, 90))

Pharos1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Pharos1.csv", header=TRUE, stringsAsFactors=FALSE)
Pharos1=Pharos1[c(1:8816),]
head(Pharos1)
plot(Pharos1, type='l', xlim= c(500,2250), ylim = c(0,90))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Crown1_interpolated <- linterp(Crown1, dt = 0.2, genplot = F)
Pharos1_interpolated <- linterp(Pharos1, dt = 0.2, genplot = F)

# Scaling the data
Crmean = Gmean(Crown1_interpolated$GR)
Crstd = Gsd(Crown1_interpolated$GR)
Crown1_scaled = (Crown1_interpolated$GR - Crmean)/Crstd
Crown1_rescaled = data.frame(Crown1_interpolated$DEPT, Crown1_scaled)

Phmean = Gmean(Pharos1_interpolated$GR)
Phstd = Gsd(Pharos1_interpolated$GR)
Pharos1_scaled = (Pharos1_interpolated$GR - Phmean)/Phstd
Pharos1_rescaled = data.frame(Pharos1_interpolated$DEPT, Pharos1_scaled)

# Resampling the data using moving window statistics
Crown1_scaled = mwStats(Crown1_rescaled, cols = 2, win=3, ends = T)
Crown1_standardized = data.frame(Crown1_scaled$Center_win, Crown1_scaled$Average)

Pharos1_scaled = mwStats(Pharos1_rescaled, cols = 2, win=3, ends = T)
Pharos1_standardized = data.frame(Pharos1_scaled$Center_win, Pharos1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Crown1_standardized, type="l", xlim = c(450, 2250), ylim = c(-20, 20), xlab = "Crown1 Resampled Depth", ylab = "Normalized GR")
plot(Pharos1_standardized, type="l", xlim = c(500, 2250), ylim = c(-20, 40), xlab = "Pharos1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_ph1_cr1_ap1 <- dtw(Pharos1_standardized$Pharos1_scaled.Average, Crown1_standardized$Crown1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_ph1_cr1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Pharos1_on_Crown1_depth = tune(Pharos1_standardized, cbind(Pharos1_standardized$Pharos1_scaled.Center_win[al_ph1_cr1_ap1$index1s], Crown1_standardized$Crown1_scaled.Center_win[al_ph1_cr1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Crown1_standardized, type = "l", ylim = c(-20, 40), xlim = c(450, 2250), xlab = "Crown1 Resampled Depth", ylab = "Normalized GR")
lines(Pharos1_on_Crown1_depth, col = "red")

# DTW Distance

al_ph1_cr1_ap1$normalizedDistance
al_ph1_cr1_ap1$distance

# Tuning Pharos1 data on Caswell1 depth
Pharos1_on_Caswell1_depth = tune(Pharos1_on_Crown1_depth, cbind(Crown1_standardized$Crown1_scaled.Center_win[al_cr1_c1_ap2$index1s], Caswell1_standardized$Caswell1_scaled.Center_win[al_cr1_c1_ap2$index2s]), extrapolate = F)

# Tuning Pharos1 data on Calliance2 depth
Pharos1_on_Calliance2_depth = tune(Pharos1_on_Caswell1_depth, cbind(Caswell1_standardized$Caswell1_scaled.Center_win[al_c2_c1_ap1$index2s], Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_c1_ap1$index1s]), extrapolate = F)

# Tuning Pharos1 data on Omar1 depth
Pharos1_on_Omar1_depth = tune(Pharos1_on_Calliance2_depth, cbind(Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_o1_ap2$index1s], Omar1_standardized$Omar1_scaled.Center_win[al_c2_o1_ap2$index2s]), extrapolate = F)

# Tuning Pharos1 data on SG1 depth
Pharos1_on_SG1_depth = tune(Pharos1_on_Omar1_depth, cbind(Omar1_standardized$Omar1_scaled.Center_win[al_o1_sg1_ap2$index1s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_o1_sg1_ap2$index2s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (Pharos-1)")
lines(Pharos1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Pharos1_originalGR_on_SG1_depth = data.frame(Pharos1_on_SG1_depth$X1, Pharos1_interpolated[1262:8536,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (Pharos-1)")
lines(Pharos1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Pharos1 

SG1Age_on_Pharos1_depth = tune(Pharos1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Pharos1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Pharos1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Pharos1_depth) <- new_column_names
SG1Age_on_Pharos1_depth[,1] = SG1Age_on_Pharos1_depth[,1] * 1000
write.csv(SG1Age_on_Pharos1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Pharos 1.csv", row.names = FALSE)
