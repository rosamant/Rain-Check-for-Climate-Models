install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Crown 1 and Poseidon1 datasets

Crown1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Crown 1.csv", header=TRUE, stringsAsFactors=FALSE)
Crown1=Crown1[c(1:8673),] # Oligocene-Miocene
head(Crown1)
plot(Crown1, type="l", xlim = c(450, 2250), ylim = c(0, 90))

Poseidon1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Poseidon 1.csv", header=TRUE, stringsAsFactors=FALSE)
Poseidon1=Poseidon1[c(1:9601),]
head(Poseidon1)
plot(Poseidon1, type='l', xlim= c(500,2450), ylim = c(0,90))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Crown1_interpolated <- linterp(Crown1, dt = 0.2, genplot = F)
Poseidon1_interpolated <- linterp(Poseidon1, dt = 0.2, genplot = F)

# Scaling the data
Crmean = Gmean(Crown1_interpolated$GR)
Crstd = Gsd(Crown1_interpolated$GR)
Crown1_scaled = (Crown1_interpolated$GR - Crmean)/Crstd
Crown1_rescaled = data.frame(Crown1_interpolated$DEPT, Crown1_scaled)

Pomean = Gmean(Poseidon1_interpolated$GR)
Postd = Gsd(Poseidon1_interpolated$GR)
Poseidon1_scaled = (Poseidon1_interpolated$GR - Pomean)/Postd
Poseidon1_rescaled = data.frame(Poseidon1_interpolated$DEPT, Poseidon1_scaled)

# Resampling the data using moving window statistics
Crown1_scaled = mwStats(Crown1_rescaled, cols = 2, win=3, ends = T)
Crown1_standardized = data.frame(Crown1_scaled$Center_win, Crown1_scaled$Average)

Poseidon1_scaled = mwStats(Poseidon1_rescaled, cols = 2, win=3, ends = T)
Poseidon1_standardized = data.frame(Poseidon1_scaled$Center_win, Poseidon1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Crown1_standardized, type="l", xlim = c(450, 2250), ylim = c(-20, 20), xlab = "Crown1 Resampled Depth", ylab = "Normalized GR")
plot(Poseidon1_standardized, type="l", xlim = c(500, 2450), ylim = c(-20, 40), xlab = "Poseidon1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_po1_cr1_ap1 <- dtw(Poseidon1_standardized$Poseidon1_scaled.Average, Crown1_standardized$Crown1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_po1_cr1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Poseidon1_on_Crown1_depth = tune(Poseidon1_standardized, cbind(Poseidon1_standardized$Poseidon1_scaled.Center_win[al_po1_cr1_ap1$index1s], Crown1_standardized$Crown1_scaled.Center_win[al_po1_cr1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Crown1_standardized, type = "l", ylim = c(-20, 40), xlim = c(450, 2250), xlab = "Crown1 Resampled Depth", ylab = "Normalized GR")
lines(Poseidon1_on_Crown1_depth, col = "red")

# DTW Distance

al_po1_cr1_ap1$normalizedDistance
al_po1_cr1_ap1$distance

# Tuning Poseidon1 data on Caswell1 depth
Poseidon1_on_Caswell1_depth = tune(Poseidon1_on_Crown1_depth, cbind(Crown1_standardized$Crown1_scaled.Center_win[al_cr1_c1_ap2$index1s], Caswell1_standardized$Caswell1_scaled.Center_win[al_cr1_c1_ap2$index2s]), extrapolate = F)

# Tuning Poseidon1 data on Calliance2 depth
Poseidon1_on_Calliance2_depth = tune(Poseidon1_on_Caswell1_depth, cbind(Caswell1_standardized$Caswell1_scaled.Center_win[al_c2_c1_ap1$index2s], Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_c1_ap1$index1s]), extrapolate = F)

# Tuning Poseidon1 data on Omar1 depth
Poseidon1_on_Omar1_depth = tune(Poseidon1_on_Calliance2_depth, cbind(Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_o1_ap2$index1s], Omar1_standardized$Omar1_scaled.Center_win[al_c2_o1_ap2$index2s]), extrapolate = F)

# Tuning Poseidon1 data on SG1 depth
Poseidon1_on_SG1_depth = tune(Poseidon1_on_Omar1_depth, cbind(Omar1_standardized$Omar1_scaled.Center_win[al_o1_sg1_ap2$index1s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_o1_sg1_ap2$index2s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (Poseidon-1)")
lines(Poseidon1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Poseidon1_originalGR_on_SG1_depth = data.frame(Poseidon1_on_SG1_depth$X1, Poseidon1_interpolated[1187:9380,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (Poseidon-1)")
lines(Poseidon1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Poseidon1 

SG1Age_on_Poseidon1_depth = tune(Poseidon1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Poseidon1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Poseidon1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Poseidon1_depth) <- new_column_names
SG1Age_on_Poseidon1_depth[,1] = SG1Age_on_Poseidon1_depth[,1] * 1000
write.csv(SG1Age_on_Poseidon1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Poseidon 1.csv", row.names = FALSE)
