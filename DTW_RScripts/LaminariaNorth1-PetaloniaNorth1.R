install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import LaminariaNorth 1 and PetaloniaNorth 1 datasets

LaminariaNorth1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Laminaria North 1.csv", header=TRUE, stringsAsFactors=FALSE)
LaminariaNorth1=LaminariaNorth1[c(1:6309),] # Oligocene-Miocene
head(LaminariaNorth1)
plot(LaminariaNorth1, type="l", xlim = c(400, 1700), ylim = c(0, 60))

PetaloniaNorth1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Petalonia North 1.csv", header=TRUE, stringsAsFactors=FALSE)
PetaloniaNorth1=PetaloniaNorth1[c(1:13971),] # Oligocene-Miocene
head(PetaloniaNorth1)
plot(PetaloniaNorth1, type="l", xlim = c(400, 1800), ylim = c(0, 60))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
LaminariaNorth1_interpolated <- linterp(LaminariaNorth1, dt = 0.2, genplot = F)
PetaloniaNorth1_interpolated <- linterp(PetaloniaNorth1, dt = 0.2, genplot = F)

# Scaling the data
Lnmean = Gmean(LaminariaNorth1_interpolated$GR)
Lnstd = Gsd(LaminariaNorth1_interpolated$GR)
LaminariaNorth1_scaled = (LaminariaNorth1_interpolated$GR - Lnmean)/Lnstd
LaminariaNorth1_rescaled = data.frame(LaminariaNorth1_interpolated$DEPT, LaminariaNorth1_scaled)

Pnmean = Gmean(PetaloniaNorth1_interpolated$GR)
Pnstd = Gsd(PetaloniaNorth1_interpolated$GR)
PetaloniaNorth1_scaled = (PetaloniaNorth1_interpolated$GR - Pnmean)/Pnstd
PetaloniaNorth1_rescaled = data.frame(PetaloniaNorth1_interpolated$DEPT, PetaloniaNorth1_scaled)

# Resampling the data using moving window statistics
LaminariaNorth1_scaled = mwStats(LaminariaNorth1_rescaled, cols = 2, win=3, ends = T)
LaminariaNorth1_standardized = data.frame(LaminariaNorth1_scaled$Center_win, LaminariaNorth1_scaled$Average)

PetaloniaNorth1_scaled = mwStats(PetaloniaNorth1_rescaled, cols = 2, win=3, ends = T)
PetaloniaNorth1_standardized = data.frame(PetaloniaNorth1_scaled$Center_win, PetaloniaNorth1_scaled$Average)

# Plotting the rescaled and resampled data
plot(LaminariaNorth1_standardized, type="l", xlim = c(400, 1700), ylim = c(-20, 20), xlab = "LaminariaNorth1 Resampled Depth", ylab = "Normalized GR")
plot(PetaloniaNorth1_standardized, type="l", xlim = c(400, 1800), ylim = c(-10, 20), xlab = "PetaloniaNorth 1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_pn1_ln1_ap1 <- dtw(PetaloniaNorth1_standardized$PetaloniaNorth1_scaled.Average, LaminariaNorth1_standardized$LaminariaNorth1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_pn1_ln1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
PetaloniaNorth1_on_LaminariaNorth1_depth = tune(PetaloniaNorth1_standardized, cbind(PetaloniaNorth1_standardized$PetaloniaNorth1_scaled.Center_win[al_pn1_ln1_ap1$index1s], LaminariaNorth1_standardized$LaminariaNorth1_scaled.Center_win[al_pn1_ln1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(LaminariaNorth1_standardized, type = "l", ylim = c(-20, 20), xlim = c(400, 1700), xlab = "LaminariaNorth1 Resampled Depth", ylab = "Normalized GR")
lines(PetaloniaNorth1_on_LaminariaNorth1_depth, col = "red")

# DTW Distance
al_pn1_ln1_ap1$normalizedDistance
al_pn1_ln1_ap1$distance

# Tuning PetaloniaNorth1 data on SG1 depth
PetaloniaNorth1_on_Claudea1_depth = tune(PetaloniaNorth1_on_LaminariaNorth1_depth, cbind(LaminariaNorth1_standardized$LaminariaNorth1_scaled.Center_win[al_ln1_cl1_ap1$index1s], Claudea1_standardized$Claudea1_scaled.Center_win[al_ln1_cl1_ap1$index2s]), extrapolate = F)

# Tuning PetaloniaNorth1 data on SG1 depth
PetaloniaNorth1_on_Brontosaurus1_depth = tune(PetaloniaNorth1_on_Claudea1_depth, cbind(Claudea1_standardized$Claudea1_scaled.Center_win[al_c1_b1_ap1$index1s], Brontosaurus1_standardized$Brontosaurus1_scaled.Center_win[al_c1_b1_ap1$index2s]), extrapolate = F)

# Tuning PetaloniaNorth1 data on SG1 depth
PetaloniaNorth1_on_SG1_depth = tune(PetaloniaNorth1_on_Brontosaurus1_depth, cbind(Brontosaurus11_standardized$Brontosaurus11_scaled.Center_win[al_sg1_b1_ap2$index2s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_sg1_b1_ap2$index1s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (PetaloniaNorth-1)")
lines(PetaloniaNorth1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
PetaloniaNorth1_originalGR_on_SG1_depth = data.frame(PetaloniaNorth1_on_SG1_depth$X1, PetaloniaNorth1_interpolated[1274:6986,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (PetaloniaNorth-1)")
lines(PetaloniaNorth1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to PetaloniaNorth1 

SG1Age_on_PetaloniaNorth1_depth = tune(PetaloniaNorth1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_PetaloniaNorth1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "PetaloniaNorth-1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_PetaloniaNorth1_depth) <- new_column_names
SG1Age_on_PetaloniaNorth1_depth[,1] = SG1Age_on_PetaloniaNorth1_depth[,1] * 1000
write.csv(SG1Age_on_PetaloniaNorth1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Petalonia North 1.csv", row.names = FALSE)
