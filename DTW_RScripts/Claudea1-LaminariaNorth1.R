install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Claudea 1 and LaminariaNorth 1 datasets

Claudea1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Claudea 1.csv", header=TRUE, stringsAsFactors=FALSE)
Claudea1=Claudea1[c(1:5764),] # Oligocene-Miocene
head(Claudea1)
plot(Claudea1, type="l", xlim = c(500, 1700), ylim = c(0, 60))

LaminariaNorth1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Laminaria North 1.csv", header=TRUE, stringsAsFactors=FALSE)
LaminariaNorth1=LaminariaNorth1[c(1:6309),] # Oligocene-Miocene
head(LaminariaNorth1)
plot(LaminariaNorth1, type="l", xlim = c(400, 1700), ylim = c(0, 60))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Claudea1_interpolated <- linterp(Claudea1, dt = 0.2, genplot = F)
LaminariaNorth1_interpolated <- linterp(LaminariaNorth1, dt = 0.2, genplot = F)

# Scaling the data
Clmean = Gmean(Claudea1_interpolated$GR)
Clstd = Gsd(Claudea1_interpolated$GR)
Claudea1_scaled = (Claudea1_interpolated$GR - Clmean)/Clstd
Claudea1_rescaled = data.frame(Claudea1_interpolated$DEPT, Claudea1_scaled)

Lnmean = Gmean(LaminariaNorth1_interpolated$GR)
Lnstd = Gsd(LaminariaNorth1_interpolated$GR)
LaminariaNorth1_scaled = (LaminariaNorth1_interpolated$GR - Lnmean)/Lnstd
LaminariaNorth1_rescaled = data.frame(LaminariaNorth1_interpolated$DEPT, LaminariaNorth1_scaled)

# Resampling the data using moving window statistics
Claudea1_scaled = mwStats(Claudea1_rescaled, cols = 2, win=3, ends = T)
Claudea1_standardized = data.frame(Claudea1_scaled$Center_win, Claudea1_scaled$Average)

LaminariaNorth1_scaled = mwStats(LaminariaNorth1_rescaled, cols = 2, win=3, ends = T)
LaminariaNorth1_standardized = data.frame(LaminariaNorth1_scaled$Center_win, LaminariaNorth1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Claudea1_standardized, type="l", xlim = c(500, 1700), ylim = c(-20, 20), xlab = "Claudea1 Resampled Depth", ylab = "Normalized GR")
plot(LaminariaNorth1_standardized, type="l", xlim = c(400, 1700), ylim = c(-10, 20), xlab = "LaminariaNorth 1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_ln1_cl1_ap1 <- dtw(LaminariaNorth1_standardized$LaminariaNorth1_scaled.Average, Claudea1_standardized$Claudea1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_ln1_cl1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
LaminariaNorth1_on_Claudea1_depth = tune(LaminariaNorth1_standardized, cbind(LaminariaNorth1_standardized$LaminariaNorth1_scaled.Center_win[al_ln1_cl1_ap1$index1s], Claudea1_standardized$Claudea1_scaled.Center_win[al_ln1_cl1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Claudea1_standardized, type = "l", ylim = c(-20, 20), xlim = c(500, 1700), xlab = "Claudea1 Resampled Depth", ylab = "Normalized GR")
lines(LaminariaNorth1_on_Claudea1_depth, col = "red")

# DTW Distance
al_ln1_cl1_ap1$normalizedDistance
al_ln1_cl1_ap1$distance

# Tuning LaminariaNorth1 data on SG1 depth
LaminariaNorth1_on_Brontosaurus1_depth = tune(LaminariaNorth1_on_Claudea1_depth, cbind(Claudea1_standardized$Claudea1_scaled.Center_win[al_c1_b1_ap1$index1s], Brontosaurus1_standardized$Brontosaurus1_scaled.Center_win[al_c1_b1_ap1$index2s]), extrapolate = F)

# Tuning LaminariaNorth1 data on SG1 depth
LaminariaNorth1_on_SG1_depth = tune(LaminariaNorth1_on_Brontosaurus1_depth, cbind(Brontosaurus11_standardized$Brontosaurus11_scaled.Center_win[al_sg1_b1_ap2$index2s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_sg1_b1_ap2$index1s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (LaminariaNorth-1)")
lines(LaminariaNorth1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
LaminariaNorth1_originalGR_on_SG1_depth = data.frame(LaminariaNorth1_on_SG1_depth$X1, LaminariaNorth1_interpolated[958:6539,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (LaminariaNorth-1)")
lines(LaminariaNorth1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to LaminariaNorth1 

SG1Age_on_LaminariaNorth1_depth = tune(LaminariaNorth1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_LaminariaNorth1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "LaminariaNorth-1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_LaminariaNorth1_depth) <- new_column_names
SG1Age_on_LaminariaNorth1_depth[,1] = SG1Age_on_LaminariaNorth1_depth[,1] * 1000
write.csv(SG1Age_on_LaminariaNorth1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Laminaria North 1.csv", row.names = FALSE)
