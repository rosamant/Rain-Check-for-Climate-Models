install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Brontosaurus 1 and Pokolbin 1 datasets

Brontosaurus1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Brontosaurus 1.csv", header=TRUE, stringsAsFactors=FALSE)
Brontosaurus1=Brontosaurus1[c(1:4527),] # Oligocene-Miocene
head(Brontosaurus1)
plot(Brontosaurus1, type="l", xlim = c(150, 1050), ylim = c(20, 60))

Pokolbin1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Pokolbin 1.csv", header=TRUE, stringsAsFactors=FALSE)
Pokolbin1=Pokolbin1[c(50:3987),] # Oligocene-Miocene
head(Pokolbin1)
plot(Pokolbin1, type="l", xlim = c(230, 1050), ylim = c(0, 60))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Brontosaurus1_interpolated <- linterp(Brontosaurus1, dt = 0.2, genplot = F)
Pokolbin1_interpolated <- linterp(Pokolbin1, dt = 0.2, genplot = F)

# Scaling the data
Bmean = Gmean(Brontosaurus1_interpolated$GR)
Bstd = Gsd(Brontosaurus1_interpolated$GR)
Brontosaurus1_scaled = (Brontosaurus1_interpolated$GR - Bmean)/Bstd
Brontosaurus1_rescaled = data.frame(Brontosaurus1_interpolated$DEPT, Brontosaurus1_scaled)

Pbmean = Gmean(Pokolbin1_interpolated$GR)
Pbstd = Gsd(Pokolbin1_interpolated$GR)
Pokolbin1_scaled = (Pokolbin1_interpolated$GR - Pbmean)/Pbstd
Pokolbin1_rescaled = data.frame(Pokolbin1_interpolated$DEPT, Pokolbin1_scaled)

# Resampling the data using moving window statistics
Brontosaurus1_scaled = mwStats(Brontosaurus1_rescaled, cols = 2, win=3, ends = T)
Brontosaurus1_standardized = data.frame(Brontosaurus1_scaled$Center_win, Brontosaurus1_scaled$Average)

Pokolbin1_scaled = mwStats(Pokolbin1_rescaled, cols = 2, win=3, ends = T)
Pokolbin1_standardized = data.frame(Pokolbin1_scaled$Center_win, Pokolbin1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Brontosaurus1_standardized, type="l", xlim = c(150, 1050), ylim = c(-20, 20), xlab = "Brontosaurus1 Resampled Depth", ylab = "Normalized GR")
plot(Pokolbin1_standardized, type="l", xlim = c(230, 1050), ylim = c(-10, 20), xlab = "Pokolbin 2 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_pb1_b1_ap1 <- dtw(Pokolbin1_standardized$Pokolbin1_scaled.Average, Brontosaurus1_standardized$Brontosaurus1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = F, open.end = F))
plot(al_pb1_b1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Pokolbin1_on_Brontosaurus1_depth = tune(Pokolbin1_standardized, cbind(Pokolbin1_standardized$Pokolbin1_scaled.Center_win[al_pb1_b1_ap1$index1s], Brontosaurus1_standardized$Brontosaurus1_scaled.Center_win[al_pb1_b1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Brontosaurus1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1050), xlab = "Brontosaurus1 Resampled Depth", ylab = "Normalized GR")
lines(Pokolbin1_on_Brontosaurus1_depth, col = "red")

# DTW Distance
al_h1_b1_ap1$normalizedDistance
al_h1_b1_ap1$distance

# Tuning Pokolbin1 data on SG1 depth
Pokolbin1_on_SG1_depth = tune(Pokolbin1_on_Brontosaurus1_depth, cbind(Brontosaurus11_standardized$Brontosaurus11_scaled.Center_win[al_sg1_b1_ap2$index2s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_sg1_b1_ap2$index1s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (Pokolbin-1)")
lines(Pokolbin1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Pokolbin1_originalGR_on_SG1_depth = data.frame(Pokolbin1_on_SG1_depth$X1, Pokolbin1_interpolated[769:3977,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (Pokolbin-1)")
lines(Pokolbin1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Pokolbin1 

SG1Age_on_Pokolbin1_depth = tune(Pokolbin1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Pokolbin1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Pokolbin-1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Pokolbin1_depth) <- new_column_names
SG1Age_on_Pokolbin1_depth[,1] = SG1Age_on_Pokolbin1_depth[,1] * 1000
write.csv(SG1Age_on_Pokolbin1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Pokolbin 1.csv", row.names = FALSE)
