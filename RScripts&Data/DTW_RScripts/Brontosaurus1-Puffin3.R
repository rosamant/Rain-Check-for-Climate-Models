install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Brontosaurus 1 and Puffin 3 datasets

Brontosaurus1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Brontosaurus 1.csv", header=TRUE, stringsAsFactors=FALSE)
Brontosaurus1=Brontosaurus1[c(1:4527),] # Oligocene-Miocene
head(Brontosaurus1)
plot(Brontosaurus1, type="l", xlim = c(150, 1050), ylim = c(20, 60))

Puffin3 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Puffin 3.csv", header=TRUE, stringsAsFactors=FALSE)
Puffin3=Puffin3[c(1:5781),] # Oligocene-Miocene
head(Puffin3)
plot(Puffin3, type="l", xlim = c(130, 1050), ylim = c(0, 60))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Brontosaurus1_interpolated <- linterp(Brontosaurus1, dt = 0.2, genplot = F)
Puffin3_interpolated <- linterp(Puffin3, dt = 0.2, genplot = F)

# Scaling the data
Bmean = Gmean(Brontosaurus1_interpolated$GR)
Bstd = Gsd(Brontosaurus1_interpolated$GR)
Brontosaurus1_scaled = (Brontosaurus1_interpolated$GR - Bmean)/Bstd
Brontosaurus1_rescaled = data.frame(Brontosaurus1_interpolated$DEPT, Brontosaurus1_scaled)

Pu3mean = Gmean(Puffin3_interpolated$GR)
Pu3std = Gsd(Puffin3_interpolated$GR)
Puffin3_scaled = (Puffin3_interpolated$GR - Pu3mean)/Pu3std
Puffin3_rescaled = data.frame(Puffin3_interpolated$DEPT, Puffin3_scaled)

# Resampling the data using moving window statistics
Brontosaurus1_scaled = mwStats(Brontosaurus1_rescaled, cols = 2, win=3, ends = T)
Brontosaurus1_standardized = data.frame(Brontosaurus1_scaled$Center_win, Brontosaurus1_scaled$Average)

Puffin3_scaled = mwStats(Puffin3_rescaled, cols = 2, win=3, ends = T)
Puffin3_standardized = data.frame(Puffin3_scaled$Center_win, Puffin3_scaled$Average)

# Plotting the rescaled and resampled data
plot(Brontosaurus1_standardized, type="l", xlim = c(150, 1050), ylim = c(-20, 20), xlab = "Brontosaurus1 Resampled Depth", ylab = "Normalized GR")
plot(Puffin3_standardized, type="l", xlim = c(130, 1050), ylim = c(-10, 20), xlab = "Puffin 2 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_p3_b1_ap1 <- dtw(Puffin3_standardized$Puffin3_scaled.Average, Brontosaurus1_standardized$Brontosaurus1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = F, open.end = T))
plot(al_p3_b1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Puffin3_on_Brontosaurus1_depth = tune(Puffin3_standardized, cbind(Puffin3_standardized$Puffin3_scaled.Center_win[al_p3_b1_ap1$index1s], Brontosaurus1_standardized$Brontosaurus1_scaled.Center_win[al_p3_b1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Brontosaurus1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1050), xlab = "Brontosaurus1 Resampled Depth", ylab = "Normalized GR")
lines(Puffin3_on_Brontosaurus1_depth, col = "red")

# DTW Distance
al_p3_b1_ap1$normalizedDistance
al_p3_b1_ap1$distance

# Tuning Puffin3 data on SG1 depth
Puffin3_on_SG1_depth = tune(Puffin3_on_Brontosaurus1_depth, cbind(Brontosaurus11_standardized$Brontosaurus11_scaled.Center_win[al_sg1_b1_ap2$index2s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_sg1_b1_ap2$index1s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (Puffin-3)")
lines(Puffin3_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Puffin3_originalGR_on_SG1_depth = data.frame(Puffin3_on_SG1_depth$X1, Puffin3_interpolated[907:4449,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (Puffin-3)")
lines(Puffin3_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Puffin3 

SG1Age_on_Puffin3_depth = tune(Puffin3_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Puffin3_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Puffin3")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Puffin3_depth) <- new_column_names
SG1Age_on_Puffin3_depth[,1] = SG1Age_on_Puffin3_depth[,1] * 1000
write.csv(SG1Age_on_Puffin3_depth, file = "RScripts&Data/Sites Data_Age-NGR/Puffin 3.csv", row.names = FALSE)
