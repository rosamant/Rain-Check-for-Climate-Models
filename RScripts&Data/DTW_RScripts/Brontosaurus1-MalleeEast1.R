install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Brontosaurus 1 and MalleeEast 1 datasets

Brontosaurus1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Brontosaurus 1.csv", header=TRUE, stringsAsFactors=FALSE)
Brontosaurus1=Brontosaurus1[c(1:4527),] # Oligocene-Miocene
head(Brontosaurus1)
plot(Brontosaurus1, type="l", xlim = c(150, 1050), ylim = c(20, 60))

MalleeEast1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Mallee East 1.csv", header=TRUE, stringsAsFactors=FALSE)
MalleeEast1=MalleeEast1[c(40:4394),] # Oligocene-Miocene
head(MalleeEast1)
plot(MalleeEast1, type="l", xlim = c(150, 1050), ylim = c(0, 60))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Brontosaurus1_interpolated <- linterp(Brontosaurus1, dt = 0.2, genplot = F)
MalleeEast1_interpolated <- linterp(MalleeEast1, dt = 0.2, genplot = F)

# Scaling the data
Bmean = Gmean(Brontosaurus1_interpolated$GR)
Bstd = Gsd(Brontosaurus1_interpolated$GR)
Brontosaurus1_scaled = (Brontosaurus1_interpolated$GR - Bmean)/Bstd
Brontosaurus1_rescaled = data.frame(Brontosaurus1_interpolated$DEPT, Brontosaurus1_scaled)

Memean = Gmean(MalleeEast1_interpolated$GR)
Mestd = Gsd(MalleeEast1_interpolated$GR)
MalleeEast1_scaled = (MalleeEast1_interpolated$GR - Memean)/Mestd
MalleeEast1_rescaled = data.frame(MalleeEast1_interpolated$DEPT, MalleeEast1_scaled)

# Resampling the data using moving window statistics
Brontosaurus1_scaled = mwStats(Brontosaurus1_rescaled, cols = 2, win=3, ends = T)
Brontosaurus1_standardized = data.frame(Brontosaurus1_scaled$Center_win, Brontosaurus1_scaled$Average)

MalleeEast1_scaled = mwStats(MalleeEast1_rescaled, cols = 2, win=3, ends = T)
MalleeEast1_standardized = data.frame(MalleeEast1_scaled$Center_win, MalleeEast1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Brontosaurus1_standardized, type="l", xlim = c(150, 1050), ylim = c(-20, 20), xlab = "Brontosaurus1 Resampled Depth", ylab = "Normalized GR")
plot(MalleeEast1_standardized, type="l", xlim = c(150, 1050), ylim = c(-10, 20), xlab = "MalleeEast 1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_me1_b1_ap1 <- dtw(MalleeEast1_standardized$MalleeEast1_scaled.Average, Brontosaurus1_standardized$Brontosaurus1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_me1_b1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
MalleeEast1_on_Brontosaurus1_depth = tune(MalleeEast1_standardized, cbind(MalleeEast1_standardized$MalleeEast1_scaled.Center_win[al_me1_b1_ap1$index1s], Brontosaurus1_standardized$Brontosaurus1_scaled.Center_win[al_me1_b1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Brontosaurus1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1050), xlab = "Brontosaurus1 Resampled Depth", ylab = "Normalized GR")
lines(MalleeEast1_on_Brontosaurus1_depth, col = "red")

# DTW Distance
al_me1_b1_ap1$normalizedDistance
al_me1_b1_ap1$distance

# Tuning MalleeEast1 data on SG1 depth
MalleeEast1_on_SG1_depth = tune(MalleeEast1_on_Brontosaurus1_depth, cbind(Brontosaurus11_standardized$Brontosaurus11_scaled.Center_win[al_sg1_b1_ap2$index2s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_sg1_b1_ap2$index1s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (MalleeEast-1)")
lines(MalleeEast1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
MalleeEast1_originalGR_on_SG1_depth = data.frame(MalleeEast1_on_SG1_depth$X1, MalleeEast1_interpolated[1284:4367,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (MalleeEast-1)")
lines(MalleeEast1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to MalleeEast1 

SG1Age_on_MalleeEast1_depth = tune(MalleeEast1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_MalleeEast1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "MalleeEast-1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_MalleeEast1_depth) <- new_column_names
SG1Age_on_MalleeEast1_depth[,1] = SG1Age_on_MalleeEast1_depth[,1] * 1000
write.csv(SG1Age_on_MalleeEast1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Mallee East 1.csv", row.names = FALSE)
