install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Finucane1 and Sable1 datasets

Finucane1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Finucane1.csv", header=TRUE, stringsAsFactors=FALSE)
Finucane1=Finucane1[c(1:9080),] # Oligocene-Miocene
head(Finucane1)
plot(Finucane1, type="l", xlim = c(150, 1500), ylim = c(0, 50))

Sable1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Sable1.csv", header=TRUE, stringsAsFactors=FALSE)
Sable1=Sable1[c(1:9102),] # Oligocene-Miocene
head(Sable1)
plot(Sable1, type="l", xlim = c(150, 1550), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Finucane1_interpolated <- linterp(Finucane1, dt = 0.2, genplot = F)
Sable1_interpolated <- linterp(Sable1, dt = 0.2, genplot = F)

# Scaling the data
Fmean = Gmean(Finucane1_interpolated$GR)
Fstd = Gsd(Finucane1_interpolated$GR)
Finucane1_scaled = (Finucane1_interpolated$GR - Fmean)/Fstd
Finucane1_rescaled = data.frame(Finucane1_interpolated$DEPT, Finucane1_scaled)

Smean = Gmean(Sable1_interpolated$GR)
Sstd = Gsd(Sable1_interpolated$GR)
Sable1_scaled = (Sable1_interpolated$GR - Smean)/Sstd
Sable1_rescaled = data.frame(Sable1_interpolated$DEPT, Sable1_scaled)

# Resampling the data using moving window statistics
Finucane1_scaled = mwStats(Finucane1_rescaled, cols = 2, win=3, ends = T)
Finucane1_standardized = data.frame(Finucane1_scaled$Center_win, Finucane1_scaled$Average)

Sable1_scaled = mwStats(Sable1_rescaled, cols = 2, win=3, ends = T)
Sable1_standardized = data.frame(Sable1_scaled$Center_win, Sable1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Finucane1_standardized, type="l", xlim = c(150, 1700), ylim = c(-20, 20), xlab = "Finucane1 Resampled Depth", ylab = "Normalized GR")
plot(Sable1_standardized, type="l", xlim = c(150, 1700), ylim = c(-20, 20), xlab = "Sable1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_s1_f1_ap1 <- dtw(Sable1_standardized$Sable1_scaled.Average, Finucane1_standardized$Finucane1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_s1_f1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Sable1_on_Finucane1_depth = tune(Sable1_standardized, cbind(Sable1_standardized$Sable1_scaled.Center_win[al_s1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_s1_f1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Finucane1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1500), xlab = "Finucane1 Resampled Depth", ylab = "Normalized GR")
lines(Sable1_on_Finucane1_depth, col = "red")

# DTW Distance 
al_s1_f1_ap1$normalizedDistance
al_s1_f1_ap1$distance

# Tuning Sable1 data on Picard1 depth
Sable1_on_Picard1_depth = tune(Sable1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Sable1_originalGR_on_Picard1_depth = data.frame(Sable1_on_Picard1_depth$X1, Sable1_interpolated[1:6934,2])

dev.off()
plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Sable-1)")
lines(Sable1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Sable1 

U1463Age_on_Sable1_depth = tune(Sable1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Sable1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Sable1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Sable1_depth) <- new_column_names
write.csv(U1463Age_on_Sable1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Sable 1.csv", row.names = FALSE)
