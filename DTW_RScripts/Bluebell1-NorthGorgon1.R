install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Bluebell1 and NorthGorgon1 datasets

Bluebell1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Bluebell1.csv", header=TRUE, stringsAsFactors=FALSE)
Bluebell1=Bluebell1[c(1:11479),] # Oligocene-Miocene
head(Bluebell1)
plot(Bluebell1, type="l", xlim = c(200, 2000), ylim = c(0, 60))

NorthGorgon1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/NorthGorgon1.csv", header=TRUE, stringsAsFactors=FALSE)
NorthGorgon1=NorthGorgon1[c(1:10129),] # Oligocene-Miocene
head(NorthGorgon1)
plot(NorthGorgon1, type="l", xlim = c(250, 1800), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Bluebell1_interpolated <- linterp(Bluebell1, dt = 0.2, genplot = F)
NorthGorgon1_interpolated <- linterp(NorthGorgon1, dt = 0.2, genplot = F)

# Scaling the data
Bmean = Gmean(Bluebell1_interpolated$GR)
Bstd = Gsd(Bluebell1_interpolated$GR)
Bluebell1_scaled = (Bluebell1_interpolated$GR - Bmean)/Bstd
Bluebell1_rescaled = data.frame(Bluebell1_interpolated$DEPT, Bluebell1_scaled)

NGmean = Gmean(NorthGorgon1_interpolated$GR)
NGstd = Gsd(NorthGorgon1_interpolated$GR)
NorthGorgon1_scaled = (NorthGorgon1_interpolated$GR - NGmean)/NGstd
NorthGorgon1_rescaled = data.frame(NorthGorgon1_interpolated$DEPT, NorthGorgon1_scaled)

# Resampling the data using moving window statistics
Bluebell1_scaled = mwStats(Bluebell1_rescaled, cols = 2, win=3, ends = T)
Bluebell1_standardized = data.frame(Bluebell1_scaled$Center_win, Bluebell1_scaled$Average)

NorthGorgon1_scaled = mwStats(NorthGorgon1_rescaled, cols = 2, win=3, ends = T)
NorthGorgon1_standardized = data.frame(NorthGorgon1_scaled$Center_win, NorthGorgon1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Bluebell1_standardized, type="l", xlim = c(200, 2000), ylim = c(-20, 20), xlab = "Bluebell1 Resampled Depth", ylab = "Normalized GR")
plot(NorthGorgon1_standardized, type="l", xlim = c(250, 1800), ylim = c(-20, 20), xlab = "NorthGorgon1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_ng1_b1_ap1 <- dtw(NorthGorgon1_standardized$NorthGorgon1_scaled.Average, Bluebell1_standardized$Bluebell1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = F))
plot(al_ng1_b1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
NorthGorgon1_on_Bluebell1_depth = tune(NorthGorgon1_standardized, cbind(NorthGorgon1_standardized$NorthGorgon1_scaled.Center_win[al_ng1_b1_ap1$index1s], Bluebell1_standardized$Bluebell1_scaled.Center_win[al_ng1_b1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Bluebell1_standardized, type = "l", ylim = c(-20, 20), xlim = c(200, 1900), xlab = "Bluebell1 Resampled Depth", ylab = "Normalized GR")
lines(NorthGorgon1_on_Bluebell1_depth, col = "red")

# DTW Distance
al_ng1_b1_ap1$normalizedDistance
al_ng1_b1_ap1$distance

# Tuning NorthGorgon1 data on Fisher1 depth
NorthGorgon1_on_Fisher1_depth = tune(NorthGorgon1_on_Bluebell1_depth, cbind(Bluebell1_standardized$Bluebell1_scaled.Center_win[al_bl1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_bl1_fi1_ap1$index2s]), extrapolate = F)

# Tuning NorthGorgon1 data on Finucane1 depth
NorthGorgon1_on_Finucane1_depth = tune(NorthGorgon1_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning NorthGorgon1 data on Picard1 depth
NorthGorgon1_on_Picard1_depth = tune(NorthGorgon1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (North Gorgon-1)")
lines(NorthGorgon1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
NorthGorgon1_originalGR_on_Picard1_depth = data.frame(NorthGorgon1_on_Picard1_depth$X1, NorthGorgon1_interpolated$GR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (North Gorgon-1)")
lines(NorthGorgon1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to NorthGorgon1 

U1463Age_on_NorthGorgon1_depth = tune(NorthGorgon1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_NorthGorgon1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "North Gorgon1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_NorthGorgon1_depth) <- new_column_names
write.csv(U1463Age_on_NorthGorgon1_depth, file = "RScripts&Data/Sites Data_Age-NGR/North Gorgon 1.csv", row.names = FALSE)
