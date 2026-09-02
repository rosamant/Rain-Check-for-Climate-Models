install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Fisher1 and NorthRankin6 datasets

Fisher1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Fisher1.csv", header=TRUE, stringsAsFactors=FALSE)
Fisher1=Fisher1[c(1:10190),] # Oligocene-Miocene
head(Fisher1)
plot(Fisher1, type="l", xlim = c(100, 1700), ylim = c(0, 50))

NorthRankin6 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/NorthRankin6.csv", header=TRUE, stringsAsFactors=FALSE)
NorthRankin6=NorthRankin6[c(1:11269),] # Oligocene-Miocene
head(NorthRankin6)
plot(NorthRankin6, type="l", xlim = c(150, 1900), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Fisher1_interpolated <- linterp(Fisher1, dt = 0.2, genplot = F)
NorthRankin6_interpolated <- linterp(NorthRankin6, dt = 0.2, genplot = F)

# Scaling the data
Fmean = Gmean(Fisher1_interpolated$GR)
Fstd = Gsd(Fisher1_interpolated$GR)
Fisher1_scaled = (Fisher1_interpolated$GR - Fmean)/Fstd
Fisher1_rescaled = data.frame(Fisher1_interpolated$DEPT, Fisher1_scaled)

NRmean = Gmean(NorthRankin6_interpolated$GR)
NRstd = Gsd(NorthRankin6_interpolated$GR)
NorthRankin6_scaled = (NorthRankin6_interpolated$GR - NRmean)/NRstd
NorthRankin6_rescaled = data.frame(NorthRankin6_interpolated$DEPT, NorthRankin6_scaled)

# Resampling the data using moving window statistics
Fisher1_scaled = mwStats(Fisher1_rescaled, cols = 2, win=3, end = T)
Fisher1_standardized = data.frame(Fisher1_scaled$Center_win, Fisher1_scaled$Average)

NorthRankin6_scaled = mwStats(NorthRankin6_rescaled, cols = 2, win=3, end = T)
NorthRankin6_standardized = data.frame(NorthRankin6_scaled$Center_win, NorthRankin6_scaled$Average)

# Plotting the rescaled and resampled data
plot(Fisher1_standardized, type="l", xlim = c(100, 1700), ylim = c(-20, 20), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
plot(NorthRankin6_standardized, type="l", xlim = c(150, 1900), ylim = c(-20, 20), xlab = "NorthRankin6 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_nr6_fi1_ap1 <- dtw(NorthRankin6_standardized$NorthRankin6_scaled.Average, Fisher1_standardized$Fisher1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_nr6_fi1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
NorthRankin6_on_Fisher1_depth = tune(NorthRankin6_standardized, cbind(NorthRankin6_standardized$NorthRankin6_scaled.Center_win[al_nr6_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_nr6_fi1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Fisher1_standardized, type = "l", ylim = c(-20, 20), xlim = c(100, 1700), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
lines(NorthRankin6_on_Fisher1_depth, col = "red")

# DTW Distance and RMSE

al_nr6_fi1_ap1$normalizedDistance
al_nr6_fi1_ap1$distance

# Tuning NorthRankin6 data on Finucane1 depth
NorthRankin6_on_Finucane1_depth = tune(NorthRankin6_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning NorthRankin6 data on Picard1 depth
NorthRankin6_on_Picard1_depth = tune(NorthRankin6_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "NorthRankin6 Resampled Depth", ylab = "Normalized GR (NorthRankin-1)")
lines(NorthRankin6_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
NorthRankin6_originalGR_on_Picard1_depth = data.frame(NorthRankin6_on_Picard1_depth$X1, NorthRankin6_interpolated[1:8586,2])

plot(Picard1_originalGR, type = "l", ylim = c(0, 60), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (NorthRankin-6)")
lines(NorthRankin6_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to NorthRankin6 

U1463Age_on_NorthRankin6_depth = tune(NorthRankin6_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_NorthRankin6_depth, type = "l", ylim = c(0, 60), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "NorthRankin6")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_NorthRankin6_depth) <- new_column_names
write.csv(U1463Age_on_NorthRankin6_depth, file = "RScripts&Data/Sites Data_Age-NGR/NorthRankin 6.csv", row.names = FALSE)
