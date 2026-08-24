install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Fisher1 and Wilcox2 datasets

Fisher1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Fisher1.csv", header=TRUE, stringsAsFactors=FALSE)
Fisher1=Fisher1[c(1:10190),] # Oligocene-Miocene
head(Fisher1)
plot(Fisher1, type="l", xlim = c(100, 1700), ylim = c(0, 50))

Wilcox2 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Wilcox2.csv", header=TRUE, stringsAsFactors=FALSE)
Wilcox2=Wilcox2[c(1:10800),] # Oligocene-Miocene
head(Wilcox2)
plot(Wilcox2, type="l", xlim = c(100, 1800), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Fisher1_interpolated <- linterp(Fisher1, dt = 0.2, genplot = F)
Wilcox2_interpolated <- linterp(Wilcox2, dt = 0.2, genplot = F)

# Scaling the data
Fimean = Gmean(Fisher1_interpolated$GR)
Fistd = Gsd(Fisher1_interpolated$GR)
Fisher1_scaled = (Fisher1_interpolated$GR - Fimean)/Fistd
Fisher1_rescaled = data.frame(Fisher1_interpolated$DEPT, Fisher1_scaled)

Wmean = Gmean(Wilcox2_interpolated$GR)
Wstd = Gsd(Wilcox2_interpolated$GR)
Wilcox2_scaled = (Wilcox2_interpolated$GR - Wmean)/Wstd
Wilcox2_rescaled = data.frame(Wilcox2_interpolated$DEPT, Wilcox2_scaled)

# Resampling the data using moving window statistics
Fisher1_scaled = mwStats(Fisher1_rescaled, cols = 2, win=3, ends = T)
Fisher1_standardized = data.frame(Fisher1_scaled$Center_win, Fisher1_scaled$Average)

Wilcox2_scaled = mwStats(Wilcox2_rescaled, cols = 2, win=3, ends = T)
Wilcox2_standardized = data.frame(Wilcox2_scaled$Center_win, Wilcox2_scaled$Average)

# Plotting the rescaled and resampled data
plot(Fisher1_standardized, type="l", xlim = c(100, 1700), ylim = c(-20, 20), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
plot(Wilcox2_standardized, type="l", xlim = c(100, 1800), ylim = c(-20, 20), xlab = "Wilcox2 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_w1_fi1_ap1 <- dtw(Wilcox2_standardized$Wilcox2_scaled.Average, Fisher1_standardized$Fisher1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_w1_fi1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Wilcox2_on_Fisher1_depth = tune(Wilcox2_standardized, cbind(Wilcox2_standardized$Wilcox2_scaled.Center_win[al_w1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_w1_fi1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Fisher1_standardized, type = "l", ylim = c(-20, 20), xlim = c(100, 1700), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
lines(Wilcox2_on_Fisher1_depth, col = "red")

# DTW Distance
al_w1_fi1_ap1$normalizedDistance
al_w1_fi1_ap1$distance

# Tuning Wilcox2 data on Finucane1 depth
Wilcox2_on_Finucane1_depth = tune(Wilcox2_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning Wilcox2 data on Picard1 depth
Wilcox2_on_Picard1_depth = tune(Wilcox2_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Wilcox-2)")
lines(Wilcox2_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Wilcox2_originalGR_on_Picard1_depth = data.frame(Wilcox2_on_Picard1_depth$X1, Wilcox2_interpolated$GR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Wilcox-2)")
lines(Wilcox2_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Wilcox2 

U1463Age_on_Wilcox2_depth = tune(Wilcox2_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Wilcox2_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Wilcox2")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Wilcox2_depth) <- new_column_names
write.csv(U1463Age_on_Wilcox2_depth, file = "RScripts&Data/Sites Data_Age-NGR/Wilcox 2.csv", row.names = FALSE)
