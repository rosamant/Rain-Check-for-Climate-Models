install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Fisher1 and Goodwyn6 datasets

Fisher1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Fisher1.csv", header=TRUE, stringsAsFactors=FALSE)
Fisher1=Fisher1[c(1:10190),] # Oligocene-Miocene
head(Fisher1)
plot(Fisher1, type="l", xlim = c(100, 1700), ylim = c(0, 50))

Goodwyn6 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Goodwyn6.csv", header=TRUE, stringsAsFactors=FALSE)
Goodwyn6=Goodwyn6[c(1:10620),] # Oligocene-Miocen
head(Goodwyn6)
plot(Goodwyn6, type="l", xlim = c(150, 1800), ylim = c(0, 50))

# Linear interpolation of datasets
Fisher1_interpolated <- linterp(Fisher1, dt = 0.2, genplot = F)
Goodwyn6_interpolated <- linterp(Goodwyn6, dt = 0.2, genplot = F)

# Scaling the data
Fimean = Gmean(Fisher1_interpolated$GR)
Fistd = Gsd(Fisher1_interpolated$GR)
Fisher1_scaled = (Fisher1_interpolated$GR - Fimean)/Fistd
Fisher1_rescaled = data.frame(Fisher1_interpolated$DEPT, Fisher1_scaled)

Gomean = Gmean(Goodwyn6_interpolated$GR)
Gostd = Gsd(Goodwyn6_interpolated$GR)
Goodwyn6_scaled = (Goodwyn6_interpolated$GR - Gomean)/Gostd
Goodwyn6_rescaled = data.frame(Goodwyn6_interpolated$DEPT, Goodwyn6_scaled)

# Resampling the data using moving window statistics
Fisher1_scaled = mwStats(Fisher1_rescaled, cols = 2, win=3, ends = T)
Fisher1_standardized = data.frame(Fisher1_scaled$Center_win, Fisher1_scaled$Average)

Goodwyn6_scaled = mwStats(Goodwyn6_rescaled, cols = 2, win=3, ends = T)
Goodwyn6_standardized = data.frame(Goodwyn6_scaled$Center_win, Goodwyn6_scaled$Average)

# Plotting the rescaled and resampled data
plot(Fisher1_standardized, type="l", xlim = c(100, 1700), ylim = c(-20, 20), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
plot(Goodwyn6_standardized, type="l", xlim = c(150, 1800), ylim = c(-20, 20), xlab = "Goodwyn6 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_g6_fi1_ap1 <- dtw(Goodwyn6_standardized$Goodwyn6_scaled.Average, Fisher1_standardized$Fisher1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_g6_fi1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Goodwyn6_on_Fisher1_depth = tune(Goodwyn6_standardized, cbind(Goodwyn6_standardized$Goodwyn6_scaled.Center_win[al_g6_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_g6_fi1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Fisher1_standardized, type = "l", ylim = c(-20, 20), xlim = c(100, 1700), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
lines(Goodwyn6_on_Fisher1_depth, col = "red")

# DTW Distance and RMSE

al_g6_fi1_ap1$normalizedDistance
al_g6_fi1_ap1$distance

# Tuning Goodwyn6 data on Finucane1 depth
Goodwyn6_on_Finucane1_depth = tune(Goodwyn6_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning Goodwyn6 data on Picard1 depth
Goodwyn6_on_Picard1_depth = tune(Goodwyn6_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Goodwyn-6)")
lines(Goodwyn6_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Goodwyn6_originalGR_on_Picard1_depth = data.frame(Goodwyn6_on_Picard1_depth$X1, Goodwyn6_interpolated$GR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Goodwyn-6)")
lines(Goodwyn6_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Goodwyn6 

U1463Age_on_Goodwyn6_depth = tune(Goodwyn6_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Goodwyn6_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Goodwyn6")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Goodwyn6_depth) <- new_column_names
write.csv(U1463Age_on_Goodwyn6_depth, file = "RScripts&Data/Sites Data_Age-NGR/Goodwyn 6.csv", row.names = FALSE)
