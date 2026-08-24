install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Fisher1 and Sultan1 datasets

Fisher1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Fisher1.csv", header=TRUE, stringsAsFactors=FALSE)
Fisher1=Fisher1[c(1:10190),] # Oligocene-Miocene
head(Fisher1)
plot(Fisher1, type="l", xlim = c(100, 1700), ylim = c(0, 50))

Sultan1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Sultan1.csv", header=TRUE, stringsAsFactors=FALSE)
Sultan1=Sultan1[c(1:11867),] # Oligocene-Miocene
head(Sultan1)
plot(Sultan1, type="l", xlim = c(150, 2000), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Fisher1_interpolated <- linterp(Fisher1, dt = 0.2, genplot = F)
Sultan1_interpolated <- linterp(Sultan1, dt = 0.2, genplot = F)

# Scaling the data
Fimean = Gmean(Fisher1_interpolated$GR)
Fistd = Gsd(Fisher1_interpolated$GR)
Fisher1_scaled = (Fisher1_interpolated$GR - Fimean)/Fistd
Fisher1_rescaled = data.frame(Fisher1_interpolated$DEPT, Fisher1_scaled)

Smean = Gmean(Sultan1_interpolated$GR)
Sstd = Gsd(Sultan1_interpolated$GR)
Sultan1_scaled = (Sultan1_interpolated$GR - Smean)/Sstd
Sultan1_rescaled = data.frame(Sultan1_interpolated$DEPT, Sultan1_scaled)

# Resampling the data using moving window statistics
Fisher1_scaled = mwStats(Fisher1_rescaled, cols = 2, win=3, ends = T)
Fisher1_standardized = data.frame(Fisher1_scaled$Center_win, Fisher1_scaled$Average)

Sultan1_scaled = mwStats(Sultan1_rescaled, cols = 2, win=3, ends = T)
Sultan1_standardized = data.frame(Sultan1_scaled$Center_win, Sultan1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Fisher1_standardized, type="l", xlim = c(100, 1700), ylim = c(-20, 20), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
plot(Sultan1_standardized, type="l", xlim = c(150, 2000), ylim = c(-20, 20), xlab = "Sultan1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_s1_fi1_ap1 <- dtw(Sultan1_standardized$Sultan1_scaled.Average, Fisher1_standardized$Fisher1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = F))
plot(al_s1_fi1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Sultan1_on_Fisher1_depth = tune(Sultan1_standardized, cbind(Sultan1_standardized$Sultan1_scaled.Center_win[al_s1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_s1_fi1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Fisher1_standardized, type = "l", ylim = c(-20, 20), xlim = c(100, 1700), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
lines(Sultan1_on_Fisher1_depth, col = "red")

# DTW Distance
al_s1_fi1_ap1$normalizedDistance
al_s1_fi1_ap1$distance

# Tuning Sultan1 data on Finucane1 depth
Sultan1_on_Finucane1_depth = tune(Sultan1_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning Sultan1 data on Picard1 depth
Sultan1_on_Picard1_depth = tune(Sultan1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Sultan-1)")
lines(Sultan1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Sultan1_originalGR_on_Picard1_depth = data.frame(Sultan1_on_Picard1_depth$X1, Sultan1_interpolated[1:9041,2])

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Sultan-1)")
lines(Sultan1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Sultan1 

U1463Age_on_Sultan1_depth = tune(Sultan1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Sultan1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Sultan1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Sultan1_depth) <- new_column_names
write.csv(U1463Age_on_Sultan1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Sultan 1.csv", row.names = FALSE)
