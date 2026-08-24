install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Plymouth1 and Athena1 datasets

Plymouth1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Plymouth 1.csv", header=TRUE, stringsAsFactors=FALSE)
Plymouth1=Plymouth1[c(1:7444),] # Oligocene-Miocene
head(Plymouth1)
plot(Plymouth1, type="l", xlim = c(200, 1650), ylim = c(0, 70))

Athena1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Athena 1.csv", header=TRUE, stringsAsFactors=FALSE)
Athena1=Athena1[c(1:8190),] # Oligocene-Miocene
head(Athena1)
plot(Athena1, type="l", xlim = c(150, 1750), ylim = c(0, 70))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Plymouth1_interpolated <- linterp(Plymouth1, dt = 0.2, genplot = F)
Athena1_interpolated <- linterp(Athena1, dt = 0.2, genplot = F)

# Scaling the data
Plmean = Gmean(Plymouth1_interpolated$GR)
Plstd = Gsd(Plymouth1_interpolated$GR)
Plymouth1_scaled = (Plymouth1_interpolated$GR - Plmean)/Plstd
Plymouth1_rescaled = data.frame(Plymouth1_interpolated$DEPT, Plymouth1_scaled)

Atmean = Gmean(Athena1_interpolated$GR)
Atstd = Gsd(Athena1_interpolated$GR)
Athena1_scaled = (Athena1_interpolated$GR - Atmean)/Atstd
Athena1_rescaled = data.frame(Athena1_interpolated$DEPT, Athena1_scaled)

# Resampling the data using moving window statistics
Plymouth1_scaled = mwStats(Plymouth1_rescaled, cols = 2, win=3, ends = T)
Plymouth1_standardized = data.frame(Plymouth1_scaled$Center_win, Plymouth1_scaled$Average)

Athena1_scaled = mwStats(Athena1_rescaled, cols = 2, win=3, ends = T)
Athena1_standardized = data.frame(Athena1_scaled$Center_win, Athena1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Plymouth1_standardized, type="l", xlim = c(200, 1650), ylim = c(-20, 20), xlab = "Plymouth1 Resampled Depth", ylab = "Normalized GR")
plot(Athena1_standardized, type="l", xlim = c(150, 1750), ylim = c(-20, 20), xlab = "Athena1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_at1_pl1_ap1 <- dtw(Athena1_standardized$Athena1_scaled.Average, Plymouth1_standardized$Plymouth1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = F, open.end = T))
plot(al_at1_pl1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Athena1_on_Plymouth1_depth = tune(Athena1_standardized, cbind(Athena1_standardized$Athena1_scaled.Center_win[al_at1_pl1_ap1$index1s], Plymouth1_standardized$Plymouth1_scaled.Center_win[al_at1_pl1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Plymouth1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1700), xlab = "Plymouth1 Resampled Depth", ylab = "Normalized GR")
lines(Athena1_on_Plymouth1_depth, col = "red")

# DTW Distance 
al_at1_pl1_ap1$normalizedDistance
al_at1_pl1_ap1$distance

# Tuning Athena1 data on Finucane1 depth
Athena1_on_Finucane1_depth = tune(Athena1_on_Plymouth1_depth, cbind(Plymouth1_standardized$Plymouth1_scaled.Center_win[al_pl1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_pl1_f1_ap1$index2s]), extrapolate = F)

# Tuning Athena1 data on Picard1 depth
Athena1_on_Picard1_depth = tune(Athena1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Athena1_originalGR_on_Picard1_depth = data.frame(Athena1_on_Picard1_depth$X1, Athena1_interpolated$GR)

dev.off()
plot(Picard1_originalGR, type = "l", ylim = c(0, 60), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Athena-1)")
lines(Athena1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Athena1 

U1463Age_on_Athena1_depth = tune(Athena1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Athena1_depth, type = "l", ylim = c(0, 70), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Athena1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Athena1_depth) <- new_column_names
write.csv(U1463Age_on_Athena1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Athena 1.csv", row.names = FALSE)
