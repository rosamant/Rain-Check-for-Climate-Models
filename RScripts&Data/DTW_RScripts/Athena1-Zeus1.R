install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Athena1 and Zeus1 datasets

Athena1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Athena 1.csv", header=TRUE, stringsAsFactors=FALSE)
Athena1=Athena1[c(1:8190),] # Oligocene-Miocene
head(Athena1)
plot(Athena1, type="l", xlim = c(150, 1750), ylim = c(0, 70))

Zeus1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Zeus1.csv", header=TRUE, stringsAsFactors=FALSE)
Zeus1=Zeus1[c(201:9889),] # Oligocene-Miocene
head(Zeus1)
plot(Zeus1, type="l", xlim = c(250, 2150), ylim = c(0, 70))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Athena1_interpolated <- linterp(Athena1, dt = 0.2, genplot = F)
Zeus1_interpolated <- linterp(Zeus1, dt = 0.2, genplot = F)

# Scaling the data
Atmean = Gmean(Athena1_interpolated$GR)
Atstd = Gsd(Athena1_interpolated$GR)
Athena1_scaled = (Athena1_interpolated$GR - Atmean)/Atstd
Athena1_rescaled = data.frame(Athena1_interpolated$DEPT, Athena1_scaled)

Zmean = Gmean(Zeus1_interpolated$GR)
Zstd = Gsd(Zeus1_interpolated$GR)
Zeus1_scaled = (Zeus1_interpolated$GR - Zmean)/Zstd
Zeus1_rescaled = data.frame(Zeus1_interpolated$DEPT, Zeus1_scaled)

# Resampling the data using moving window statistics
Athena1_scaled = mwStats(Athena1_rescaled, cols = 2, win=3, ends = T)
Athena1_standardized = data.frame(Athena1_scaled$Center_win, Athena1_scaled$Average)

Zeus1_scaled = mwStats(Zeus1_rescaled, cols = 2, win=3, ends = T)
Zeus1_standardized = data.frame(Zeus1_scaled$Center_win, Zeus1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Athena1_standardized, type="l", xlim = c(150, 1750), ylim = c(-20, 20), xlab = "Athena1 Resampled Depth", ylab = "Normalized GR")
plot(Zeus1_standardized, type="l", xlim = c(250, 2150), ylim = c(-20, 20), xlab = "Zeus1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_z1_at1_ap1 <- dtw(Zeus1_standardized$Zeus1_scaled.Average, Athena1_standardized$Athena1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_z1_at1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Zeus1_on_Athena1_depth = tune(Zeus1_standardized, cbind(Zeus1_standardized$Zeus1_scaled.Center_win[al_z1_at1_ap1$index1s], Athena1_standardized$Athena1_scaled.Center_win[al_z1_at1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Athena1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1800), xlab = "Angel2 Resampled Depth", ylab = "Normalized GR")
lines(Zeus1_on_Athena1_depth, col = "red")

# DTW Distance
al_z1_at1_ap1$normalizedDistance
al_z1_at1_ap1$distance

# Tuning Zeus1 data on Plymouth1 depth
Zeus1_on_Plymouth1_depth = tune(Zeus1_on_Athena1_depth, cbind(Athena1_standardized$Athena1_scaled.Center_win[al_at1_pl1_ap1$index1s], Plymouth1_standardized$Plymouth1_scaled.Center_win[al_at1_pl1_ap1$index2s]), extrapolate = F)

# Tuning Zeus1 data on Finucane1 depth
Zeus1_on_Finucane1_depth = tune(Zeus1_on_Plymouth1_depth, cbind(Plymouth1_standardized$Plymouth1_scaled.Center_win[al_pl1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_pl1_f1_ap1$index2s]), extrapolate = F)

# Tuning Athena1 data on Picard1 depth
Zeus1_on_Picard1_depth = tune(Zeus1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Zeus-1)")
lines(Zeus1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Zeus1_originalGR_on_Picard1_depth = data.frame(Zeus1_on_Picard1_depth$X1, Zeus1_interpolated$GR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 70), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Zeus-1)")
lines(Zeus1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Zeus1 

U1463Age_on_Zeus1_depth = tune(Zeus1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Zeus1_depth, type = "l", ylim = c(0, 70), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Zeus1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Zeus1_depth) <- new_column_names
write.csv(U1463Age_on_Zeus1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Zeus 1.csv", row.names = FALSE)
