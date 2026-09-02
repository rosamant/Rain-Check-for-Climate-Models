install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Gorgon1 and Ramillies1 datasets

Gorgon1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Gorgon1.csv", header=TRUE, stringsAsFactors=FALSE)
Gorgon1=Gorgon1[c(1:9280),] # Data required till Eocene-Miocene Unconformity
head(Gorgon1)
plot(Gorgon1, type="l", xlim = c(300, 1700), ylim = c(0, 50))

Ramillies1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Ramillies1.csv", header=TRUE, stringsAsFactors=FALSE)
Ramillies1 = Ramillies1[c(1:6862),] # Data required till Eocene-Miocene Unconformity
head(Ramillies1)
plot(Ramillies1, type="l", xlim = c(150, 1550), ylim = c(0, 40))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Gorgon1_interpolated <- linterp(Gorgon1, dt = 0.2, genplot = F)
Ramillies1_interpolated <- linterp(Ramillies1, dt = 0.2, genplot = F)

# Scaling the data
Gmean = Gmean(Gorgon1_interpolated$GR)
Gstd = Gsd(Gorgon1_interpolated$GR)
Gorgon1_scaled = (Gorgon1_interpolated$GR - Gmean)/Gstd
Gorgon1_rescaled = data.frame(Gorgon1_interpolated$DEPT, Gorgon1_scaled)

Rmean = Gmean(Ramillies1_interpolated$GR)
Rstd = Gsd(Ramillies1_interpolated$GR)
Ramillies1_scaled = (Ramillies1_interpolated$GR - Rmean)/Rstd
Ramillies1_rescaled = data.frame(Ramillies1_interpolated$DEPT, Ramillies1_scaled)

# Resampling the data using moving window statistics
Gorgon1_scaled = mwStats(Gorgon1_rescaled, cols = 2, win = 3, ends = T)
Gorgon1_standardized = data.frame(Gorgon1_scaled$Center_win, Gorgon1_scaled$Average)

Ramillies1_scaled = mwStats(Ramillies1_rescaled, cols = 2, win = 3, ends = T)
Ramillies1_standardized = data.frame(Ramillies1_scaled$Center_win, Ramillies1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Gorgon1_standardized, type="l", xlim = c(300, 1700), ylim = c(-20, 20), xlab = "Gorgon1 Resampled Depth", ylab = "Normalized GR")
plot(Ramillies1_standardized, type="l", xlim = c(150, 1550), ylim = c(-20, 20), xlab = "Ramillies1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_r1_g1_ap1 <- dtw(Ramillies1_standardized$Ramillies1_scaled.Average, Gorgon1_standardized$Gorgon1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = F, open.end = T))
plot(al_r1_g1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Ramillies1_on_Gorgon1_depth = tune(Ramillies1_standardized, cbind(Ramillies1_standardized$Ramillies1_scaled.Center_win[al_r1_g1_ap1$index1s], Gorgon1_standardized$Gorgon1_scaled.Center_win[al_r1_g1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Gorgon1_standardized, type = "l", ylim = c(-20, 20), xlim = c(300, 1700), xlab = "Gorgon1 Resampled Depth", ylab = "Normalized GR")
lines(Ramillies1_on_Gorgon1_depth, col = "red")

# DTW Distance measure
al_r1_g1_ap1$normalizedDistance
al_r1_g1_ap1$distance

# Tuning Ramillies1 data on Bluebell1 depth
Ramillies1_on_Bluebell1_depth = tune(Ramillies1_on_Gorgon1_depth, cbind(Gorgon1_standardized$Gorgon1_scaled.Center_win[al_g1_b1_ap1$index1s], Bluebell1_standardized$Bluebell1_scaled.Center_win[al_g1_b1_ap1$index2s]), extrapolate = F)

# Tuning Ramillies1 data on Fisher1 depth
Ramillies1_on_Fisher1_depth = tune(Ramillies1_on_Bluebell1_depth, cbind(Bluebell1_standardized$Bluebell1_scaled.Center_win[al_bl1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_bl1_fi1_ap1$index2s]), extrapolate = F)

# Tuning Ramillies1 data on Finucane1 depth
Ramillies1_on_Finucane1_depth = tune(Ramillies1_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning Ramillies1 data on Picard1 depth
Ramillies1_on_Picard1_depth = tune(Ramillies1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Bowers-1)")
lines(Ramillies1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Ramillies1_originalGR_on_Picard1_depth = data.frame(Ramillies1_on_Picard1_depth$X1, Ramillies1_interpolated[1:6838,2])

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Bowers-1)")
lines(Ramillies1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Ramillies1 

U1463Age_on_Ramillies1_depth = tune(Ramillies1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Ramillies1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Ramillies1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Ramillies1_depth) <- new_column_names
write.csv(U1463Age_on_Ramillies1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Ramillies 1.csv", row.names = FALSE)
