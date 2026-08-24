install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Fisher1 and Brigadier1 datasets

Fisher1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Fisher1.csv", header=TRUE, stringsAsFactors=FALSE)
Fisher1=Fisher1[c(1:10190),] # Oligocene-Miocene
head(Fisher1)
plot(Fisher1, type="l", xlim = c(100, 1700), ylim = c(0, 50))

Brigadier1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Brigadier1.csv", header=TRUE, stringsAsFactors=FALSE)
Brigadier1=Brigadier1[c(1:12660),] # Oligocene-Miocene
head(Brigadier1)
plot(Brigadier1, type="l", xlim = c(350, 2250), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Fisher1_interpolated <- linterp(Fisher1, dt = 0.2, genplot = F)
Brigadier1_interpolated <- linterp(Brigadier1, dt = 0.2, genplot = F)

# Scaling the data
Fimean = Gmean(Fisher1_interpolated$GR)
Fistd = Gsd(Fisher1_interpolated$GR)
Fisher1_scaled = (Fisher1_interpolated$GR - Fimean)/Fistd
Fisher1_rescaled = data.frame(Fisher1_interpolated$DEPT, Fisher1_scaled)

Bmean = Gmean(Brigadier1_interpolated$GR)
Bstd = Gsd(Brigadier1_interpolated$GR)
Brigadier1_scaled = (Brigadier1_interpolated$GR - Bmean)/Bstd
Brigadier1_rescaled = data.frame(Brigadier1_interpolated$DEPT, Brigadier1_scaled)

# Resampling the data using moving window statistics
Fisher1_scaled = mwStats(Fisher1_rescaled, cols = 2, win=3, ends = T)
Fisher1_standardized = data.frame(Fisher1_scaled$Center_win, Fisher1_scaled$Average)

Brigadier1_scaled = mwStats(Brigadier1_rescaled, cols = 2, win=3, ends = T)
Brigadier1_standardized = data.frame(Brigadier1_scaled$Center_win, Brigadier1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Fisher1_standardized, type="l", xlim = c(100, 1700), ylim = c(-20, 20), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
plot(Brigadier1_standardized, type="l", xlim = c(350, 2250), ylim = c(-20, 20), xlab = "Brigadier1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_b1_fi1_ap1 <- dtw(Brigadier1_standardized$Brigadier1_scaled.Average, Fisher1_standardized$Fisher1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_b1_fi1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Brigadier1_on_Fisher1_depth = tune(Brigadier1_standardized, cbind(Brigadier1_standardized$Brigadier1_scaled.Center_win[al_b1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_b1_fi1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Fisher1_standardized, type = "l", ylim = c(-20, 20), xlim = c(100, 1700), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
lines(Brigadier1_on_Fisher1_depth, col = "red")

# DTW Distance
al_b1_fi1_ap1$normalizedDistance
al_b1_fi1_ap1$distance

# Tuning Brigadier1 data on Finucane1 depth
Brigadier1_on_Finucane1_depth = tune(Brigadier1_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning Brigadier1 data on Picard1 depth
Brigadier1_on_Picard1_depth = tune(Brigadier1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Brigadier-1)")
lines(Brigadier1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Brigadier1_originalGR_on_Picard1_depth = data.frame(Brigadier1_on_Picard1_depth$X1, Brigadier1_interpolated[1:9647,2])

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Brigadier-1)")
lines(Brigadier1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Brigadier1 

U1463Age_on_Brigadier1_depth = tune(Brigadier1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Brigadier1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Brigadier1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Brigadier1_depth) <- new_column_names
write.csv(U1463Age_on_Brigadier1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Brigadier 1.csv", row.names = FALSE)


Brigadier1_age_depth <- approx(x = Brigadier1_interpolated$GR,
                               y = Brigadier1_interpolated$DEPT,
                               xout = U1463Age_on_Brigadier1_depth$GR,
                               rule = 1)$y

Brigadier1_agemodel = data.frame(Brigadier1_age_depth, U1463Age_on_Brigadier1_depth$AGE)
Brigadier1_agemodel <- na.omit(Brigadier1_agemodel)
new_column_names1 <- c("Depth", "Age")
colnames(Brigadier1_agemodel) <- new_column_names1
plot(Brigadier1_agemodel)
write.csv(Brigadier1_agemodel, file = "RScripts&Data/Sites Data_Age-NGR/Brigadier1_DepthAge.csv", row.names = FALSE)
