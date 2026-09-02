install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Fisher1 and Bluebell1 datasets

Fisher1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Fisher1.csv", header=TRUE, stringsAsFactors=FALSE)
Fisher1=Fisher1[c(1:10190),] # Oligocene-Miocene
head(Fisher1)
plot(Fisher1, type="l", xlim = c(100, 1700), ylim = c(0, 50))

Bluebell1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Bluebell1.csv", header=TRUE, stringsAsFactors=FALSE)
Bluebell1=Bluebell1[c(1:11479),] # Oligocene-Miocene
head(Bluebell1)
plot(Bluebell1, type="l", xlim = c(200, 2000), ylim = c(0, 60))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Fisher1_interpolated <- linterp(Fisher1, dt = 0.2, genplot = F)
Bluebell1_interpolated <- linterp(Bluebell1, dt = 0.2, genplot = F)

# Scaling the data
Fimean = Gmean(Fisher1_interpolated$GR)
Fistd = Gsd(Fisher1_interpolated$GR)
Fisher1_scaled = (Fisher1_interpolated$GR - Fimean)/Fistd
Fisher1_rescaled = data.frame(Fisher1_interpolated$DEPT, Fisher1_scaled)

Blmean = Gmean(Bluebell1_interpolated$GR)
Blstd = Gsd(Bluebell1_interpolated$GR)
Bluebell1_scaled = (Bluebell1_interpolated$GR - Blmean)/Blstd
Bluebell1_rescaled = data.frame(Bluebell1_interpolated$DEPT, Bluebell1_scaled)

# Resampling the data using moving window statistics
Fisher1_scaled = mwStats(Fisher1_rescaled, cols = 2, win=3, ends = T)
Fisher1_standardized = data.frame(Fisher1_scaled$Center_win, Fisher1_scaled$Average)

Bluebell1_scaled = mwStats(Bluebell1_rescaled, cols = 2, win=3, ends = T)
Bluebell1_standardized = data.frame(Bluebell1_scaled$Center_win, Bluebell1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Fisher1_standardized, type="l", xlim = c(100, 1700), ylim = c(-20, 20), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
plot(Bluebell1_standardized, type="l", xlim = c(200, 2000), ylim = c(-20, 20), xlab = "Bluebell1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_bl1_fi1_ap1 <- dtw(Bluebell1_standardized$Bluebell1_scaled.Average, Fisher1_standardized$Fisher1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_bl1_fi1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Bluebell1_on_Fisher1_depth = tune(Bluebell1_standardized, cbind(Bluebell1_standardized$Bluebell1_scaled.Center_win[al_bl1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_bl1_fi1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Fisher1_standardized, type = "l", ylim = c(-20, 20), xlim = c(100, 1700), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
lines(Bluebell1_on_Fisher1_depth, col = "red")

# DTW Distance
al_bl1_fi1_ap1$normalizedDistance
al_bl1_fi1_ap1$distance

# Tuning Bluebell1 data on Finucane1 depth
Bluebell1_on_Finucane1_depth = tune(Bluebell1_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning Bluebell1 data on Picard1 depth
Bluebell1_on_Picard1_depth = tune(Bluebell1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Bluebell-1)")
lines(Bluebell1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Bluebell1_originalGR_on_Picard1_depth = data.frame(Bluebell1_on_Picard1_depth$X1, Bluebell1_interpolated[1:8746,2])

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Bluebell-1)")
lines(Bluebell1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Bluebell1 

U1463Age_on_Bluebell1_depth = tune(Bluebell1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Bluebell1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Bluebell1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Bluebell1_depth) <- new_column_names
write.csv(U1463Age_on_Bluebell1_depth, file = "C:/Users/Rohit/OneDrive - Universität Münster/Maps/Base Map/U1463_Age-Site_GR/Bluebell 1.csv", row.names = FALSE)

Bluebell1_age_depth <- approx(x = Bluebell1_interpolated$GR,
                             y = Bluebell1_interpolated$DEPT,
                             xout = U1463Age_on_Bluebell1_depth$GR,
                             rule = 1)$y

Bluebell1_agemodel = data.frame(Bluebell1_age_depth, U1463Age_on_Bluebell1_depth$AGE)
Bluebell1_agemodel <- na.omit(Bluebell1_agemodel)
new_column_names1 <- c("Depth", "Age")
colnames(Bluebell1_agemodel) <- new_column_names1
plot(Bluebell1_agemodel)
write.csv(Bluebell1_agemodel, file = "RScripts&Data/Sites Data_Age-NGR/Bluebell1_DepthAge.csv", row.names = FALSE)
