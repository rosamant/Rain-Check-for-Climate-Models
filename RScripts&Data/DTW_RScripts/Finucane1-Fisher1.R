install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Finucane1 and Fisher1 datasets

Finucane1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Finucane1.csv", header=TRUE, stringsAsFactors=FALSE)
Finucane1=Finucane1[c(1:9080),] # Oligocene-Miocene
head(Finucane1)
plot(Finucane1, type="l", xlim = c(150, 1550), ylim = c(0, 50))

Fisher1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Fisher1.csv", header=TRUE, stringsAsFactors=FALSE)
Fisher1=Fisher1[c(1:10190),] # Oligocene-Miocene
head(Fisher1)
plot(Fisher1, type="l", xlim = c(100, 1700), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Finucane1_interpolated <- linterp(Finucane1, dt = 0.2, genplot = F)
Fisher1_interpolated <- linterp(Fisher1, dt = 0.2, genplot = F)

# Scaling the data
Fmean = Gmean(Finucane1_interpolated$GR)
Fstd = Gsd(Finucane1_interpolated$GR)
Finucane1_scaled = (Finucane1_interpolated$GR - Fmean)/Fstd
Finucane1_rescaled = data.frame(Finucane1_interpolated$DEPT, Finucane1_scaled)

Fimean = Gmean(Fisher1_interpolated$GR)
Fistd = Gsd(Fisher1_interpolated$GR)
Fisher1_scaled = (Fisher1_interpolated$GR - Fimean)/Fistd
Fisher1_rescaled = data.frame(Fisher1_interpolated$DEPT, Fisher1_scaled)

# Resampling the data using moving window statistics
Finucane1_scaled = mwStats(Finucane1_rescaled, cols = 2, win=3, ends = T)
Finucane1_standardized = data.frame(Finucane1_scaled$Center_win, Finucane1_scaled$Average)

Fisher1_scaled = mwStats(Fisher1_rescaled, cols = 2, win=3, ends = T)
Fisher1_standardized = data.frame(Fisher1_scaled$Center_win, Fisher1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Finucane1_standardized, type="l", xlim = c(150, 1550), ylim = c(-20, 20), xlab = "Finucane1 Resampled Depth", ylab = "Normalized GR")
plot(Fisher1_standardized, type="l", xlim = c(100, 1700), ylim = c(-20, 20), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_fi1_f1_ap1 <- dtw(Fisher1_standardized$Fisher1_scaled.Average, Finucane1_standardized$Finucane1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_fi1_f1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Fisher1_on_Finucane1_depth = tune(Fisher1_standardized, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Finucane1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1550), xlab = "Athena1 Resampled Depth", ylab = "Normalized GR")
lines(Fisher1_on_Finucane1_depth, col = "red")

# DTW Distance
al_fi1_f1_ap1$normalizedDistance
al_fi1_f1_ap1$distance

# Tuning Fisher1 data on Picard1 depth
Fisher1_on_Picard1_depth = tune(Fisher1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Fisher-1)")
lines(Fisher1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Fisher1_originalGR_on_Picard1_depth = data.frame(Fisher1_on_Picard1_depth$X1, Fisher1_interpolated[1:8105,2])

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Fisher-1)")
lines(Fisher1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Fisher1 

U1463Age_on_Fisher1_depth = tune(Fisher1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Fisher1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Fisher1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Fisher1_depth) <- new_column_names
write.csv(U1463Age_on_Fisher1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Fisher 1.csv", row.names = FALSE)


Fisher1_age_depth <- approx(x = Fisher1_interpolated$GR,
                          y = Fisher1_interpolated$DEPT,
                          xout = U1463Age_on_Fisher1_depth$GR,
                          rule = 1)$y

Fisher1_agemodel = data.frame(Fisher1_age_depth, U1463Age_on_Fisher1_depth$AGE)
Fisher1_agemodel <- na.omit(Fisher1_agemodel)
new_column_names1 <- c("Depth", "Age")
colnames(Fisher1_agemodel) <- new_column_names1
plot(Fisher1_agemodel)
write.csv(Fisher1_agemodel, file = "RScripts&Data/Sites Data_Age-NGR/Fisher1_DepthAge.csv", row.names = FALSE)
