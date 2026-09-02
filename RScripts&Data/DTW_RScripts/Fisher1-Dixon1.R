install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Fisher1 and Dixon1 datasets

Fisher1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Fisher1.csv", header=TRUE, stringsAsFactors=FALSE)
Fisher1=Fisher1[c(1:10190),] # Oligocene-Miocene
head(Fisher1)
plot(Fisher1, type="l", xlim = c(100, 1700), ylim = c(0, 50))

Dixon1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Dixon1.csv", header=TRUE, stringsAsFactors=FALSE)
Dixon1=Dixon1[c(1:10220),] # Oligocene-Miocene
head(Dixon1)
plot(Dixon1, type="l", xlim = c(100, 1700), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Fisher1_interpolated <- linterp(Fisher1, dt = 0.2, genplot = F)
Dixon1_interpolated <- linterp(Dixon1, dt = 0.2, genplot = F)

# Scaling the data
Fimean = Gmean(Fisher1_interpolated$GR)
Fistd = Gsd(Fisher1_interpolated$GR)
Fisher1_scaled = (Fisher1_interpolated$GR - Fimean)/Fistd
Fisher1_rescaled = data.frame(Fisher1_interpolated$DEPT, Fisher1_scaled)

Dmean = Gmean(Dixon1_interpolated$GR)
Dstd = Gsd(Dixon1_interpolated$GR)
Dixon1_scaled = (Dixon1_interpolated$GR - Dmean)/Dstd
Dixon1_rescaled = data.frame(Dixon1_interpolated$DEPT, Dixon1_scaled)

# Resampling the data using moving window statistics
Fisher1_scaled = mwStats(Fisher1_rescaled, cols = 2, win=3, ends = T)
Fisher1_standardized = data.frame(Fisher1_scaled$Center_win, Fisher1_scaled$Average)

Dixon1_scaled = mwStats(Dixon1_rescaled, cols = 2, win=3, ends = T)
Dixon1_standardized = data.frame(Dixon1_scaled$Center_win, Dixon1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Fisher1_standardized, type="l", xlim = c(100, 1700), ylim = c(-20, 20), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
plot(Dixon1_standardized, type="l", xlim = c(100, 1700), ylim = c(-20, 20), xlab = "Dixon1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_d1_fi1_ap1 <- dtw(Dixon1_standardized$Dixon1_scaled.Average, Fisher1_standardized$Fisher1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = F, open.end = T))
plot(al_d1_fi1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Dixon1_on_Fisher1_depth = tune(Dixon1_standardized, cbind(Dixon1_standardized$Dixon1_scaled.Center_win[al_d1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_d1_fi1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Fisher1_standardized, type = "l", ylim = c(-20, 20), xlim = c(100, 1700), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
lines(Dixon1_on_Fisher1_depth, col = "red")

# DTW Distance
al_d1_fi1_ap1$normalizedDistance
al_d1_fi1_ap1$distance

# Tuning Dixon1 data on Finucane1 depth
Dixon1_on_Finucane1_depth = tune(Dixon1_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning Dixon1 data on Picard1 depth
Dixon1_on_Picard1_depth = tune(Dixon1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Dixon-1)")
lines(Dixon1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Dixon1_originalGR_on_Picard1_depth = data.frame(Dixon1_on_Picard1_depth$X1, Dixon1_interpolated[1:7786,2])

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Dixon-1)")
lines(Dixon1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Dixon1 

U1463Age_on_Dixon1_depth = tune(Dixon1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Dixon1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Dixon1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Dixon1_depth) <- new_column_names
write.csv(U1463Age_on_Dixon1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Dixon 1.csv", row.names = FALSE)
