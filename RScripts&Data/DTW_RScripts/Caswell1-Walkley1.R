install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Caswell1 and Walkley1 datasets

Caswell1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Caswell1.csv", header=TRUE, stringsAsFactors=FALSE)
Caswell1=Caswell1[c(1:6434),] # Oligocene-Miocene
head(Caswell1)
plot(Caswell1, type="l", xlim = c(350, 1750), ylim = c(0, 50))

Walkley1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Walkley 1.csv", header=TRUE, stringsAsFactors=FALSE)
Walkley1=Walkley1[c(6:6312),] # Oligocene-Miocene
head(Walkley1)
plot(Walkley1, type="l", xlim = c(350, 1600), ylim = c(0, 60))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Caswell1_interpolated <- linterp(Caswell1, dt = 0.2, genplot = F)
Walkley1_interpolated <- linterp(Walkley1, dt = 0.2, genplot = F)

# Scaling the data
Cmean = Gmean(Caswell1_interpolated$GR)
Cstd = Gsd(Caswell1_interpolated$GR)
Caswell1_scaled = (Caswell1_interpolated$GR - Cmean)/Cstd
Caswell1_rescaled = data.frame(Caswell1_interpolated$DEPT, Caswell1_scaled)

Wmean = Gmean(Walkley1_interpolated$GR)
Wstd = Gsd(Walkley1_interpolated$GR)
Walkley1_scaled = (Walkley1_interpolated$GR - Wmean)/Wstd
Walkley1_rescaled = data.frame(Walkley1_interpolated$DEPT, Walkley1_scaled)

# Resampling the data using moving window statistics
Caswell1_scaled = mwStats(Caswell1_rescaled, cols = 2, win=3, ends = T)
Caswell1_standardized = data.frame(Caswell1_scaled$Center_win, Caswell1_scaled$Average)

Walkley1_scaled = mwStats(Walkley1_rescaled, cols = 2, win=3, ends = T)
Walkley1_standardized = data.frame(Walkley1_scaled$Center_win, Walkley1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Caswell1_standardized, type="l", xlim = c(350, 1700), ylim = c(-20, 20), xlab = "Caswell1 Resampled Depth", ylab = "Normalized GR")
plot(Walkley1_standardized, type="l", xlim = c(350, 1600), ylim = c(-20, 20), xlab = "Walkley1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_w1_c1_ap1 <- dtw(Walkley1_standardized$Walkley1_scaled.Average, Caswell1_standardized$Caswell1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_w1_c1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Walkley1_on_Caswell1_depth = tune(Walkley1_standardized, cbind(Walkley1_standardized$Walkley1_scaled.Center_win[al_w1_c1_ap1$index1s], Caswell1_standardized$Caswell1_scaled.Center_win[al_w1_c1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Caswell1_standardized, type = "l", ylim = c(-20, 20), xlim = c(350, 1700), xlab = "Caswell1 Resampled Depth", ylab = "Normalized GR")
lines(Walkley1_on_Caswell1_depth, col = "red")

# DTW Distance

al_w1_c1_ap1$normalizedDistance
al_w1_c1_ap1$distance

# Tuning Walkley1 data on Calliance2 depth
Walkley1_on_Calliance2_depth = tune(Walkley1_on_Caswell1_depth, cbind(Caswell1_standardized$Caswell1_scaled.Center_win[al_c2_c1_ap1$index2s], Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_c1_ap1$index1s]), extrapolate = F)

# Tuning Walkley1 data on Omar1 depth
Walkley1_on_Omar1_depth = tune(Walkley1_on_Calliance2_depth, cbind(Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_o1_ap2$index1s], Omar1_standardized$Omar1_scaled.Center_win[al_c2_o1_ap2$index2s]), extrapolate = F)

# Tuning Walkley1 data on SG1 depth
Walkley1_on_SG1_depth = tune(Walkley1_on_Omar1_depth, cbind(Omar1_standardized$Omar1_scaled.Center_win[al_o1_sg1_ap2$index1s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_o1_sg1_ap2$index2s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 20), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (Walkley-1)")
lines(Walkley1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Walkley1_originalGR_on_SG1_depth = data.frame(Walkley1_on_SG1_depth$X1, Walkley1_interpolated[1527:6142,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (Walkley-1)")
lines(Walkley1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Walkley1 

SG1Age_on_Walkley1_depth = tune(Walkley1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Walkley1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Walkley1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Walkley1_depth) <- new_column_names
SG1Age_on_Walkley1_depth[,1] = SG1Age_on_Walkley1_depth[,1] * 1000
write.csv(SG1Age_on_Walkley1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Walkley 1.csv", row.names = FALSE)


Walkley1_age_depth <- approx(x = Walkley1_interpolated$GR,
                          y = Walkley1_interpolated$DEPT,
                          xout = SG1Age_on_Walkley1_depth$GR,
                          rule = 1)$y

Walkley1_agemodel = data.frame(Walkley1_age_depth, SG1Age_on_Walkley1_depth$AGE)
Walkley1_agemodel <- na.omit(Walkley1_agemodel)
new_column_names1 <- c("Depth", "Age")
colnames(Walkley1_agemodel) <- new_column_names1
plot(Walkley1_agemodel)
write.csv(Walkley1_agemodel, file = "RScripts&Data/Sites Data_Age-NGR/Walkley1_DepthAge.csv", row.names = FALSE)
