install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Kalyptea 1 and GreatEastern 1 datasets

Kalyptea1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Kalyptea 1.csv", header=TRUE, stringsAsFactors=FALSE)
Kalyptea1=Kalyptea1[c(110:6429),] # Oligocene-Miocene
head(Kalyptea1)
plot(Kalyptea1, type="l", xlim = c(230, 1200), ylim = c(0, 60))

GreatEastern1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Great Eastern 1.csv", header=TRUE, stringsAsFactors=FALSE)
GreatEastern1=GreatEastern1[c(1:6633),] # Oligocene-Miocene
head(GreatEastern1)
plot(GreatEastern1, type="l", xlim = c(150, 1150), ylim = c(0, 60))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Kalyptea1_interpolated <- linterp(Kalyptea1, dt = 0.2, genplot = F)
GreatEastern1_interpolated <- linterp(GreatEastern1, dt = 0.2, genplot = F)

# Scaling the data
Kamean = Gmean(Kalyptea1_interpolated$GR)
Kastd = Gsd(Kalyptea1_interpolated$GR)
Kalyptea1_scaled = (Kalyptea1_interpolated$GR - Kamean)/Kastd
Kalyptea1_rescaled = data.frame(Kalyptea1_interpolated$DEPT, Kalyptea1_scaled)

GEmean = Gmean(GreatEastern1_interpolated$GR)
GEstd = Gsd(GreatEastern1_interpolated$GR)
GreatEastern1_scaled = (GreatEastern1_interpolated$GR - GEmean)/GEstd
GreatEastern1_rescaled = data.frame(GreatEastern1_interpolated$DEPT, GreatEastern1_scaled)

# Resampling the data using moving window statistics
Kalyptea1_scaled = mwStats(Kalyptea1_rescaled, cols = 2, win=3, ends = T)
Kalyptea1_standardized = data.frame(Kalyptea1_scaled$Center_win, Kalyptea1_scaled$Average)

GreatEastern1_scaled = mwStats(GreatEastern1_rescaled, cols = 2, win=3, ends = T)
GreatEastern1_standardized = data.frame(GreatEastern1_scaled$Center_win, GreatEastern1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Kalyptea1_standardized, type="l", xlim = c(230, 1200), ylim = c(-20, 20), xlab = "Kalyptea1 Resampled Depth", ylab = "Normalized GR")
plot(GreatEastern1_standardized, type="l", xlim = c(150, 1150), ylim = c(-10, 20), xlab = "GreatEastern 1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_ge1_ka1_ap1 <- dtw(GreatEastern1_standardized$GreatEastern1_scaled.Average, Kalyptea1_standardized$Kalyptea1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_ge1_ka1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
GreatEastern1_on_Kalyptea1_depth = tune(GreatEastern1_standardized, cbind(GreatEastern1_standardized$GreatEastern1_scaled.Center_win[al_ge1_ka1_ap1$index1s], Kalyptea1_standardized$Kalyptea1_scaled.Center_win[al_ge1_ka1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Kalyptea1_standardized, type = "l", ylim = c(-20, 20), xlim = c(230, 1200), xlab = "Kalyptea1 Resampled Depth", ylab = "Normalized GR")
lines(GreatEastern1_on_Kalyptea1_depth, col = "red")

# DTW Distance
al_ge1_ka1_ap1$normalizedDistance
al_ge1_ka1_ap1$distance

# Tuning GreatEastern1 data on Prelude1 depth
GreatEastern1_on_Prelude1_depth = tune(GreatEastern1_on_Kalyptea1_depth, cbind(Kalyptea1_standardized$Kalyptea1_scaled.Center_win[al_ka1_p1_ap1$index1s], Prelude1_standardized$Prelude1_scaled.Center_win[al_ka1_p1_ap1$index2s]), extrapolate = F)

# Tuning GreatEastern1 data on Gorgonichthys1 depth
GreatEastern1_on_Gorgonichthys1_depth = tune(GreatEastern1_on_Prelude1_depth, cbind(Prelude1_standardized$Prelude1_scaled.Center_win[al_p1_g1_ap1$index1s], Gorgonichthys1_standardized$Gorgonichthys1_scaled.Center_win[al_p1_g1_ap1$index2s]), extrapolate = F)

# Tuning GreatEastern1 data on Walkley1 depth
GreatEastern1_on_Walkley1_depth = tune(GreatEastern1_on_Gorgonichthys1_depth, cbind(Gorgonichthys1_standardized$Gorgonichthys1_scaled.Center_win[al_g1_w1_ap2$index1s], Walkley1_standardized$Walkley1_scaled.Center_win[al_g1_w1_ap2$index2s]), extrapolate = F)

# Tuning GreatEastern1 data on Caswell1 depth
GreatEastern1_on_Caswell1_depth = tune(GreatEastern1_on_Walkley1_depth, cbind(Walkley1_standardized$Walkley1_scaled.Center_win[al_w1_c1_ap1$index1s], Caswell1_standardized$Caswell1_scaled.Center_win[al_w1_c1_ap1$index2s]), extrapolate = F)

# Tuning GreatEastern1 data on Calliance2 depth
GreatEastern1_on_Calliance2_depth = tune(GreatEastern1_on_Caswell1_depth, cbind(Caswell1_standardized$Caswell1_scaled.Center_win[al_c2_c1_ap1$index2s], Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_c1_ap1$index1s]), extrapolate = F)

# Tuning GreatEastern1 data on Omar1 depth
GreatEastern1_on_Omar1_depth = tune(GreatEastern1_on_Calliance2_depth, cbind(Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_o1_ap2$index1s], Omar1_standardized$Omar1_scaled.Center_win[al_c2_o1_ap2$index2s]), extrapolate = F)

# Tuning GreatEastern1 data on SG1 depth
GreatEastern1_on_SG1_depth = tune(GreatEastern1_on_Omar1_depth, cbind(Omar1_standardized$Omar1_scaled.Center_win[al_o1_sg1_ap2$index1s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_o1_sg1_ap2$index2s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 20), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (GreatEastern-1)")
lines(GreatEastern1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
GreatEastern1_originalGR_on_SG1_depth = data.frame(GreatEastern1_on_SG1_depth$X1, GreatEastern1_interpolated[1949:5052,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (GreatEastern-1)")
lines(GreatEastern1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to GreatEastern1 

SG1Age_on_GreatEastern1_depth = tune(GreatEastern1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_GreatEastern1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "GreatEastern1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_GreatEastern1_depth) <- new_column_names
SG1Age_on_GreatEastern1_depth[,1] = SG1Age_on_GreatEastern1_depth[,1] * 1000
write.csv(SG1Age_on_GreatEastern1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Great Eastern 1.csv", row.names = FALSE)

GreatEastern1_age_depth <- approx(x = GreatEastern1_interpolated$GR,
                              y = GreatEastern1_interpolated$DEPT,
                              xout = SG1Age_on_GreatEastern1_depth$GR,
                              rule = 1)$y

GreatEastern1_agemodel = data.frame(GreatEastern1_age_depth, SG1Age_on_GreatEastern1_depth$AGE)
GreatEastern1_agemodel <- na.omit(GreatEastern1_agemodel)
new_column_names1 <- c("Depth", "Age")
colnames(GreatEastern1_agemodel) <- new_column_names1
plot(GreatEastern1_agemodel)
write.csv(GreatEastern1_agemodel, file = "RScripts&Data/Sites Data_Age-NGR/GreatEastern1_DepthAge.csv", row.names = FALSE)
