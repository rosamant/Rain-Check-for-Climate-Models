install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Hadrian 1 and Mandorah 1 datasets

Hadrian1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Hadrian 1.csv", header=TRUE, stringsAsFactors=FALSE)
Hadrian1=Hadrian1[c(28:6786),] # Oligocene-Miocene
head(Hadrian1)
plot(Hadrian1, type="l", xlim = c(300, 1350), ylim = c(0, 60))

Mandorah1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Mandorah 1.csv", header=TRUE, stringsAsFactors=FALSE)
Mandorah1=Mandorah1[c(1:14315),] # Oligocene-Miocene
head(Mandorah1)
plot(Mandorah1, type="l", xlim = c(370, 1850), ylim = c(0, 60))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Hadrian1_interpolated <- linterp(Hadrian1, dt = 0.2, genplot = F)
Mandorah1_interpolated <- linterp(Mandorah1, dt = 0.2, genplot = F)

# Scaling the data
Hmean = Gmean(Hadrian1_interpolated$GR)
Hstd = Gsd(Hadrian1_interpolated$GR)
Hadrian1_scaled = (Hadrian1_interpolated$GR - Hmean)/Hstd
Hadrian1_rescaled = data.frame(Hadrian1_interpolated$DEPT, Hadrian1_scaled)

Mhmean = Gmean(Mandorah1_interpolated$GR)
Mhstd = Gsd(Mandorah1_interpolated$GR)
Mandorah1_scaled = (Mandorah1_interpolated$GR - Mhmean)/Mhstd
Mandorah1_rescaled = data.frame(Mandorah1_interpolated$DEPT, Mandorah1_scaled)

# Resampling the data using moving window statistics
Hadrian1_scaled = mwStats(Hadrian1_rescaled, cols = 2, win=3, ends = T)
Hadrian1_standardized = data.frame(Hadrian1_scaled$Center_win, Hadrian1_scaled$Average)

Mandorah1_scaled = mwStats(Mandorah1_rescaled, cols = 2, win=3, ends = T)
Mandorah1_standardized = data.frame(Mandorah1_scaled$Center_win, Mandorah1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Hadrian1_standardized, type="l", xlim = c(300, 1350), ylim = c(-20, 20), xlab = "Hadrian1 Resampled Depth", ylab = "Normalized GR")
plot(Mandorah1_standardized, type="l", xlim = c(370, 1850), ylim = c(-10, 20), xlab = "Mandorah 1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_m1_h1_ap1 <- dtw(Mandorah1_standardized$Mandorah1_scaled.Average, Hadrian1_standardized$Hadrian1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_m1_h1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Mandorah1_on_Hadrian1_depth = tune(Mandorah1_standardized, cbind(Mandorah1_standardized$Mandorah1_scaled.Center_win[al_m1_h1_ap1$index1s], Hadrian1_standardized$Hadrian1_scaled.Center_win[al_m1_h1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Hadrian1_standardized, type = "l", ylim = c(-20, 20), xlim = c(300, 1350), xlab = "Hadrian1 Resampled Depth", ylab = "Normalized GR")
lines(Mandorah1_on_Hadrian1_depth, col = "red")

# DTW Distance
al_m1_h1_ap1$normalizedDistance
al_m1_h1_ap1$distance

# Tuning Mandorah1 data on SG1 depth
Mandorah1_on_Brontosaurus1_depth = tune(Mandorah1_on_Hadrian1_depth, cbind(Hadrian1_standardized$Hadrian1_scaled.Center_win[al_h1_b1_ap1$index1s], Brontosaurus1_standardized$Brontosaurus1_scaled.Center_win[al_h1_b1_ap1$index2s]), extrapolate = F)

# Tuning Mandorah1 data on SG1 depth
Mandorah1_on_SG1_depth = tune(Mandorah1_on_Brontosaurus1_depth, cbind(Brontosaurus11_standardized$Brontosaurus11_scaled.Center_win[al_sg1_b1_ap2$index2s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_sg1_b1_ap2$index1s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 40), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (Mandorah-1)")
lines(Mandorah1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Mandorah1_originalGR_on_SG1_depth = data.frame(Mandorah1_on_SG1_depth$X1, Mandorah1_interpolated[1672:7312,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (Mandorah-1)")
lines(Mandorah1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Mandorah1 

SG1Age_on_Mandorah1_depth = tune(Mandorah1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Mandorah1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Mandorah-1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Mandorah1_depth) <- new_column_names
SG1Age_on_Mandorah1_depth[,1] = SG1Age_on_Mandorah1_depth[,1] * 1000
write.csv(SG1Age_on_Mandorah1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Mandorah 1.csv", row.names = FALSE)

Mandorah1_age_depth <- approx(x = Mandorah1_interpolated$GR,
                             y = Mandorah1_interpolated$DEPT,
                             xout = SG1Age_on_Mandorah1_depth$GR,
                             rule = 1)$y

Mandorah1_agemodel = data.frame(Mandorah1_age_depth, SG1Age_on_Mandorah1_depth$AGE)
Mandorah1_agemodel <- na.omit(Mandorah1_agemodel)
new_column_names1 <- c("Depth", "Age")
colnames(Mandorah1_agemodel) <- new_column_names1
plot(Mandorah1_agemodel)
write.csv(Mandorah1_agemodel, file = "RScripts&Data/Sites Data_Age-NGR/Mandorah1_DepthAge.csv", row.names = FALSE)
