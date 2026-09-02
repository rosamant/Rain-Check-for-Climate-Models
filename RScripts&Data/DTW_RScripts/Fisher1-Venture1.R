install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Fisher1 and Venture1 datasets

Fisher1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Fisher1.csv", header=TRUE, stringsAsFactors=FALSE)
Fisher1=Fisher1[c(1:10190),] # Oligocene-Miocene
head(Fisher1)
plot(Fisher1, type="l", xlim = c(100, 1700), ylim = c(0, 50))

Venture1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Venture1.csv", header=TRUE, stringsAsFactors=FALSE)

# Recorrecting attenuated signal
V1 = Gmean(Venture1[c(1:1766),2])
V2 = Gmean(Venture1[c(800:1500),2])
SD1 = Gsd(Venture1[c(1:1766),2])
SD2 = Gsd(Venture1[c(800:1500),2])
Venture1[c(1:1766),2]=(Venture1[c(1:1766),2]+(V2-V1))*(SD1/SD2)

Venture1=Venture1[c(1:7223),] # Oligocene-Miocene
head(Venture1)
plot(Venture1, type="l", xlim = c(100, 1550), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Fisher1_interpolated <- linterp(Fisher1, dt = 0.2, genplot = F)
Venture1_interpolated <- linterp(Venture1, dt = 0.2, genplot = F)

# Scaling the data
Fimean = Gmean(Fisher1_interpolated$GR)
Fistd = Gsd(Fisher1_interpolated$GR)
Fisher1_scaled = (Fisher1_interpolated$GR - Fimean)/Fistd
Fisher1_rescaled = data.frame(Fisher1_interpolated$DEPT, Fisher1_scaled)

Vmean = Gmean(Venture1_interpolated$GR)
Vstd = Gsd(Venture1_interpolated$GR)
Venture1_scaled = (Venture1_interpolated$GR - Vmean)/Vstd
Venture1_rescaled = data.frame(Venture1_interpolated$DEPT, Venture1_scaled)

# Resampling the data using moving window statistics
Fisher1_scaled = mwStats(Fisher1_rescaled, cols = 2, win=3, ends = T)
Fisher1_standardized = data.frame(Fisher1_scaled$Center_win, Fisher1_scaled$Average)

Venture1_scaled = mwStats(Venture1_rescaled, cols = 2, win=3, ends = T)
Venture1_standardized = data.frame(Venture1_scaled$Center_win, Venture1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Fisher1_standardized, type="l", xlim = c(100, 1700), ylim = c(-20, 20), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
plot(Venture1_standardized, type="l", xlim = c(100, 1550), ylim = c(-20, 20), xlab = "Venture1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_v1_fi1_ap1 <- dtw(Venture1_standardized$Venture1_scaled.Average, Fisher1_standardized$Fisher1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_v1_fi1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Venture1_on_Fisher1_depth = tune(Venture1_standardized, cbind(Venture1_standardized$Venture1_scaled.Center_win[al_v1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_v1_fi1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Fisher1_standardized, type = "l", ylim = c(-20, 20), xlim = c(100, 1700), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
lines(Venture1_on_Fisher1_depth, col = "red")

# DTW Distance
al_v1_fi1_ap1$normalizedDistance
al_v1_fi1_ap1$distance

# Tuning Venture1 data on Finucane1 depth
Venture1_on_Finucane1_depth = tune(Venture1_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning Venture1 data on Picard1 depth
Venture1_on_Picard1_depth = tune(Venture1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Venture-1)")
lines(Venture1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Venture1_originalGR_on_Picard1_depth = data.frame(Venture1_on_Picard1_depth$X1, Venture1_interpolated$GR)

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Venture-1)")
lines(Venture1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Venture1 

U1463Age_on_Venture1_depth = tune(Venture1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Venture1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Venture1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Venture1_depth) <- new_column_names
write.csv(U1463Age_on_Venture1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Venture 1.csv", row.names = FALSE)


Venture1_age_depth <- approx(x = Venture1_interpolated$GR,
                                y = Venture1_interpolated$DEPT,
                                xout = U1463Age_on_Venture1_depth$GR,
                                rule = 1)$y

Venture1_agemodel = data.frame(Venture1_age_depth, U1463Age_on_Venture1_depth$AGE)
Venture1_agemodel <- na.omit(Venture1_agemodel)
new_column_names1 <- c("Depth", "Age")
colnames(Venture1_agemodel) <- new_column_names1
plot(Venture1_agemodel)
write.csv(Venture1_agemodel, file = "RScripts&Data/Sites Data_Age-NGR/Venture1_DepthAge.csv", row.names = FALSE)
