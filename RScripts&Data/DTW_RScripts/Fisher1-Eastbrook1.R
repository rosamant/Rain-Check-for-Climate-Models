install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Fisher1 and Eastbrook1 datasets

Fisher1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Fisher1.csv", header=TRUE, stringsAsFactors=FALSE)
Fisher1=Fisher1[c(1:10190),] # Oligocene-Miocene
head(Fisher1)
plot(Fisher1, type="l", xlim = c(100, 1700), ylim = c(0, 50))

Eastbrook1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Eastbrook1.csv", header=TRUE, stringsAsFactors=FALSE)

# Recorrecting attenuated signal
E1 = Gmean(Eastbrook1[c(3000:4200),2])
E2 = Gmean(Eastbrook1[c(4292:4915),2])
SD1 = Gsd(Eastbrook1[c(3000:4200),2])
SD2 = Gsd(Eastbrook1[c(4292:4915),2])
Eastbrook1[c(4292:4915),2]=(Eastbrook1[c(4292:4915),2]+(E1-E2))*(SD2/SD1)

Eastbrook1=Eastbrook1[c(1:12029),] # Oligocene-Miocene
head(Eastbrook1)
plot(Eastbrook1, type="l", xlim = c(200, 2100), ylim = c(0, 50))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Fisher1_interpolated <- linterp(Fisher1, dt = 0.2, genplot = F)
Eastbrook1_interpolated <- linterp(Eastbrook1, dt = 0.2, genplot = F)

# Scaling the data
Fimean = Gmean(Fisher1_interpolated$GR)
Fistd = Gsd(Fisher1_interpolated$GR)
Fisher1_scaled = (Fisher1_interpolated$GR - Fimean)/Fistd
Fisher1_rescaled = data.frame(Fisher1_interpolated$DEPT, Fisher1_scaled)

Emean = Gmean(Eastbrook1_interpolated$GR)
Estd = Gsd(Eastbrook1_interpolated$GR)
Eastbrook1_scaled = (Eastbrook1_interpolated$GR - Emean)/Estd
Eastbrook1_rescaled = data.frame(Eastbrook1_interpolated$DEPT, Eastbrook1_scaled)

# Resampling the data using moving window statistics
Fisher1_scaled = mwStats(Fisher1_rescaled, cols = 2, win=3, ends = T)
Fisher1_standardized = data.frame(Fisher1_scaled$Center_win, Fisher1_scaled$Average)

Eastbrook1_scaled = mwStats(Eastbrook1_rescaled, cols = 2, win=3, ends = T)
Eastbrook1_standardized = data.frame(Eastbrook1_scaled$Center_win, Eastbrook1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Fisher1_standardized, type="l", xlim = c(100, 1700), ylim = c(-20, 20), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
plot(Eastbrook1_standardized, type="l", xlim = c(200, 2000), ylim = c(-20, 20), xlab = "Eastbrook1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_e1_fi1_ap1 <- dtw(Eastbrook1_standardized$Eastbrook1_scaled.Average, Fisher1_standardized$Fisher1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_e1_fi1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Eastbrook1_on_Fisher1_depth = tune(Eastbrook1_standardized, cbind(Eastbrook1_standardized$Eastbrook1_scaled.Center_win[al_e1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_e1_fi1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Fisher1_standardized, type = "l", ylim = c(-20, 20), xlim = c(100, 1700), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
lines(Eastbrook1_on_Fisher1_depth, col = "red")

# DTW Distance
al_e1_fi1_ap1$normalizedDistance
al_e1_fi1_ap1$distance

# Tuning Eastbrook1 data on Finucane1 depth
Eastbrook1_on_Finucane1_depth = tune(Eastbrook1_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning Eastbrook1 data on Picard1 depth
Eastbrook1_on_Picard1_depth = tune(Eastbrook1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Eastbrook-1)")
lines(Eastbrook1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Eastbrook1_originalGR_on_Picard1_depth = data.frame(Eastbrook1_on_Picard1_depth$X1, Eastbrook1_interpolated[1:9166,2])

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Eastbrook-1)")
lines(Eastbrook1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Eastbrook1 

U1463Age_on_Eastbrook1_depth = tune(Eastbrook1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Eastbrook1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Eastbrook1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Eastbrook1_depth) <- new_column_names
write.csv(U1463Age_on_Eastbrook1_depth, file = "C:/Users/Rohit/OneDrive - Universität Münster/Maps/Base Map/U1463_Age-Site_GR/Eastbrook 1.csv", row.names = FALSE)


Eastbrook1_age_depth <- approx(x = Eastbrook1_interpolated$GR,
                                   y = Eastbrook1_interpolated$DEPT,
                                   xout = U1463Age_on_Eastbrook1_depth$GR,
                                   rule = 1)$y

Eastbrook1_agemodel = data.frame(Eastbrook1_age_depth, U1463Age_on_Eastbrook1_depth$AGE)
Eastbrook1_agemodel <- na.omit(Eastbrook1_agemodel)
new_column_names1 <- c("Depth", "Age")
colnames(Eastbrook1_agemodel) <- new_column_names1
plot(Eastbrook1_agemodel)
write.csv(Eastbrook1_agemodel, file = "RScripts&Data/Sites Data_Age-NGR/Eastbrook1_DepthAge.csv", row.names = FALSE)
