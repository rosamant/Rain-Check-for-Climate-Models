install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Hyde1 and Ermine1 datasets

Hyde1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Hyde1.csv", header=TRUE, stringsAsFactors=FALSE)

# Recorrecting attenuated signal
H1 = Gmean(Hyde1[c(1:300),2])
H2 = Gmean(Hyde1[c(300:700),2])
SD1 = Gsd(Hyde1[c(1:300),2])
SD2 = Gsd(Hyde1[c(300:700),2])
Hyde1[c(1:300),2]=(Hyde1[c(1:300),2]+(H2-H1))*(SD1/SD2)

Hyde1=Hyde1[c(1:8335),] # Oligocene-Miocene
head(Hyde1)
plot(Hyde1, type="l", xlim = c(400, 2150), ylim = c(0, 60))

Ermine1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Ermine1.csv", header=TRUE, stringsAsFactors=FALSE)
Ermine1=Ermine1[c(559:14259),] # Oligocene-Miocene
head(Ermine1)
plot(Ermine1, type="l", xlim = c(600, 2000), ylim = c(0, 70))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Hyde1_interpolated <- linterp(Hyde1, dt = 0.2, genplot = F)
Ermine1_interpolated <- linterp(Ermine1, dt = 0.2, genplot = F)

# Scaling the data
Hmean = Gmean(Hyde1_interpolated$GR)
Hstd = Gsd(Hyde1_interpolated$GR)
Hyde1_scaled = (Hyde1_interpolated$GR - Hmean)/Hstd
Hyde1_rescaled = data.frame(Hyde1_interpolated$DEPT, Hyde1_scaled)

Emean = Gmean(Ermine1_interpolated$GR)
Estd = Gsd(Ermine1_interpolated$GR)
Ermine1_scaled = (Ermine1_interpolated$GR - Emean)/Estd
Ermine1_rescaled = data.frame(Ermine1_interpolated$DEPT, Ermine1_scaled)

# Resampling the data using moving window statistics
Hyde1_scaled = mwStats(Hyde1_rescaled, cols = 2, win=3, ends = T)
Hyde1_standardized = data.frame(Hyde1_scaled$Center_win, Hyde1_scaled$Average)

Ermine1_scaled = mwStats(Ermine1_rescaled, cols = 2, win=3, ends = T)
Ermine1_standardized = data.frame(Ermine1_scaled$Center_win, Ermine1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Hyde1_standardized, type="l", xlim = c(400, 2150), ylim = c(-20, 20), xlab = "Hyde1 Resampled Depth", ylab = "Normalized GR")
plot(Ermine1_standardized, type="l", xlim = c(600, 2000), ylim = c(-20, 30), xlab = "Ermine1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_e1_h1_ap1 <- dtw(Ermine1_standardized$Ermine1_scaled.Average, Hyde1_standardized$Hyde1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = F, open.end = F))
plot(al_e1_h1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Ermine1_on_Hyde1_depth = tune(Ermine1_standardized, cbind(Ermine1_standardized$Ermine1_scaled.Center_win[al_e1_h1_ap1$index1s], Hyde1_standardized$Hyde1_scaled.Center_win[al_e1_h1_ap1$index2s]), extrapolate = T)

dev.off()

# Plotting the data

plot(Hyde1_standardized, type = "l", ylim = c(-20, 20), xlim = c(400, 2100), xlab = "Hyde1 Resampled Depth", ylab = "Normalized GR")
lines(Ermine1_on_Hyde1_depth, col = "red")

# DTW Distance
al_e1_h1_ap1$normalizedDistance
al_e1_h1_ap1$distance

# Tuning Ermine1 data on Andromeda1 depth
Ermine1_on_Andromeda1_depth = tune(Ermine1_on_Hyde1_depth, cbind(Hyde1_standardized$Hyde1_scaled.Center_win[al_h1_am1_ap1$index1s], Andromeda1_standardized$Andromeda1_scaled.Center_win[al_h1_am1_ap1$index2s]), extrapolate = F)

# Tuning Ermine1 data on Fisher1 depth
Ermine1_on_Fisher1_depth = tune(Ermine1_on_Andromeda1_depth, cbind(Andromeda1_standardized$Andromeda1_scaled.Center_win[al_am1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_am1_fi1_ap1$index2s]), extrapolate = F)

# Tuning Ermine1 data on Finucane1 depth
Ermine1_on_Finucane1_depth = tune(Ermine1_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning Ermine1 data on Picard1 depth
Ermine1_on_Picard1_depth = tune(Ermine1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Ermine-1)")
lines(Ermine1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Ermine1_originalGR_on_Picard1_depth = data.frame(Ermine1_on_Picard1_depth$X1, Ermine1_interpolated[1:6850,2])

plot(Picard1_originalGR, type = "l", ylim = c(0, 70), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Ermine-1)")
lines(Ermine1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Ermine1 

U1463Age_on_Ermine1_depth = tune(Ermine1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Ermine1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Ermine1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Ermine1_depth) <- new_column_names
write.csv(U1463Age_on_Ermine1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Ermine 1.csv", row.names = FALSE)

Ermine1_age_depth <- approx(x = Ermine1_interpolated$GR,
                               y = Ermine1_interpolated$DEPT,
                               xout = U1463Age_on_Ermine1_depth$GR,
                               rule = 1)$y

Ermine1_agemodel = data.frame(Ermine1_age_depth, U1463Age_on_Ermine1_depth$AGE)
Ermine1_agemodel <- na.omit(Ermine1_agemodel)
new_column_names1 <- c("Depth", "Age")
colnames(Ermine1_agemodel) <- new_column_names1
plot(Ermine1_agemodel)
write.csv(Ermine1_agemodel, file = "RScripts&Data/Sites Data_Age-NGR/Ermine1_DepthAge.csv", row.names = FALSE)
