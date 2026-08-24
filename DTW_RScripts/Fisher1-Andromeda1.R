install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtwMultiAlign)
library(dtw)
library(DescTools)
library(astrochron)

# Import Fisher1 and Andromeda1 datasets

Fisher1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Fisher1.csv", header=TRUE, stringsAsFactors=FALSE)
Fisher1=Fisher1[c(1:10190),] # Oligocene-Miocene
head(Fisher1)
plot(Fisher1, type="l", xlim = c(100, 1700), ylim = c(0, 50))

Andromeda1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Andromeda1.csv", header=TRUE, stringsAsFactors=FALSE)

# Recorrecting signal
An1 = Gmean(Andromeda1[c(5997:9740),2])
An2 = Gmean(Andromeda1[c(3000:5000),2])
SD1 = Gsd(Andromeda1[c(5997:9740),2])
SD2 = Gsd(Andromeda1[c(3000:7000),2])
Andromeda1[c(5997:9740),2]=(Andromeda1[c(5997:9740),2]+(An2-An1))*(SD1/SD2)

# Recorrecting signal
An1 = Gmean(Andromeda1[c(7099:9740),2])
An2 = Gmean(Andromeda1[c(3000:7000),2])
SD1 = Gsd(Andromeda1[c(7099:9740),2])
SD2 = Gsd(Andromeda1[c(3000:7000),2])
Andromeda1[c(7099:9740),2]=(Andromeda1[c(7099:9740),2]+(An2-An1))*(SD1/SD2)

Andromeda1=Andromeda1[c(1:9573),] # Oligocene-Miocene
head(Andromeda1)
plot(Andromeda1, type="l", xlim = c(400, 2200), ylim = c(0, 60))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Fisher1_interpolated <- linterp(Fisher1, dt = 0.2, genplot = F)
Andromeda1_interpolated <- linterp(Andromeda1, dt = 0.2, genplot = F)

# Scaling the data
Fimean = Gmean(Fisher1_interpolated$GR)
Fistd = Gsd(Fisher1_interpolated$GR)
Fisher1_scaled = (Fisher1_interpolated$GR - Fimean)/Fistd
Fisher1_rescaled = data.frame(Fisher1_interpolated$DEPT, Fisher1_scaled)

Ammean = Gmean(Andromeda1_interpolated$GR)
Amstd = Gsd(Andromeda1_interpolated$GR)
Andromeda1_scaled = (Andromeda1_interpolated$GR - Ammean)/Amstd
Andromeda1_rescaled = data.frame(Andromeda1_interpolated$DEPT, Andromeda1_scaled)

# Resampling the data using moving window statistics
Fisher1_scaled = mwStats(Fisher1_rescaled, cols = 2, win=3, ends = T)
Fisher1_standardized = data.frame(Fisher1_scaled$Center_win, Fisher1_scaled$Average)

Andromeda1_scaled = mwStats(Andromeda1_rescaled, cols = 2, win=3, ends = T)
Andromeda1_standardized = data.frame(Andromeda1_scaled$Center_win, Andromeda1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Fisher1_standardized, type="l", xlim = c(100, 1700), ylim = c(-20, 20), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
plot(Andromeda1_standardized, type="l", xlim = c(350, 2200), ylim = c(-20, 20), xlab = "Andromeda1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_am1_fi1_ap1 <- dtw(Andromeda1_standardized$Andromeda1_scaled.Average, Fisher1_standardized$Fisher1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_am1_fi1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Andromeda1_on_Fisher1_depth = tune(Andromeda1_standardized, cbind(Andromeda1_standardized$Andromeda1_scaled.Center_win[al_am1_fi1_ap1$index1s], Fisher1_standardized$Fisher1_scaled.Center_win[al_am1_fi1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Fisher1_standardized, type = "l", ylim = c(-20, 20), xlim = c(100, 1700), xlab = "Fisher1 Resampled Depth", ylab = "Normalized GR")
lines(Andromeda1_on_Fisher1_depth, col = "red")

# DTW Distance and RMSE

al_am1_fi1_ap1$normalizedDistance
al_am1_fi1_ap1$distance

# Tuning Andromeda1 data on Finucane1 depth
Andromeda1_on_Finucane1_depth = tune(Andromeda1_on_Fisher1_depth, cbind(Fisher1_standardized$Fisher1_scaled.Center_win[al_fi1_f1_ap1$index1s], Finucane1_standardized$Finucane1_scaled.Center_win[al_fi1_f1_ap1$index2s]), extrapolate = F)

# Tuning Andromeda1 data on Picard1 depth
Andromeda1_on_Picard1_depth = tune(Andromeda1_on_Finucane1_depth, cbind(Finucane1_standardized$Finucane1_scaled.Center_win[al_f1_p1_ap1$index1s], Picard1_standardized$Picard1_scaled.Center_win[al_f1_p1_ap1$index2s]), extrapolate = F)

dev.off()
plot(Picard1_standardized, type = "l", ylim = c(-20, 20), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Andromeda-1)")
lines(Andromeda1_on_Picard1_depth, col = "red")

# Changing the GR values to original and reploting

Picard1_originalGR = data.frame(Picard1_standardized$Picard1_scaled.Center_win, Picard1_interpolated$GR)
Andromeda1_originalGR_on_Picard1_depth = data.frame(Andromeda1_on_Picard1_depth$X1, Andromeda1_interpolated[1:8906,2])

plot(Picard1_originalGR, type = "l", ylim = c(0, 50), xlim = c(150, 1300), xlab = "Picard1 Resampled Depth", ylab = "Normalized GR (Andromeda-1)")
lines(Andromeda1_originalGR_on_Picard1_depth, col = "red")

# Age Model
AgeModelPicard <-read.csv("RScripts&Data/Sites Data_Depth-NGR/Picard1-U1463_AgeModel.csv", header=TRUE, stringsAsFactors=FALSE)
plot(AgeModelPicard, type="l")

# Tuning the age model data to Andromeda1 

U1463Age_on_Andromeda1_depth = tune(Andromeda1_originalGR_on_Picard1_depth, AgeModelPicard, extrapolate = F)
dev.off()

plot(U1463Age_on_Andromeda1_depth, type = "l", ylim = c(0, 50), xlim = c(500, 21000), xaxt = "n", xlab = "Age (ka)", ylab = "Andromeda1")
axis(1, at = c(440,5000,10000,15000,20000), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(U1463Age_on_Andromeda1_depth) <- new_column_names
write.csv(U1463Age_on_Andromeda1_depth, file = "RScripts&Data/Sites Data_Age-NGR/Andromeda 1.csv", row.names = FALSE)
