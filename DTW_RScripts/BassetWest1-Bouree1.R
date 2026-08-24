install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import BassetWest 1 and Bouree 1 datasets

BassetWest1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/BassetWest1.csv", header=TRUE, stringsAsFactors=FALSE)
BassetWest1=BassetWest1[c(1:3384),]
head(BassetWest1)
plot(BassetWest1, type='l', xlim= c(450,1150), ylim = c(0,90))

Bouree1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Bouree 1.csv", header=TRUE, stringsAsFactors=FALSE)
Bouree1=Bouree1[c(1:2689),]
head(Bouree1)
plot(Bouree1, type='l', xlim= c(350,950), ylim = c(0,90))


#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
BassetWest1_interpolated <- linterp(BassetWest1, dt = 0.2, genplot = F)
Bouree1_interpolated <- linterp(Bouree1, dt = 0.2, genplot = F)

# Scaling the data
Bmean = Gmean(BassetWest1_interpolated$GR)
Bstd = Gsd(BassetWest1_interpolated$GR)
BassetWest1_scaled = (BassetWest1_interpolated$GR - Bmean)/Bstd
BassetWest1_rescaled = data.frame(BassetWest1_interpolated$DEPT, BassetWest1_scaled)

Bomean = Gmean(Bouree1_interpolated$GR)
Bostd = Gsd(Bouree1_interpolated$GR)
Bouree1_scaled = (Bouree1_interpolated$GR - Bomean)/Bostd
Bouree1_rescaled = data.frame(Bouree1_interpolated$DEPT, Bouree1_scaled)

# Resampling the data using moving window statistics
BassetWest1_scaled = mwStats(BassetWest1_rescaled, cols = 2, win=3, ends = T)
BassetWest1_standardized = data.frame(BassetWest1_scaled$Center_win, BassetWest1_scaled$Average)

Bouree1_scaled = mwStats(Bouree1_rescaled, cols = 2, win=3, ends = T)
Bouree1_standardized = data.frame(Bouree1_scaled$Center_win, Bouree1_scaled$Average)

# Plotting the rescaled and resampled data
plot(BassetWest1_standardized, type="l", xlim = c(450, 1150), ylim = c(-20, 30), xlab = "BassetWest 1 Resampled Depth", ylab = "Normalized GR")
plot(Bouree1_standardized, type="l", xlim = c(350, 950), ylim = c(-20, 20), xlab = "Bouree 1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_b1_b1_ap1 <- dtw(Bouree1_standardized$Bouree1_scaled.Average, BassetWest1_standardized$BassetWest1_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = T))
plot(al_b1_b1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
Bouree1_on_BassetWest1_depth = tune(Bouree1_standardized, cbind(Bouree1_standardized$Bouree1_scaled.Center_win[al_b1_b1_ap1$index1s], BassetWest1_standardized$BassetWest1_scaled.Center_win[al_b1_b1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(BassetWest1_standardized, type = "l", ylim = c(-20, 30), xlim = c(450, 1150), xlab = "BassetWest1 Resampled Depth", ylab = "Normalized GR")
lines(Bouree1_on_BassetWest1_depth, col = "red")

# DTW Distance
al_b1_b1_ap1$normalizedDistance
al_b1_b1_ap1$distance

# Tuning BassetWest1 data on Bouree1 depth
BassetWest1_on_Bouree1_depth = tune(BassetWest1_standardized, cbind(BassetWest1_standardized$BassetWest1_scaled.Center_win[al_b1_b1_ap1$index2s], Bouree1_standardized$Bouree1_scaled.Center_win[al_b1_b1_ap1$index1s]), extrapolate = F)

# Tuning BassetWest1 data on Gorgonichthys1 depth
BassetWest1_on_Gorgonichthys1_depth = tune(BassetWest1_on_Bouree1_depth, cbind(Bouree1_standardized$Bouree1_scaled.Center_win[al_b1_g1_ap1$index1s], Gorgonichthys1_standardized$Gorgonichthys1_scaled.Center_win[al_b1_g1_ap1$index2s]), extrapolate = F)

# Tuning BassetWest1 data on Walkley1 depth
BassetWest1_on_Walkley1_depth = tune(BassetWest1_on_Gorgonichthys1_depth, cbind(Gorgonichthys1_standardized$Gorgonichthys1_scaled.Center_win[al_g1_w1_ap2$index1s], Walkley1_standardized$Walkley1_scaled.Center_win[al_g1_w1_ap2$index2s]), extrapolate = F)

# Tuning BassetWest1 data on Caswell1 depth
BassetWest1_on_Caswell1_depth = tune(BassetWest1_on_Walkley1_depth, cbind(Walkley1_standardized$Walkley1_scaled.Center_win[al_w1_c1_ap1$index1s], Caswell1_standardized$Caswell1_scaled.Center_win[al_w1_c1_ap1$index2s]), extrapolate = F)

# Tuning BassetWest1 data on Calliance2 depth
BassetWest1_on_Calliance2_depth = tune(BassetWest1_on_Caswell1_depth, cbind(Caswell1_standardized$Caswell1_scaled.Center_win[al_c2_c1_ap1$index2s], Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_c1_ap1$index1s]), extrapolate = F)

# Tuning BassetWest1 data on Omar1 depth
BassetWest1_on_Omar1_depth = tune(BassetWest1_on_Calliance2_depth, cbind(Calliance2_standardized$Calliance2_scaled.Center_win[al_c2_o1_ap2$index1s], Omar1_standardized$Omar1_scaled.Center_win[al_c2_o1_ap2$index2s]), extrapolate = F)

# Tuning BassetWest1 data on SG1 depth
BassetWest1_on_SG1_depth = tune(BassetWest1_on_Omar1_depth, cbind(Omar1_standardized$Omar1_scaled.Center_win[al_o1_sg1_ap2$index1s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_o1_sg1_ap2$index2s]), extrapolate = F)

dev.off()
plot(SouthGalapagos1_standardized, type = "l", ylim = c(-20, 20), xlim = c(500, 1200), xlab = "SG1 Resampled Depth", ylab = "Normalized GR (BassetWest-1)")
lines(BassetWest1_on_SG1_depth, col = "red")

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
BassetWest1_originalGR_on_SG1_depth = data.frame(BassetWest1_on_SG1_depth$X1, BassetWest1_interpolated[567:3145,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (BassetWest-1)")
lines(BassetWest1_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to BassetWest1 

SG1Age_on_BassetWest1_depth = tune(BassetWest1_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_BassetWest1_depth, type = "l", ylim = c(0, 90), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "BassetWest1")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_BassetWest1_depth) <- new_column_names
SG1Age_on_BassetWest1_depth[,1] = SG1Age_on_BassetWest1_depth[,1] * 1000
write.csv(SG1Age_on_BassetWest1_depth, file = "RScripts&Data/Sites Data_Age-NGR/BassetWest 1.csv", row.names = FALSE)
