install.packages(setdiff(c("DescTools", "astrochron", "dtw"), rownames(installed.packages())))

# Import packages

library(dtw)
library(DescTools)
library(astrochron)

# Import Brontosaurus 1 and SouthGalapagos1 datasets

Brontosaurus11 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/Brontosaurus 1.csv", header=TRUE, stringsAsFactors=FALSE)
Brontosaurus11=Brontosaurus11[c(1241:4527),] # Oligocene-Miocene
head(Brontosaurus11)
plot(Brontosaurus11, type="l", xlim = c(400, 1050), ylim = c(20, 60))

SouthGalapagos1 <- read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1.csv", header=TRUE, stringsAsFactors=FALSE)
SouthGalapagos1=SouthGalapagos1[c(30:3710),] # Oligocene-Miocene
head(SouthGalapagos1)
plot(SouthGalapagos1, type="l", xlim = c(500, 1200), ylim = c(0, 70))

#### Rescaling and resampling of the data ####

# Linear interpolation of datasets
Brontosaurus11_interpolated <- linterp(Brontosaurus11, dt = 0.2, genplot = F)
SouthGalapagos1_interpolated <- linterp(SouthGalapagos1, dt = 0.2, genplot = F)

# Scaling the data
Bmean = Gmean(Brontosaurus11_interpolated$GR)
Bstd = Gsd(Brontosaurus11_interpolated$GR)
Brontosaurus11_scaled = (Brontosaurus11_interpolated$GR - Bmean)/Bstd
Brontosaurus11_rescaled = data.frame(Brontosaurus11_interpolated$DEPT, Brontosaurus11_scaled)

Smean = Gmean(SouthGalapagos1_interpolated$GR)
Sstd = Gsd(SouthGalapagos1_interpolated$GR)
SouthGalapagos1_scaled = (SouthGalapagos1_interpolated$GR - Smean)/Sstd
SouthGalapagos1_rescaled = data.frame(SouthGalapagos1_interpolated$DEPT, SouthGalapagos1_scaled)

# Resampling the data using moving window statistics
Brontosaurus11_scaled = mwStats(Brontosaurus11_rescaled, cols = 2, win=3, ends = T)
Brontosaurus11_standardized = data.frame(Brontosaurus11_scaled$Center_win, Brontosaurus11_scaled$Average)

SouthGalapagos1_scaled = mwStats(SouthGalapagos1_rescaled, cols = 2, win=3, ends = T)
SouthGalapagos1_standardized = data.frame(SouthGalapagos1_scaled$Center_win, SouthGalapagos1_scaled$Average)

# Plotting the rescaled and resampled data
plot(Brontosaurus11_standardized, type="l", xlim = c(400, 1100), ylim = c(-20, 20), xlab = "Brontosaurus11 Resampled Depth", ylab = "Normalized GR")
plot(SouthGalapagos1_standardized, type="l", xlim = c(500, 1200), ylim = c(-20, 20), xlab = "SouthGalapagos1 Resampled Depth", ylab = "Normalized GR")

#### DTW with custom step pattern asymmetricP1.1 but no custom window ####

# Perform dtw
system.time(al_sg1_b1_ap1 <- dtw(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Average, Brontosaurus11_standardized$Brontosaurus11_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, open.begin = T, open.end = F))
plot(al_sg1_b1_ap1, "threeway")

# Tuning the standardized data on reference depth scale
SouthGalapagos1_on_Brontosaurus11_depth = tune(SouthGalapagos1_standardized, cbind(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_sg1_b1_ap1$index1s], Brontosaurus11_standardized$Brontosaurus11_scaled.Center_win[al_sg1_b1_ap1$index2s]), extrapolate = F)

dev.off()

# Plotting the data

plot(Brontosaurus11_standardized, type = "l", ylim = c(-20, 20), xlim = c(400, 1050), xlab = "Brontosaurus11 Resampled Depth", ylab = "Normalized GR")
lines(SouthGalapagos1_on_Brontosaurus11_depth, col = "red")

# DTW Distance
al_sg1_b1_ap1$normalizedDistance
al_sg1_b1_ap1$distance


#### DTW with custom step pattern asymmetricP1.1 and custom window ####

# create matrix for the custom window

compare.window <- matrix(data=TRUE,nrow=nrow(SouthGalapagos1_standardized),ncol=nrow(Brontosaurus11_standardized))
image(x=Brontosaurus11_standardized[,1],y=SouthGalapagos1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Assigning stratigraphic depth locations for reference and query sites

# Depth values for first datum
base_1_x <- which.min(abs(Brontosaurus11_standardized[,1] - 450))
base_1_y <- which.min(abs(SouthGalapagos1_standardized[,1] - 525))

# Depth values for second datum
base_2_x <- which.min(abs(Brontosaurus11_standardized[,1] - 620))
base_2_y <- which.min(abs(SouthGalapagos1_standardized[,1] - 655))

# Assigning depth uncertainty "slack" to the tie-points

# Create a matrix to store the comparison window
compare.window <- matrix(data = TRUE, nrow = nrow(SouthGalapagos1_standardized), ncol = nrow(Brontosaurus11_standardized))

# Slack provided based on specific indices 

compare.window[(base_1_y+5):nrow(SouthGalapagos1_standardized),1:(base_1_x-5)] <- 0
compare.window[1:(base_1_y-5),(base_1_x+5):ncol(compare.window)] <- 0

compare.window[(base_2_y+5):nrow(SouthGalapagos1_standardized),1:(base_2_x-5)] <- 0
compare.window[1:(base_2_y-5),(base_2_x+5):ncol(compare.window)] <- 0

# Visualize the comparison window
image(x=Brontosaurus11_standardized[,1],y=SouthGalapagos1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Convert the comparison window matrix to logical values
compare.window <- sapply(as.data.frame(compare.window), as.logical)
compare.window <- unname(as.matrix(compare.window))

image(x=Brontosaurus11_standardized[,1],y=SouthGalapagos1_standardized[,1],z=t(compare.window),useRaster=TRUE)

# Define a custom window function for use in DTW
win.f <- function(iw,jw,query.size, reference.size, window.size, ...) compare.window >0

# Perform dtw with custom window
system.time(al_sg1_b1_ap2 <- dtw(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Average, Brontosaurus11_standardized$Brontosaurus11_scaled.Average, keep.internals = T, step.pattern = asymmetricP1.1, window.type = win.f, open.end = F, open.begin = F))
plot(al_sg1_b1_ap2, type = "threeway")

# DTW Distance measure
al_sg1_b1_ap2$normalizedDistance
al_sg1_b1_ap2$distance

image(y = Brontosaurus11_standardized[,1], x = SouthGalapagos1_standardized[,1], z = compare.window, useRaster = T)
lines(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_sg1_b1_ap2$index1], Brontosaurus11_standardized$Brontosaurus11_scaled.Center_win[al_sg1_b1_ap2$index2], col = "white", lwd = 2)

# Tuning the standardized data on reference depth scale
SouthGalapagos1_on_Brontosaurus11_depth1 = tune(SouthGalapagos1_standardized, cbind(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_sg1_b1_ap2$index1s], Brontosaurus11_standardized$Brontosaurus11_scaled.Center_win[al_sg1_b1_ap2$index2s]), extrapolate = F)

dev.off()

# Plotting the data
plot(Brontosaurus11_standardized, type = "l", ylim = c(-20, 20), xlim = c(400, 1050), xlab = "Brontosaurus-1 Resampled Depth", ylab = "Normalized GR (South Galapagos-1)")
lines(SouthGalapagos1_on_Brontosaurus11_depth1, col = "red")

# Tuning Brontosaurus11 data on SouthGalapagos1 depth
Brontosaurus11_on_SouthGalapagos1_depth = tune(Brontosaurus11_standardized, cbind(Brontosaurus11_standardized$Brontosaurus11_scaled.Center_win[al_sg1_b1_ap2$index2s], SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_sg1_b1_ap2$index1s]), extrapolate = F)

# Changing the GR values to original and reploting

SouthGalapagos1_originalGR = data.frame(SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win, SouthGalapagos1_interpolated$GR)
Brontosaurus11_originalGR_on_SG1_depth = data.frame(Brontosaurus11_on_SouthGalapagos1_depth$X1, Brontosaurus11_interpolated[1:3287,2])

plot(SouthGalapagos1_originalGR, type = "l", ylim = c(0, 80), xlim = c(500, 1200), xlab = "South Galapagos-1 Depth", ylab = "GR (Brontosaurus-1)")
lines(Brontosaurus11_originalGR_on_SG1_depth, col = "red")

# Age Model
AgeModelSouthGalapagos <-read.csv("RScripts&Data/Sites Data_Depth-NGR/SouthGalapagos1_DepthAge.csv", header=TRUE, stringsAsFactors=FALSE)
AgeModelSouthGalapagos = data.frame(AgeModelSouthGalapagos$Depth, AgeModelSouthGalapagos$Time_Ma)
plot(AgeModelSouthGalapagos, type="l")

# Tuning the age model data to Brontosaurus11 

SG1Age_on_Brontosaurus11_depth = tune(Brontosaurus11_originalGR_on_SG1_depth, AgeModelSouthGalapagos, extrapolate = F)
dev.off()

plot(SG1Age_on_Brontosaurus11_depth, type = "l", ylim = c(20, 60), xlim = c(2.5, 22), xaxt = "n", xlab = "Age (Ma)", ylab = "Brontosaurus11")
axis(1, at = c(2.5,5,10,15,20), cex.axis = 1.0, las = 1)

new_column_names <- c("AGE", "GR")
colnames(SG1Age_on_Brontosaurus11_depth) <- new_column_names
write.csv(SG1Age_on_Brontosaurus11_depth, file = "RScripts&Data/Sites Data_Age-NGR/Brontosaurus 1.csv", row.names = FALSE)

SG1_age_aligned <- approx(x = AgeModelSouthGalapagos$AgeModelSouthGalapagos.Depth,
                          y = AgeModelSouthGalapagos$AgeModelSouthGalapagos.Time_Ma,
                          xout = SouthGalapagos1_standardized$SouthGalapagos1_scaled.Center_win[al_sg1_b1_ap2$index1s],
                          rule = 1)$y

Bronto1_agemodel = data.frame(SG1_age_aligned, Brontosaurus11_standardized$Brontosaurus11_scaled.Center_win[al_sg1_b1_ap2$index2s])
Bronto1_agemodel <- na.omit(Bronto1_agemodel)
new_column_names1 <- c("SG1_Age", "Brontosaurus11_Depth")
colnames(Bronto1_agemodel) <- new_column_names1

Bronto1_SF = data.frame(SG1_Age = 0, Brontosaurus11_Depth = 152.2)
Bronto1_agemodel <- rbind(Bronto1_SF, Bronto1_agemodel)
Bronto1_sedrate=Bronto1_agemodel[c(1:2912),]
for (i in c(1:2911)){
  Bronto1_sedrate[i,2]=(Bronto1_agemodel[i+1,2]-Bronto1_agemodel[i,2])*100/(Bronto1_agemodel[i+1,1]-Bronto1_agemodel[i,1])
} 
new_column_names2 <- c("Age", "Sedimentation Rate")
colnames(Bronto1_sedrate) <- new_column_names2


Bronto1_sedrate = mwStats(Bronto1_sedrate, cols = 2, win=1, ends = T)
Bronto1_sedrate = data.frame(Bronto1_sedrate$Center_win, Bronto1_sedrate$Average)

plot(Bronto1_sedrate, lwd = 1.2, col = "black", type = "l", xlim = c(2.54, 21.52), ylim = c(0,22000))

SG1_sedrate = data.frame(AgeModelSouthGalapagos$AgeModelSouthGalapagos.Time_Ma, AgeModelSouthGalapagos$AgeModelSouthGalapagos.Depth)
for (i in c(1:7215)){
  SG1_sedrate[i,2]=(SG1_sedrate[i+1,2]-SG1_sedrate[i,2])*100/(SG1_sedrate[i+1,1]-SG1_sedrate[i,1])
} 
plot(SG1_sedrate, lwd = 1.2, col = "black", type = "l", xlim = c(2.54, 21.52), ylim = c(2000,12000))

