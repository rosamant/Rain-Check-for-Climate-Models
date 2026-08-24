rm(list=ls())

library(astrochron)

# Required color scales

CENOGRID = read.csv("RScripts&Data/CENOGRID_revised2025.csv")

CENO1 = data.frame(CENOGRID[,1],CENOGRID[,3])
CENO1 = CENO1[c(2458:11756),]
CENO1 = linterp(CENO1, dt = 0.02, start = 2.49)
dev.off()

plot(CENO1, type = "l")



setwd("RScripts&Data/Time_Slice_Map")

png(filename = "CENO.png", width = 2500, height = 6000, res = 600)

par(mar=c(5,5,3,1))

plot(CENO1$CENOGRID...3.,CENO1$CENOGRID...1., type = "l", xlim = c(5.0,1.0), ylim = c(22.0,2.0), xaxs = "i", yaxs = "i", col = "grey40", lwd = 2, axes = F, xaxt = "n", xlab = "", ylab = "")
axis(1, at = c(4.5, seq(1.5,4.5,1)), cex.axis = 1.6, col = "grey40", col.axis = "grey40")
axis(2, at = c(seq(2,22,2)), cex.axis = 1.6)
mtext(expression("CENOGRID"~delta^18*"O (‰)"), side = 1, line = 3, cex = 1.5, col = "grey40")
mtext("Age (Ma)", side = 2, line = 2.5, cex = 1.5)

dev.off()
