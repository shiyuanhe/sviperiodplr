## Plot some mira examples for Simulation 1

# ------- simu1Plot ---------
load("./result/simu1/SVIresult.RData")
source("./params.R")
# read in the true frequency
trueFreq = read.table(paste0(PARAMS_DATAPATH,"/id_map_O.dat"), 
                      header = TRUE, stringsAsFactors = FALSE)
allQfFiles = list.files("./result/simu1/SVIqf", full.names = TRUE)

diff = abs(trueFreq$period -  allPL[1,]) / trueFreq$period
plotQf = function(starI, Phat, trueP){
    qf = read.table(allQfFiles[starI], header = FALSE)
    plot(qf, type = "l", xlab = "", ylab = "")
    points(1 / Phat, max(qf[,2])*0.05, pch = 25, col = "red", lwd = 2)
    points(1 / trueP, -max(qf[,2])*0.02, pch = 24, col = "blue", lwd = 2)
    #abline(v = 1 / allPL[1, starI], lwd = 2, lty = 2, col = "red")
    #abline(v = 1 / trueP, col = "blue", lty = 3, lwd = 2)
    mtext(expression(paste("Frequency [", day^-1, "]")),
          side = 1, line = 2, cex = 1)
    mtext("q(f)", side = 2, line = 2, cex = 1.1)
}

trueP = trueFreq$period[c(1,2,221,720)]


if(SAVEPLOT) pdf("./figures/simu1-1q.pdf", height = 3.5, width = 5)
par(mar = c(3,3,1,1) )
starI = 1
plotQf(starI, allPL[1, starI], trueFreq$period[starI])
legend(0.006,40000, 
       legend = c("True freq", "Estimated freq"),
       pch = c(24, 25), lwd = 2, lty = 0,
       col = c("blue", "red"), cex = 0.9)
if(SAVEPLOT) dev.off()

if(SAVEPLOT) pdf("./figures/simu1-2q.pdf", height = 3.5, width = 5)
par(mar = c(3,3,1,1) )
starI = 2
plotQf(starI, allPL[1, starI], trueFreq$period[starI])
if(SAVEPLOT) dev.off()

if(SAVEPLOT) pdf("./figures/simu1-1084q.pdf", height = 3.5, width = 5)
par(mar = c(3,3,1,1) )
starI = 1084
plotQf(starI, allPL[1, starI], trueFreq$period[starI])
if(SAVEPLOT) dev.off()

if(SAVEPLOT) pdf("./figures/simu1-1771q.pdf", height = 3.5, width = 5)
par(mar = c(3,3,1,1) )
starI = 1771
plotQf(starI, allPL[1, starI], trueFreq$period[starI])
if(SAVEPLOT) dev.off()

## Check convergance for gamma
## Each column is the parameter for each round of iteration. 
# x = modelObj$get_gammaTrace()
# plot(x[4,], pch = 20, type = "l")
# plot(x[3,], pch = 20, type = "l")
# plot(x[2,], pch = 20, type = "l")
# plot(x[1,], pch = 20, type = "l")
# 
# # Check convergence for alpha
# y = modelObj$get_alphaTrace()
# plot(y[4,], type = "l")
# plot(y[5,], type = "l")
# plot(y[6,], type = "l")