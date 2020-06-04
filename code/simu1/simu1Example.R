## Draw four light curves for the first M33 simulation

# ------- simu1Example ---------

source("./params.R")

dataPath = paste0(PARAMS_DATAPATH,"simu1/")
allFiles = list.files(dataPath, full.names = TRUE)
nSample = length(allFiles)


library(readr)
library(tidyverse)

drawOneSimu1 = function(ii){
    # Data for the simulated dataset
    starObs = read.table(allFiles[ii], stringsAsFactors = FALSE)
    colors = c("black", "red", "green", "blue")
    shape = c(20, 15, 17, 18)
    starObs$V2 = factor(starObs$V2, levels = c("I","J","H","K"))
    with(starObs, plot(V1, V3,  ylim = rev(range(V3)), type = "n", 
                       xlab = "", ylab = ""))
    mtext("Date [day]", side = 1, line = 2, cex = 0.9)
    mtext("Magnitude", side = 2, line = 2, cex = 0.9)
    with(starObs, arrows(V1, V3 + V4, V1, V3 - V4, code = 3, length = 0.02,  angle = 90, col = "grey"))
    with(starObs, points(V1, V3, pch = shape[as.numeric(V2)], col = colors[as.numeric(V2)]))
    if(ii ==1){
        legend(1000, 17, legend = c(expression(italic(I)),
                                    expression(italic(J)),
                                    expression(italic(H)),
                                    expression(italic(K)[s])), 
               pch = shape, col = colors)
    }
}

if(SAVEPLOT) pdf("./figures/simu1-1.pdf", height = 3.5, width = 5)
par(mar = c(3,3,1,1) )
drawOneSimu1(1)
if(SAVEPLOT) dev.off()

if(SAVEPLOT) pdf("./figures/simu1-2.pdf", height = 3.5, width = 5)
par(mar = c(3,3,1,1) )
drawOneSimu1(2)
if(SAVEPLOT) dev.off()

if(SAVEPLOT) pdf("./figures/simu1-1084.pdf", height = 3.5, width = 5)
par(mar = c(3,3,1,1) )
drawOneSimu1(1084)
if(SAVEPLOT) dev.off()

if(SAVEPLOT) pdf("./figures/simu1-1771.pdf", height = 3.5, width = 5)
par(mar = c(3,3,1,1) )
drawOneSimu1(1771)
if(SAVEPLOT) dev.off()

