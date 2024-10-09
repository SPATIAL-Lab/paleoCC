source("helpers.R")

d = read.data()
plot(d$Time - 20000, d$PCO2, type = "l")
plot(d$Time - 20000, d$d13CAtm, type = "l")
plot(d$Time - 20000, d$d13CSurf, type = "l")
plot(d$Time - 20000, d$BioC, type = "l")
plot(d$Time - 20000, d$OrgBurial, type = "l")

fb = read.data()
base = read.data()

png("fernSpike/figs/d13C.png", width = 6, height = 8, units = "in", res = 600)
layout(matrix(1:3, nrow = 3))
par(mar = c(5, 5, 1, 5))
plot(fb$Time - 20000, fb$d13CAtm, type = "l", lwd = 3, axes = FALSE, 
     xlab = "", ylab = "")
lines(base$Time - 20000, base$d13CAtm)
axis(2)
mtext(expression(delta^{13}*"C"[atm]), 2, 3)
plot(fb$Time - 20000, fb$d13CSurf, type = "l", lwd = 3, axes = FALSE, 
     xlab = "", ylab = "")
lines(base$Time - 20000, base$d13CSurf)
axis(4)
mtext(expression(delta^{13}*"C"[surf]), 4, 3)
plot(fb$Time - 20000, fb$d13CDeep, type = "l", lwd = 3, axes = FALSE, 
     xlab = "", ylab = "")
lines(base$Time - 20000, base$d13CDeep)
axis(2)
axis(1)
mtext("Time (kyr)", 1, 3)
mtext(expression(delta^{13}*"C"[deep]), 2, 3)
dev.off()


png("fernSpike/figs/CC.png", width = 6, height = 8, units = "in", res = 600)
layout(matrix(1:3, nrow = 3))
par(mar = c(5, 5, 1, 5))
plot(fb$Time - 20000, fb$Inject * 12e3, type = "l", lwd = 3, axes = FALSE, 
     xlab = "", ylab = "")
lines(base$Time - 20000, base$Inject * 12e3)
axis(2)
mtext("C injection (Pg/yr)", 2, 3)
plot(fb$Time - 20000, fb$PCO2, type = "l", lwd = 3, axes = FALSE, 
     xlab = "", ylab = "")
lines(base$Time - 20000, base$PCO2)
axis(4)
mtext(expression("pCO"[2]*" (ppmv)"), 4, 3)
plot(fb$Time - 20000, fb$BioC * 12e3, type = "l", lwd = 3, axes = FALSE, 
     xlab = "", ylab = "")
lines(base$Time - 20000, base$BioC * 12e3)
axis(2)
axis(1)
mtext("Time (kyr)", 1, 3)
mtext("Terrestrial biosphere (Pg)", 2, 3)
dev.off()



d = read.data()
nw = d

d = read.data()
fpoc = d

d = read.data()
ow = d

d = read.data()
allw = d

plot(allw$Time - 20000, allw$PCO2, type = "l")
lines(nw$Time - 20000, nw$PCO2, col = "blue")
lines(ow$Time - 20000, ow$PCO2, col = "brown")
lines(fpoc$Time - 20000, fpoc$PCO2, col = "darkgreen")

plot(allw$Time - 20000, allw$BioC, type = "l")
lines(nw$Time - 20000, nw$BioC, col = "blue")
lines(fpoc$Time - 20000, fpoc$BioC, col = "darkgreen")
