source("helpers.R")

d = read.data()
plot(d$Time, d$Inject, type = "l")
plot(d$Time, d$PCO2, type = "l")
plot(d$Time, d$BioC, type = "l")
plot(d$Time, d$Stemp, type = "l")

plot(d$Time, exp(-d$Inject*2e4), type = "l")
