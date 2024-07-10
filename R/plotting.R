exo = read.table("exogenic.txt", header = TRUE, fill = TRUE)
exo = head(exo, -2)

plot(exo$Time, exo$PCO2, type = "l")
plot(exo$Time, exo$d13CBio, type = "l")
plot(exo$Time, exo$d13CSurf, type = "l")
