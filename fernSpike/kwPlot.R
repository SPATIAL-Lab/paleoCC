
# Biosphere size
B = seq(0, 2, by = 0.01)

# Biosphere factor
Bf = (B - 0.7) / 0.3

# Kerogen weathering
kw = 4 - ((sign(Bf) * abs(Bf)^(1/3)) + 2)

png("fernSpike/figs/kwFeedback.png", width = 4, height = 4, units = "in", 
    res = 600)
par(mar = c(5, 5, 1, 1))
plot(B, kw, type = "l", lwd = 3, xlab = expression("B"[t]*"/B"[0]),
     ylab = expression("kw"[t]*"/kw"[0]))
dev.off()
