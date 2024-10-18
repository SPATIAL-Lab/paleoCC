write.config = function(parms){
  str = paste0(parms$inj, "\n",
               parms$fb_bio, "\n",
               parms$fb_oc, "\n",
               parms$fb_fpoc, "\n",
               parms$fb_ow, "\n",
               parms$duration, "\n",
               parms$injmass, "\n",
               parms$assfb, "\n",
               parms$casename, "\n") 
  write(str, "config.txt")
}

run.model = function(parms){
  if(!dir.exists(parms$casename)){
    dir.create(parms$casename)
  }
  
  system("fernSpike.exe")
}

read.data = function(parms){
  ifname = file.path(parms$casename, "exogenic.txt")
  d = read.table(ifname, header = TRUE, fill = TRUE)
  d = head(d, nrow(d) - 3)
  for(i in 1:3){
    d[, i] = as.numeric(d[, i])
  }
  return(d)
}

sim = function(parms){
  write.config(parms)
  run.model(parms)
  d = read.data(parms)
  return(d)
}

plot.case = function(case, base, casename){
  fname = paste0(casename, "_d13C.png")
  png(file.path("fernSpike", "figs", fname), width = 6, height = 8, units = "in", res = 600)
  layout(matrix(1:3, nrow = 3))
  par(mar = c(5, 5, 1, 5))
  plot(case$Time - 20000, case$d13CAtm, type = "l", lwd = 3, axes = FALSE, 
       xlab = "", ylab = "")
  lines(base$Time - 20000, base$d13CAtm)
  axis(2)
  mtext(expression(delta^{13}*"C"[atm]), 2, 3)
  plot(case$Time - 20000, case$d13CSurf, type = "l", lwd = 3, axes = FALSE, 
       xlab = "", ylab = "")
  lines(base$Time - 20000, base$d13CSurf)
  axis(4)
  mtext(expression(delta^{13}*"C"[surf]), 4, 3)
  plot(case$Time - 20000, case$d13CDeep, type = "l", lwd = 3, axes = FALSE, 
       xlab = "", ylab = "")
  lines(base$Time - 20000, base$d13CDeep)
  axis(2)
  axis(1)
  mtext("Time (kyr)", 1, 3)
  mtext(expression(delta^{13}*"C"[deep]), 2, 3)
  dev.off()
  
  fname = paste0(casename, "_CC.png")
  png(file.path("fernSpike", "figs", fname), width = 6, height = 8, units = "in", res = 600)
  layout(matrix(1:3, nrow = 3))
  par(mar = c(5, 5, 1, 5))
  plot(case$Time - 20000, case$Inject * 12e3, type = "l", lwd = 3, axes = FALSE, 
       xlab = "", ylab = "")
  lines(base$Time - 20000, base$Inject * 12e3)
  axis(2)
  mtext("C injection (Pg/yr)", 2, 3)
  plot(case$Time - 20000, case$PCO2, type = "l", lwd = 3, axes = FALSE, 
       xlab = "", ylab = "")
  lines(base$Time - 20000, base$PCO2)
  axis(4)
  mtext(expression("pCO"[2]*" (x preindustrial)"), 4, 3)
  plot(case$Time - 20000, case$BioC * 12e3, type = "l", lwd = 3, axes = FALSE, 
       xlab = "", ylab = "")
  lines(base$Time - 20000, base$BioC * 12e3)
  axis(2)
  axis(1)
  mtext("Time (kyr)", 1, 3)
  mtext("Terrestrial biosphere (Pg)", 2, 3)
  dev.off()
}


