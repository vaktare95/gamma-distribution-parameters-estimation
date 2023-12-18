plot_theo_sim = function(par, theo, sim, title, xlab, ylab) {
  plot(par, theo, type="l", main=title, col="red", xlab=xlab, ylab=ylab)
  points(par, sim, col="blue")
  legend("topright", inset=c(0,-0.2), xpd=TRUE, legend=c("theo","sim"), col=c("red", "blue"), lty=c(1,NA), pch=c(NA,'O'))
}