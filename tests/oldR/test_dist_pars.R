library(ggplot2)
xfit<-seq(0,10,length=40)
df = data.frame(x = seq(-10,10,length=40),
                y = log(dnorm(xfit,mean=0,sd=1)),
                z = lognorm_dens(xfit, 0, 1))
ggplot(df) + aes(x, y) + geom_line(color = "blue") + 
  aes(x, z) + geom_line(color = "red")
