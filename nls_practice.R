#using nls() instead of nls.multstart()

x <- c(32,64,96,118,126,144,152.5,158)
y <- c(99.5,104.8,108.5,100,86,64,35.3,15)
#we will make y the response variable and x the predictor
#the response variable is usually on the y-axis
plot(x,y,pch=19)

#fit first degree polynomial equation:
fit  <- lm(y~x)
#second degree
fit2 <- lm(y~poly(x,2,raw=TRUE))
#third degree
fit3 <- lm(y~poly(x,3,raw=TRUE))
#fourth degree
fit4 <- lm(y~poly(x,4,raw=TRUE))
#generate range of 50 numbers starting from 30 and ending at 160
xx <- seq(30,160, length=50)
plot(x,y,pch=19,ylim=c(0,150))
lines(xx, predict(fit, data.frame(x=xx)), col="red")
lines(xx, predict(fit2, data.frame(x=xx)), col="green")
lines(xx, predict(fit3, data.frame(x=xx)), col="blue")
lines(xx, predict(fit4, data.frame(x=xx)), col="purple")

summary(fit)

anova(fit,fit2)

#y=ax + b
coef(fit)
#y=ax^2 + bx + c
coef(fit2)
#y=ax^3 + bx^2 + cx + d
coef(fit3)
#y=ax^4 + bx^3 + cx^2 + dx + e
coef(fit4)

#create function for the third order polynomial
third_order <- function(newdist, model) {
  coefs <- coef(model)
  #y = d + cx + bx^2 + ax^3
  res <- coefs[1] + (coefs[2] * newdist) + (coefs[3] * newdist^2) + (coefs[4] * newdist^3)
  return(res)
}
yy <- third_order(xx,fit3)
plot(xx,yy,type="l")
lines(x,y,col="red")

xx<-seq(-1000,1000,length=500)
yy<-third_order(xx,fit3)
plot(xx,yy,type="l")

