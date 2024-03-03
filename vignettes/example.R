# This is our external R script called example.R
# We're adding two chunks variablesXY and plotXY

## @knitr variablesXY

x<-1:100
y<-x+rnorm(100)
head(data.frame(x,y))

## @knitr plotXY

plot(x,y)
