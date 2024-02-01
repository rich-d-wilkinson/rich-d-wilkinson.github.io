
################################################
# Q1.
###################################################################################################

plot(iris$Sepal.Length, iris$Sepal.Width)
library(ggplot2)
qplot(Sepal.Length, Sepal.Width, data=iris)


library(dplyr)

# Note there is already a `filter` command. Loading dplyr masks the other one. If in doubt, use
dplyr::filter

iris[,1]+iris[,2]

iris |> select(Sepal.Length, Sepal.Width)





iris |> mutate(Sepal.sum=Sepal.Length+Sepal.Width)
# ---

mean(iris[iris$Species=='setosa', 'Petal.Length'])


mean(iris[iris$Species=='setosa', 'Petal.Length'])


iris |> filter(Species=='setosa') |>
  select(Petal.Length) |> 
  summarise(mean(Petal.Length))

iris |>
  group_by(Species) |>
  summarise( mean(Petal.Length))

# ----------
iris[iris$Petal.Length>5,]
summary(iris[iris$Petal.Length>5,'Species'])/dim(iris[iris$Petal.Length>5,])[1]

iris |>
  filter(Petal.Length > 5)

# Lots of ways to do this.
iris |>
  filter(Petal.Length > 5)|>
  group_by(Species) |>
  summarise(n = n()) |>
  mutate(freq = n / sum(n))

iris |>
  filter(Petal.Length > 5) |>
  count(Species) |>
  mutate(freq = n / sum(n))

# -----

iris[iris$Sepal.Length>6& iris$Petal.Length<5,]
summary(iris[iris$Sepal.Length>6& iris$Petal.Length<5,'Species'])/dim(iris[iris$Sepal.Length>6& iris$Petal.Length<5,])[1]


iris |> filter(Sepal.Length>6, Petal.Length<5) |>
  count(Species) |>
  mutate(freq =n/sum(n))


# ---------------------
# Data frame vs matrix is an annoyance we have to deal with.

as.matrix(iris[,1:4]) %*% diag(1:4)

iris |> select_if(is.numeric) |> as.matrix %*% diag(1:4)

###################################################################################################
# 2.
###################################################################################################
library(dplyr)
Ex1 <- data.frame(
  Student=LETTERS[1:5],
  P = c(41,72,46,77,59),
  S = c(63,82,38,57,85)
)

Ex1 |> select_if(is.numeric) |> colMeans()

(mu = colMeans(Ex1[,2:3]))

Ex1 |> select_if(is.numeric) |> cov()
cov(Ex1[,2:3])

# Covariance by hand
X<- as.matrix(Ex1[,2:3])
Xcent<- sweep(X, 2, mu) #remove column means



t(Xcent)%*%Xcent/4


###########################################################
#3.
###################################################################################################

# https://rpubs.com/chrisbrunsdon/spacetime_1a
mtcars2 <- within(mtcars, {
  vs <- factor(vs, labels = c("V", "S"))
  am <- factor(am, labels = c("automatic", "manual"))
  cyl  <- ordered(cyl)
  gear <- ordered(gear)
  carb <- ordered(carb)
})

library(ggplot2)
library(GGally)

pairs(mtcars2)
ggpairs(mtcars2, 
        columns=c('mpg', 'cyl', 'disp', 'hp', 'drat', 'wt', 'qsec'), 
        mapping=ggplot2::aes(colour = am))


ggpairs(mtcars2, columns=c('mpg', 'cyl', 'disp', 'hp', 'drat', 'wt', 'qsec'))

ggpairs(mtcars2, columns=c('mpg', 'cyl', 'disp', 'hp', 'drat', 'wt', 'qsec'), 
        mapping=ggplot2::aes(colour = gear))



ggplot(mtcars2,aes(x=cyl,y=mpg)) + geom_boxplot() + 
  xlab('Cylinders') + ylab('Miles per Gallon')

ggplot(data=mtcars2, aes(x=gear, y=mpg))+geom_boxplot()

#qplot(disp, mpg,data=mtcars2,facets=cyl)

ggplot(mtcars2,aes(x=wt,y=mpg)) + geom_point()+ 
  xlab('Weight (x 1000lbs)') + ylab('Miles per Gallon') + geom_smooth()

ggplot(mtcars2,aes(x=wt,y=mpg,col=cyl)) + geom_point() + 
  labs(x='Weight (x1000lbs)',y='Miles per Gallon',colour='Number of\n Cylinders')

ggplot(mtcars2,aes(x=wt,y=mpg,col=cyl)) + geom_point() + 
  facet_grid(~cyl)

ggplot(mtcars2,aes(x=wt,y=mpg,col=cyl)) + geom_point() + facet_grid(am~cyl)

###########################################################
# 4.
###################################################################################################

library(mvtnorm)
mu = c(1,0)
Sigma=matrix(c(2,1,1,2), nr=2)

X <- rmvnorm(n=100, mean=mu, sigma=Sigma)
colMeans(X)
cov(X)

set.seed(1)
X <- rmvnorm(n=100, mean=mu, sigma=Sigma)
colMeans(X)
cov(X)



X <- rmvnorm(n=10^6, mean=mu, sigma=Sigma)
colMeans(X)
cov(X)


########################################################################
# 5.
###################################################################################################


load('mnist.rda')
mnist$train ## a training set of 60000 images
mnist$test ## a test set of 10000 images
mnist$train$x # image intensities
mnist$train$y # image labels

mnist$train$x[1,]


library(reshape2)
library(ggplot2)


plot.mnist <- function(im){
  #im[im<0]<-0 # set any negative intensities to zero
  #im[im>1]<-1 # set an intensities bigger than 1 to 1.


  if(is.vector(im)){ # a single image

    A<-matrix(im, nr=28, byrow=F)
    C<- melt(A, varnames = c("x", "y"), value.name = "intensity")
    p<-ggplot(C, aes(x = x, y = y, fill = intensity))+
      geom_tile(aes(fill=intensity))+
      scale_fill_gradient(low='white', high='black')+
      scale_y_reverse()+theme(
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )
  }
  else{
    if (dim(im)[2]!=784){
      im = t(im)
    }
    n <- dim(im)[1]
    As <- array(im, dim = c(n, 28, 28))

    Cs<- melt(As, varnames = c("image","x", "y"), value.name = "intensity")
    p<-ggplot(Cs, aes(x = x, y = y, fill = intensity))+
      geom_tile(aes(fill=intensity))+
      scale_fill_gradient(low='white', high='black')+
      facet_wrap(~ image, nrow = floor(sqrt(n))+1, ncol = floor(sqrt(n))+1)+
      scale_y_reverse()+theme(
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )

  }
  return(p)
}

plot.mnist(mnist$train$x[77,])

plot.mnist(mnist$train$x[1:10,])

mnist5 <- mnist$train$x[mnist$train$y==5,]

plot.mnist(mnist5[20:30,])
