###################################################################
### jsp data
###################################################################

library(faraway)
library(lme4)
library(ggplot2)
data(jsp)
head(jsp)


####
# select final year only
jspr <- jsp[jsp$year==2,]

qplot(social, math, data = jspr, geom='boxplot')

qplot(school, math, data = jspr, geom='boxplot')


## Lets account for school

lmer(math~ social + (1|school), data=jspr, contrasts=list(social=contr.sum))


# To account for class we need to nest effects

lmer(math~ social + (1|school)+(1| school:class), data=jspr, contrasts=list(social=contr.sum))

# or equivalently
lmer(math~ social + (1| school/class), data=jspr, contrasts=list(social=contr.sum))

####################################################################



load('MAS473.RData')
#########################################
####### Chapter 3 Section 3.1.1 #########

# Inspect the data
Oxide

# Note that Wafer isn't a 3-level factor; Wafer is nested within Lot
# eg Y_{121} and Y_{221} do not refer to the same Wafer.

library(ggplot2)
qplot(data=Oxide, Wafer, Thickness, facets=~Lot)
qplot(data=Oxide, Wafer, Thickness, facets=~Lot, geom='boxplot')
# which do you prefer?

qplot(data=Oxide,Lot,Thickness,geom='boxplot')

# We can see variation between Wafers in the same Lot
# and variation between Lots


fm1<-lmer(Thickness~1+(1|Lot/Wafer),data=Oxide)
summary(fm1)

# equivalent command
fm2 <- lmer(Thickness~1 +(1|Lot)+(1|Lot:Wafer), data =Oxide)
summary(fm2)
# Estimated random effects
ranef(fm1)



###############################
####### Section 3.1.2 #########


Oats
attach(Oats)
qplot(Block, yield, geom='boxplot')
qplot(nitro, yield)
qplot(Variety,yield, geom='boxplot')

qplot(Variety:Block, yield, geom='boxplot', col=Variety)

# Can get the plot random effect via Block/Variety
# Variety included as a fixed effect too
(fm1<-lmer(yield~nitro+Variety+(1|Block/Variety),Oats))
# doesn't mean that variety is a random effect, just that it indicates the different plots within the block
# note that we don't need a random effect for the subplot - only 1 observation in each so represented by epsilon_{ijk}

# If this confuses you, can specify a plot factor:
Plt<-gl(18,4)
Plt
Oats<-data.frame(Oats,Plt)

# Same model as fm1
(fm1b<-lmer(yield~nitro+Variety+(1|Block)+(1|Plt),Oats))

# Model is
# Y_{ijk} = mu + tau_{v(i,j)} + beta x_{ijk} + b_i + b_{ij} + epsilon_{ijk}
# i=1,...,6: block
# j=1,2,3: plot, 
# k=1,2,3,4: subplot
# v(i,j) = 1,2,3 corresponding to variety (eg v(1,2)=2 for "Golden Rain") 
# x_{ijk}: nitro level


summary(fm1b)
ranef(fm1b)

# Note: can't have an ordinary linear model with plot effect, as confounded with Variety
(lm1<-lm(yield~nitro+Variety+Plt ,Oats, contrasts=list(Plt=contr.sum)))

# The 'equivalent' fixed effects model can be achieved by having a Variety:Block interaction
# Gives same estimates for Variety effects, but with smaller standard errors.
lm1<-lm(yield~nitro+Variety*Block ,Oats, contrasts=list(Block=contr.sum))
summary(lm1)

