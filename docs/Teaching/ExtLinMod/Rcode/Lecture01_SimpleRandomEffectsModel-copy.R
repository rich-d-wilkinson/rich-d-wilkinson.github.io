## ---- echo=FALSE, warning=FALSE, message=FALSE---------------------------
library(lme4)
Cropdata <- Dyestuff
colnames(Cropdata) <- c('Farm', 'Yield')

## ------------------------------------------------------------------------
str(Cropdata)
head(Cropdata)

## ---- fig.width=5, fig.height=3------------------------------------------
library(ggplot2)
qplot(Farm, Yield, geom='boxplot', data = Cropdata)

## ---- fig.width=5, fig.height=3------------------------------------------
qplot(reorder(Farm, Yield), Yield, geom='boxplot', data = Cropdata, xlab = 'Farm')

## ------------------------------------------------------------------------
(fit1 <- lm(Yield~Farm-1, data=Cropdata))
#model.matrix(fit1) # useful for checking what model we have fit.

## ---- eval=FALSE---------------------------------------------------------
## lm(Yield~Farm, data=Cropdata)

## ------------------------------------------------------------------------
(fit2 <- lm(Yield~Farm, data=Cropdata, contrasts=list(Farm=contr.sum)))

## ------------------------------------------------------------------------
library(lme4)
(fm02 <- lmer(Yield ~ 1+(1|Farm), Cropdata))

## ------------------------------------------------------------------------
summary(fm02)
summary(fit2)

## ------------------------------------------------------------------------
ranef(fm02)

## ------------------------------------------------------------------------
library(lme4)
str(sleepstudy)
head(sleepstudy)

## ---- fig.width=7, fig.height=4------------------------------------------
library(ggplot2)
qplot(Days, Reaction, facets=~Subject, data = sleepstudy)

