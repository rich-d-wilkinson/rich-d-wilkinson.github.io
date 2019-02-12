## ------------------------------------------------------------------------
library(lme4)
str(sleepstudy)
head(sleepstudy)

## ------------------------------------------------------------------------
library(ggplot2)
qplot(Days, Reaction, facets=~Subject, data = sleepstudy)

## ------------------------------------------------------------------------
c <- ggplot(sleepstudy, aes(Days, Reaction))+geom_point() + facet_wrap(~ Subject, nc=5)
c+stat_smooth(method='lm')

## ------------------------------------------------------------------------
(fm06 <- lmer(Reaction ~ 1 + Days + (1 +Days | Subject), sleepstudy))

## ------------------------------------------------------------------------
(fm07 <- lmer(Reaction ~ 1 + (1 | Subject)+ Days + (Days-1|Subject), sleepstudy))

