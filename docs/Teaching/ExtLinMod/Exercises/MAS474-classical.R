# Section 3.3. Example to show that standard fixed effects model
# will not handle variation in data properly, if one factor is 'random'

# Choose number of patients
n.patients<-2

# Choose number of replicate observations per patient/drug combination. Try something large.
n.reps<-100

# Create a vector of factor levels for patient
patient<-gl(n.patients,2*n.reps)

# Create a vector of factor levels for drug
drug<-gl(2,n.reps,2*n.reps*n.patients,labels=c("A","B"))

# Choose the true values of the mean response on each drug
drug.A<-5
drug.B<-0

# Generate some random patient effects, then arrange into a suitable vector for
# simulating the data. 
patient.effect<-rnorm(n.patients,0,2)
patient.effect.vector<-rep(patient.effect,each=2*n.reps)

# Generate some random patient-drug interaction effects, then arrange into a suitable vector for
# simulating the data
interaction.effect<-rnorm(n.patients*2,0,1)
interaction.effect.vector<-rep(interaction.effect,each=n.reps)

# Generate the random errors
errors<-rnorm(n.patients*n.reps*2,0,0.5)

# Calculate the vector of observations
response<-drug.A*(drug=="A")+drug.B*(drug=="B")+patient.effect.vector+interaction.effect.vector+errors

# Fit a standard fixed effects two-way anova model
lm1<-lm(response~patient*drug,contrasts=list(drug=contr.sum,patient=contr.sum))

# Calculate confidence intervals for the model parameters
# Note that the"drug1" term corresponds to (drug.A-drug.B)/2, so the the true value is 2.5. 
confint(lm1)
# You should observe that even with only 2 patients, if n.reps is large
# you'll get a narrow confidence interval around the 'wrong' value.
# This is because the drug1 parameter is the average for the 2 patients, not the population. 


#### Calculating a confidence interval for the population drug effect

# Equation 2.16, (Same as the 'drug1' term in the fixed effects model)
effect<-0.5*(mean(response[drug=="A"])-mean(response[drug=="B"]))

# We can get the variance estimates from the anova table
anova.table<-anova(lm1)

(sigmasq<-anova.table[4,3]) # equation (2.13)
(sigmasq.patient.drug<-(anova.table[3,3]-sigmasq)/n.reps) # equation (2.15)

# Calculate variance (then take square root) of estimator given in (2.17)
N<-n.patients*n.reps*2
(effect.se<-(sigmasq/N+ 0.5*sigmasq.patient.drug/n.patients )^0.5)

# The 95% confidence interval - much wider if only 2 patients
c(effect-2*effect.se, effect+2*effect.se)

##### One you've understood how to fit mixed effects models in R, you can try
##### the equivalent analysis with a mixed effects model

library(lme4)

(fm1<-lmer(response~drug +(1|patient/drug),contrasts=list(drug=contr.sum)))


# Should get the same std error for the estimate of the drug1 term, and
# the same estimates of the variance components