# Supplementary R code for lecture 4, STA303 summer 2018
# Alex Stringer
# 

library(tidyverse)
library(faraway)


# Load and plot the orings data

data(orings)
orings_tbl <- as_data_frame(orings) %>%
  mutate(prop_damaged = damage / 6) # Maximum likelihood estimate

glimpse(orings_tbl)

orings_tbl %>%
  ggplot(aes(x = temp,y = prop_damaged)) +
  theme_classic() +
  geom_point(pch=21) +
  geom_smooth(method = "lm",se = FALSE,colour="blue") +
  labs(title = "Orings Data",
       subtitle = "Probability of o-ring failure for 23 space shuttle missions as a function of temperature on launch day",
       x = "Temperature",
       y = "Probability of Damage") +
  scale_y_continuous(labels = scales::percent_format())


# Fit a binomial GLM
# family = binomial ==> binomial regression
# family = poisson ==> poisson regression
# family = gaussian ==> same as lm()
# 
# Standard errors got from likelihood theory
# NOT exact
# They are off in same cases as deviance is off
# ...including binomial regression when m_i are small.
# So, here.
# 
# Residual deviance: deviance as defined in lecture.
# -2log(LRT) between full and SATURATED models
# 
# Null deviance: deviance of the null model,
# So -2log(LRT) between null and residual
#
# Goal was: find a model that fits "as well as"
# saturated model, or at least better than null.
# So, look at (null deviance) - (residual deviance)-
# "amount of deviance explained" by the fitted model
#
# How big is big enough?
#
# (null) - (residual) approx chi squared, df = dim(null) - dim(fitted model)
# Here, this equals 22 - 21 = 1
# This also corresponds to the fact that we fit a model with one parameter
# more than the intercept
#
# Mean of chisquare(n) = n
# Variance = 2n
# So quick check: if model fits as well as saturated,
# Residual deviance won't be more than 2sqrt(2df) away from df
# and if fits better than null, then (null - residual) will be
# more than 2sqrt(2df) from its df

glm1 <- glm(cbind(damage,6-damage) ~ temp,data=orings_tbl,family=binomial)
summary(glm1)

# Don't believe me? (Actually it's because I'm not sure...)
glm1null <- glm(cbind(damage,6-damage) ~ 1,data=orings_tbl,family=binomial)
summary(glm1null)

# Plot the fitted curve
orings_tbl %>%
  mutate(predicted_prob = ilogit(coef(glm1)[1] + coef(glm1)[2]*temp)) %>%
  ggplot(aes(x = temp,y = prop_damaged)) +
  theme_classic() +
  geom_point(pch=21) +
  geom_line(aes(y = predicted_prob),colour="blue",size = 1) +
  labs(title = "Orings Data - Fitted Binomial Regression Model",
       subtitle = "Probability of o-ring failure for 23 space shuttle missions as a function of temperature on launch day",
       x = "Temperature",
       y = "Probability of Damage") +
  scale_y_continuous(labels = scales::percent_format())

# Looks better. Can we extrapolate below the range of observed data?

orings_tbl %>%
  ggplot(aes(x = temp,y = prop_damaged)) +
  theme_classic() +
  geom_point(pch=21) +
  geom_line(data = data_frame(
    temp = seq(25,90,by=0.1),
    prop_damaged = ilogit(coef(glm1)[1] + coef(glm1)[2]*temp)
    ),colour="blue",size = 1) +
  labs(title = "Orings Data - Fitted Binomial Regression Model",
       subtitle = "Probability of o-ring failure for 23 space shuttle missions as a function of temperature on launch day",
       x = "Temperature",
       y = "Probability of Damage") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = seq(30,80,by=10))


# Test the hypothesis that temperature as no effect on probability of failure
# Compare the model deviance with the null deviance
# 
# In R, type str(glm1) to see all the stuff available
# in the object glm1

D_temp <- glm1$null.deviance - deviance(glm1)
D_temp
# Compare to a chisq(2 - 1) distribution
1 - pchisq(D_temp,1)
# full better

# Confidence interval for beta(temp)
betatemp_se <- sqrt(diag(vcov(glm1)))[2]
cint_temp <- c(coef(glm1)[2] - 1.96 * betatemp_se,coef(glm1)[2] + 1.96 * betatemp_se)
cint_temp
# What about a confidence interval for the factor by which odds increases for a unit increase in temp?
# That's exp(beta), by the way
exp(cint_temp)


# Predict at a new temp
# E.g. 31F, the launch temperature on the day Challenger exploded

xnew <- cbind(c(1,31)) # Intercept
etapred <- t(xnew) %*% cbind(coef(glm1))
# What is that? That's eta!
etapred
# Predicted probability: calculate pi
ilogit(etapred)
# About 99%- yikes
# Confidence interval for the prediction?
etapred_sd <- sqrt(t(xnew) %*% vcov(glm1) %*% xnew)
etapred_sd
eta_confint <- c(etapred - 1.96 * etapred_sd,etapred + 1.96 * etapred_sd)
eta_confint
ilogit(eta_confint)

# Compare this with the output of the predict() function
predict(glm1,newdata = data_frame(temp = 31),se.fit = TRUE,type = "link")
predict(glm1,newdata = data_frame(temp = 31),se.fit = TRUE,type = "response")


### Binary logistic regression: complete example ###

# We'll analyze the Wisconsin Breast Cancer Data, exercise 2 chapter 2 page 58 from ELMR
# Note the data is misspelled in the Faraway package
data(wbca) # Called wbcd in the book
wbca_tbl <- as_data_frame(wbca)

wbca_tbl

# First step: take a look at the response, Class, which is an indicator of malignancy of
# tumours (1 == benign, 0 == malignant)
with(wbca_tbl,table(Class))
443 / (443 + 238)
238 / (443 + 238)


# We have 9 predictors, so it's feasible to plot the data separately for each predictor
# We'll do this in a clever, "tidy" way

wbca_tbl %>%
  gather(variable,value,Adhes:USize) %>%
  ggplot(aes(x = value,y = Class)) +
  theme_classic() +
  facet_wrap(~variable) +
  geom_jitter(width=0.5,height=0.1,alpha = 0.3) +
  scale_y_continuous(breaks = c(0,1),labels = c("Malignant","Benign")) +
  labs(title = "Malignancy of Tumours, by Predictor",
       subtitle = "WBCA Data",
       x = "Predictor Value",
       y = "Malignant or Benign?")

# We have several predictors, so it is a good idea to look at their correlation matrix
round(cor(wbca_tbl %>% dplyr::select(-Class)),2)
corrplot::corrplot(cor(wbca_tbl %>% dplyr::select(-Class)),order="AOE")
# There are some highly correlated variables, especially UShap and USize

# Fit an initial binary logistic regression model
# with all x
wbca_glm1 <- glm(Class ~ .,data = wbca_tbl,family = binomial)
summary(wbca_glm1)

# The coefficients all look reasonable
# What about the residual deviance?
# ALL THESE QUANTITIES ARE USELESS FOR BINARY REGRESSION!
# But you can still compute them... but you shouldn't
Dstat <- wbca_glm1$null.deviance - deviance(wbca_glm1)
Dstat
# Compare to a chisq(680 - 671) distribution
1 - pchisq(Dstat,9)

# So the model fits better than the null model
# NO! CAN'T SAY THAT
# Compare to the saturated model
# But wait... what is the deviance of the saturated model for BINARY logistic regression? 
# 0 since likelihood =1
# deviance = -2log(L)
y <- wbca_tbl %>% pull(Class)
dev_sat <- sum( y*log(y) + (1 - y)*log(1 - y))
dev_sat
# What's the problem?
# 
# Take dev_sat = 0 (again, why?)
# Then a test of whether the model fits as well as the saturated model is simply a test of comparing
# the residual deviance to its degrees of freedom
deviance(wbca_glm1)
1 - pchisq(deviance(wbca_glm1),671)
# Small residual deviance ==> model fits the data "as well" as the saturated model... usually
# But here, we cannot use residual deviance, because it doesn't depend on the data. Only on p-hat
# Note though that the chisquare approximation is terrible for binary data, and there are other problems in using deviance with binary data
# Let's make sure the null model doesn't also fit the data as well as the saturated model
1 - pchisq(wbca_glm1$null.deviance,680)
# But don't trust these comparisons, for BINARY data
# The above would be a good thing to do for binomial data with bigger n per observation
# Or for count data
# But NOT for BINARY data
# 
# Now. Can we get a simpler model that fits just as well?

wbca_glm2 <- step(wbca_glm1,
                  lower = glm(Class~1,data=wbca_tbl,family=binomial))
summary(wbca_glm2)

# We removed Epith and USize. Remember that USize was highly correlated with UShap.
# It's likely that either of these variables would have been fine; e.g. if there are business/scientific
# reasons for wanting to include one over the other, you probably could choose
# 
# What next? Let's look at the actual predicted probabilities. The variables in the model were all pretty skewed, so what do
# we expect the distribution of predicted probabilities to look like?

# FYI, you should be able to compute the below from scratch, i.e. not using predict()
wbca_predicted_probs <- predict(wbca_glm2,type="response")

data_frame(x = wbca_predicted_probs) %>%
  ggplot(aes(x = x)) +
  theme_classic() +
  geom_histogram(bins = 100,colour = "black",fill = "orange") +
  labs(title = "Histogram of Predicted Probabilities",
       subtitle = "Predicted probability of tumour being benign, logistic regression on WBCA data",
       x = "Predicted Probability",
       y = "# of tumours") +
  scale_x_continuous(labels = scales::percent_format())

# Most of the predicted probabilities are near 0, or near 1
# How to make a hard 0/1 prediction for each case? The textbook suggests using a cutoff of 0.5.
# Let's look at the classification error if we do that

wbca_tbl %>%
  dplyr::select(Class) %>%
  bind_cols(data_frame(predprob = wbca_predicted_probs)) %>%
  mutate(ypred = as.numeric(predprob > .5)) %>%
  group_by(Class,ypred) %>%
  summarize(cnt = n(),
            pct = scales::percent(cnt / nrow(wbca_tbl)))

# What about using .9 as a cutoff?

wbca_tbl %>%
  dplyr::select(Class) %>%
  bind_cols(data_frame(predprob = wbca_predicted_probs)) %>%
  mutate(ypred = as.numeric(predprob > .9)) %>%
  group_by(Class,ypred) %>%
  summarize(cnt = n(),
            pct = scales::percent(cnt / nrow(wbca_tbl)))
# FP, FN bu neng tong shi bian xiao
# False Positive: ACTUAL negative, PREDICT positive
# False Negative: ACTUAL positive, PREDICT negative
# True Positive: ACTUAL positive, PREDICT positive
# True Negative: ACTUAL negative, PREDICT negative
# 
# False positive rate: FP / (All actual negatives)
# = FP / (FP + TN)
# True positive rate: TP / (All actual positives)
# = TP / (TP + FN)
# 
# Strongly recommend you read the wikipedia article too.
# And then copy this out 100 times.

# We see there is a tradeoff- When we increase the cutoff, we get less false positives, and more false negatives
# Since a positive here means benign, a false positive is classifying a malignant tumour as being benign, which seems worse
# We can plot a parametric curve of the false positives vs true positives for cutoffs ranging from 0 to 1- the ROC curve
# A wide ROC curve means there exists favourable cutoffs, i.e. the model is good
# 
# I won't require you to be able to create the below from scratch like I am here
# But you do need to know how to interpret the below graph.

fpr <- function(y,ypred) sum(ypred==1 & y==0) / sum(y==0)
tpr <- function(y,ypred) sum(ypred==1 & y==1) / sum(y==1)

wbca_with_predprob <- wbca_tbl %>%
  dplyr::select(Class) %>%
  bind_cols(data_frame(predprob = wbca_predicted_probs))

roc_data <- wbca_with_predprob %>%
  arrange(predprob) %>%
  pull(predprob) %>%
  round(3) %>%
  unique() %>%
  map(~c(.x,tpr(y = wbca_with_predprob$Class,ypred = as.numeric(wbca_with_predprob$predprob >= .x)),fpr(y = wbca_with_predprob$Class,ypred = as.numeric(wbca_with_predprob$predprob >= .x)))) %>%
  map(~.x %>% as_data_frame() %>% t() %>% as_data_frame()) %>%
  reduce(bind_rows) %>%
  rename(h=V1,tpr=V2,fpr=V3)

roc_data %>%
  ggplot(aes(x=fpr,y=tpr)) +
  theme_light() +
  geom_path(colour="#195FFF",size=1) +
  scale_x_continuous(labels=scales::percent_format()) +
  scale_y_continuous(labels=scales::percent_format()) +
  labs(title="ROC Curve",
       subtitle = "Wisconsin Breast Cancer Data - Logistic Regression Model",
       x="FPR",
       y="TPR") +
  geom_line(data = data_frame(x = c(0,1),y = c(0,1)),aes(x=x,y=y),linetype="dashed",colour="red")

# Every point on this curve is a (FPR,TPR)
# combination for a given cutoff
# Answers: can we pick a cutoff that separates the 0s and 1s?
# And, what FPR/TPR correspond to a given cutoff?
