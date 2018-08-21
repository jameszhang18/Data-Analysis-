rm(list = ls())
library(tidyverse)
library(ggplot2)



#######################QUESTION 1 ###############
data = read.table(file="https://www.stat.berkeley.edu/~statlabs/data/gauge.data", header=T)
data_t <- as_data_frame(data)
glimpse(data_t)

# Check the observation counts in each level of the factors
data_t %>%
  group_by(gain,density) %>%
  summarize(cnt = n())


## linear model
ggplot(data_t, aes(x=gain, y=density)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE) +
  theme_classic() +
  labs(title = "Linear Regression Model",
       x = "Gain",
       y = "Density"
  )


fit1<-lm(density~gain,data=data_t)
summary(fit1)
anova(fit1)
library(ggfortify)
#diagnostic plots
autoplot(fit1)

#boxplot
ggplot(data_t, aes(x = gain,y = density,group=1)) +
  theme_classic() +
  geom_boxplot() +
  labs(title = "Boxplot of density by gain",
       x = "gain",
       y = "density"
  )


#boxplot for residuals
boxplot(residuals(fit1), main="boxplot for residual")


#####TODO###########
library(MASS)
boxcox(fit1, family="yjPower", plotit = TRUE)
density_boxcox <- boxcox(fit1,lambda = seq(-2,2,by=0.01))
density_boxcox$x[which.max(density_boxcox$y)]

density_boxcox <- boxcox(density ~ gain,data=data_t)
density_boxcox$x[which(density_boxcox$y == max(density_boxcox$y))]


data_tbl_transform <- data_t %>%
  mutate(log_density = log(density))
##################

# try to plot the data and transfer gain to log(gain), then fit it
fit2 <-lm(density~I(log(gain)),data=data_t)
ggplot(data_t, aes(x=I(log(gain)), y=density)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm,   # Add linear regression line
              se=FALSE) +
  theme_classic() +
  labs(title = "Log Regression Model",
       x = "Gain",
       y = "Density"
  )

anova(fit2)
summary(fit2)
#boxplot for residuals
boxplot(residuals(fit2), main="boxplot for residual")
# standardized residuals vs. fitted values 
autoplot(fit2)

#
#boxplot
ggplot(data_t, aes(x = I(log(gain)),y = density,group=1)) +
  theme_classic() +
  geom_boxplot() +
  labs(title = "Boxplot of density by gain",
       x = "gain",
       y = "density"
  )
##
##

## Regression for polynomial regression degree 2
fit3 <-lm(density~poly(gain,2),data=data_t)
anova(fit3)
summary(fit3)
autoplot(fit3)




#plot(gain,density,main="regression for polynomial regression degree 2")
pts=seq(min(gain),max(gain),len=100)
val=predict(fit3,data.frame(gain=pts))
lines(pts,val,col="red",lwd=3)
#boxplot for residuals
boxplot(residuals(fit3), main="boxplot for residual")
# standardized residuals vs. fitted values 
plot(fit3, which=3)




##################QUESTION 2######################
crabs<-readr::read_table("https://www.stat.berkeley.edu/users/statlabs/data/crabpop.data",col_types = "dc")
crabs_t<-as_data_frame(crabs)
glimpse(crabs_t)

# Check the observation counts in each level of the factors
crabs_t %>%
  group_by(shell) %>%
  summarize(cnt = n())

#boxplot
ggplot(crabs_t, aes(x = shell,y = size)) +
  theme_classic() +
  geom_boxplot() +
  labs(title = "Boxplot of size by shell",
       x = "shell",
       y = "size"
  )
## Have different mean



# Grand mean and standard deviation
summarize(crabs_t,grand_mean = mean(size),grand_sd = sd(size))


# Group mean
group_by(crabs_t, shell) %>%
  summarise(
    count = n(),
    median = median(size, na.rm = TRUE),
    mean = mean(size, na.rm = TRUE),
    sd = sd(size, na.rm = TRUE)
  )

#anova table
crabs_aov <-aov(size ~ shell, data=crabs_t)
summary(crabs_aov)

# As the p-value is less than the significance level 0.05,
# we can conclude that there are significant differences between the groups highlighted with "*" in the model summary.


# Normal QQ
fit4 <-lm(size ~ shell,data=crabs_t)
e <- fit4$residuals
crabs_t %>%
  mutate_at("e",funs( (. - mean(.)) / sd(.))) %>%
  arrange(e) %>%
  mutate(q = qnorm(1:n() / (n() + 1))) %>%
  ggplot(aes(x = q,y = e)) +
  theme_classic() +
  geom_point() +
  geom_abline(slope = 1,intercept = 0) +
  labs(title = "Normal QQ-plot",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")

## not normal -> two-sample t-test

# 1. Homogeneity of variances
plot(res.aov, 1)
## The ANOVA test assumes that, the data are normally distributed and the variance across groups are homogeneous.

# two-sample t-test
# Assumption 1: Are the two samples independents?
# Assumtion 2: Are the data from each of the 2 groups follow a normal distribution?
# Shapiro-Wilk normality test for shell=0
with(crabs_t, shapiro.test(size[shell == 0]))# p = 0.1
# Shapiro-Wilk normality test for shell=1
with(crabs_t, shapiro.test(size[shell == 1])) # p = 0.6
# Assumption 3. Do the two populations have the same variances?
var.test(size ~ shell, data = crabs_t)

t.test(size ~ shell, data = crabs_t, var.equal = TRUE)

t.test(size ~ shell, data = crabs_t,
       var.equal = TRUE, alternative = "greater")


