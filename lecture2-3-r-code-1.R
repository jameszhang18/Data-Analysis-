# Supplementary R code for lectures 2 and 3, STA303 summer 2018
# Alex Stringer
# 

library(tidyverse)
library(faraway)

# Analyze the "pvc" data from chapter 15
data(pvc)
glimpse(pvc)

# Put it in a dplyr tbl
pvc_t <- as_data_frame(pvc)

glimpse(pvc_t)

summary(pvc$psize)


# Pairwise bloxplot
pvc_t %>%
  ggplot(aes(x = operator,y = psize)) +
  theme_classic() +
  geom_boxplot() +
  labs(title = "Boxplot of Particle Size by Operator",
       x = "Operator",
       y = "Particle Size"
  )



# Summary statistics

# Grand mean and standard deviation
grand_mean <- pvc_t %>%
  summarize(grand_mean = mean(psize),grand_sd = sd(psize))
grand_mean

# Operator means and standard deviations
group_means <- pvc_t %>%
  group_by(operator) %>%
  summarize(group_mean = mean(psize),
            group_median = median(psize),
            group_sd = sd(psize),
            group_size = n())
group_means

# Sums of squares
# Here it gets a bit tricky
sums_of_squares <- pvc_t %>%
  # Add in the group mean
  left_join(group_means,by="operator") %>%
  summarize(total = sum( (psize - grand_mean$grand_mean)^2 ),
            error = sum( (psize - group_mean)^2  ),
            model = total - error
  ) %>%
  gather(type,SS,total:model)
sums_of_squares

# Create the ANOVA table
m <- pvc_t %>% pull(operator) %>% unique() %>% length()
n <- nrow(pvc_t)
sums_of_squares$df <- c(n-1,n-m,m-1)

anova_table <- sums_of_squares %>%
  mutate(MS = SS / df)
anova_table

Fstat <- (anova_table %>% filter(type == 'model') %>% pull(MS)) / (anova_table %>% filter(type == 'error') %>% pull(MS))
Fstat

critval <- qf(.95,m-1,n-m)
critval
# For a dataset of this size, we'd need to see an average between group sum-of-squares about 3.2 times as large as the average
# within-group sum-of-squares in order to conclude that the group means really were different

pval <- 1 - pf(Fstat,m-1,n-m)
pval

# ...or instead of all that, we can just use the aov function in R

pvc_anova <- aov(psize ~ operator,data = pvc_t)
summary(pvc_anova)

# Total sum of squares:
sum( (pvc_t$psize - mean(pvc_t$psize))^2 )
var(pvc_t$psize) * (nrow(pvc_t) - 1)

# Homework: calculate Df column yourself
# Then calculate Mean Sq and F value yourself


### Assumptions
### 

# Already looked at standard deviations of groups, and the boxplot
# 
# Normality:
# 
# Two random variables are equal in distribution
# if and only if their cumulative distribution functions
# are equal.
# Easier to compare their INVERSE CDFs- the so-called 
# "quantile functions"
pvc_t %>%
  mutate_at("psize",funs( (. - mean(.)) / sd(.))) %>%
  arrange(psize) %>%
  mutate(q = qnorm(1:n() / (n() + 1))) %>%
  ggplot(aes(x = q,y = psize)) +
  theme_classic() +
  geom_point() +
  geom_abline(slope = 1,intercept = 0,colour = "red") +
  labs(title = "Normal QQ-plot",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")

# Some minor deviations, but nothing alarming
# We care more about deviations in the tails than
# in the centre, for testing normality of the residuals (why?)
# - because observations in the tails are more influential to
# sample statistics (i.e. estimates of population quantites)
# 
# To get better understanding: recreate the above, but for data
# that is definitely not normally distributed
# E.g. rexp(100)
# or rpois(100)
# or even better: rt(100,df = 5)
# 

### Using linear regression
### 

# Create the model matrix
# Specify the contrasts
# 
# contr.treatment = "Effects coding" from lecture
# This is actually the default
# This happens automatically when you include an intercept in the model
# ...and including an intercept is the default.
X1 <- model.matrix(psize ~ operator,data=pvc_t,contrasts = "contr.treatment")
#View(X1)

t(X1) %*% X1 # Not orthogonal


# To get cell means coding, remove the intercept
# Add a -1 into the formula
X2 <- model.matrix(psize ~ operator -1,data=pvc_t)
#View(X2)

t(X2) %*% X2 # Orthogonal

# Define a standard linear model
# 
lm1 <- lm(psize ~ operator,data=pvc_t)
summary(lm1)

# What are those parameters estimating?
coef(lm1)[1]
coef(lm1)[1] + coef(lm1)[2]
coef(lm1)[1] + coef(lm1)[3]

group_means

# Now fit the null model, and compare them
lmnull <- lm(psize ~ 1,data=pvc_t)
summary(lmnull)

# What is that parameter estimating?
coef(lmnull)
grand_mean$grand_mean

# Compare them. VERY CONFUSING: in R, the anova() function is used to perform a likelihood ratio test
# that two models fit the data equally well. For a standard linear regression, this happens to be equivalent to 
# the ANOVA we have been discussing (and what is done using the aov() function)
# 
anova(lmnull,lm1)

# Same F statistic! This is provably true


### Two-Way ANOVA ###

# The PVC data has two factors: operator and resin
# Want to estimate group means of operator, group means of resin, and the interaction between the two
# 
# First need to check: is there data for all combinations of levels?
# Need that to be true in order for the sum of squares decomposition
# as learned in class to hold
# Also, we need replication- more than one datapoint for each combination

pvc_t %>%
  group_by(operator,resin) %>%
  summarize(count = n()) %>%
  print(n = Inf) # Prints the whole table

# What about just by resin?
pvc_t %>%
  group_by(resin) %>%
  summarize(count = n(),stddev = sd(psize)) %>%
  print(n = Inf) # Prints the whole table


# Yes, and furthermore, the design is balanced- each possible combination of factors has an equal number of levels

# Start in the same way, with boxplots, this time of both factors
boxplt1 <- pvc_t %>%
  ggplot(aes(x = operator,y = psize)) +
  theme_classic() +
  geom_boxplot() +
  labs(title = "Boxplot of Particle Size by Operator",
       x = "Operator",
       y = "Particle Size"
  )

boxplt2 <- pvc_t %>%
  ggplot(aes(x = resin,y = psize)) +
  theme_classic() +
  geom_boxplot() +
  labs(title = "Boxplot of Particle Size by Resin",
       x = "Resin",
       y = "Particle Size"
  )

# Plot these plots on a grid, together
# We use the plot_grid function in the cowplot package
# This is the only function we use from this package,
# so don't load the whole package with library()
# Just reference a single function, with cowplot::
# Still need to install with install.packages("cowplot")
cowplot::plot_grid(boxplt1,boxplt2,nrow=1)

# Unfortunately, only 6 observations per resin level
# Not enough to accurately estimate a variance
# Constant variance assumption? INCONCLUSIVE
# 


# Interaction plots
# You can do these two ways...
# 
# Note: if you're trying to plot a line in ggplot
# and you get an error about the "group aesthetic",
# try setting group=1 in aes()
plt1 <- pvc_t %>%
  group_by(operator,resin) %>%
  summarize(psize = mean(psize)) %>%
  ggplot(aes(x = operator,y = psize,group=resin)) +
  facet_grid(~resin) + 
  theme_classic() +
  geom_point() +
  geom_path() +
  labs(title = "Interaction Plot",
       subtitle = "Facetted by levels of resin",
       x = "Operator",
       y = "Particle Size")

plt1

plt2 <- pvc_t %>%
  group_by(operator,resin) %>%
  summarize(psize = mean(psize)) %>%
  ggplot(aes(x = resin,y = psize,group=operator)) +
  facet_grid(~operator) + 
  theme_classic() +
  geom_point() +
  geom_path() +
  labs(title = "Interaction Plot",
       subtitle = "Facetted by levels of operator",
       x = "Resin",
       y = "Particle Size")

plt2

cowplot::plot_grid(plt1,plt2,ncol=1)

# Does the mean of operator differ among levels of resin?
# Does the mean of resin differ among levels of operator?
# 
# For marginal mean: look at slope of line (different from zero?)
# and whether lines have same vertical location across plots
# 
# For interactions, look at shape of lines. Same across plots?
#


# Do the ANOVA using aov()
twoway_anova <- aov(psize ~ operator + resin + operator:resin,data=pvc_t)
summary(twoway_anova)

# Or, more common:
twoway_anova <- aov(psize ~ operator * resin,data=pvc_t)
summary(twoway_anova)

# Compare the sums of squares to the ANOVA without interactions
twoway_anova_mainonly <- aov(psize ~ operator + resin,data=pvc_t)
summary(twoway_anova_mainonly)

# The sums of squares are the same, but the MSE went down, which makes the F-stats for these
# effects larger, and the corresponding p-values smaller

# Equivalent regression model
lm2 <- lm(psize ~ operator + resin + operator:resin,data=pvc_t)
summary(lm2)

# What are all these parameters estimating?
group_means_int <- pvc_t %>%
  group_by(operator,resin) %>%
  summarize(psize = mean(psize))
group_means_int
# That's the mean of each group. Now, the mean should equal regression intercept plus main effect plus interaction
# For example...
# Mean of operator 1, resin 1 is 36.25
coef(lm2)['(Intercept)'] # Operator 1 and Resin 1 are the reference levels

# Mean of operator 2, resin 2 is 35.35
coef(lm2)['(Intercept)'] + coef(lm2)['operator2'] + coef(lm2)['resin2'] + coef(lm2)['operator2:resin2']

# ...and so on

# Likelihood ratio test
# We can test any hypothesis we want by specifying full model and reduced model.
# Let's first test the interaction hypothesis
# This corresponds to the gamma_ij = 0 for all ij in the linear model we wrote down
# In the equivalent linear model fit here, this corresponds to jointly testing all the operator:resin interaction
# coefficients equalling 0

lm3 <- lm(psize ~ operator + resin,data=pvc_t)
summary(lm3)

anova(lm3,lm2)

# Sum of squares, F and p all the same
#
# Let's test whether resin is useful at all- jointly test main effect of resin = 0 and interaction = 0

lm4 <- lm(psize ~ operator,data=pvc_t)
summary(lm4)

anova(lm4,lm2)


### Full Example ###
#
# The warpbreaks data (Faraway LMR, chapter 15, exercise 1) contains data on breaks in yarn woven
# by 9 looms, for each combination of levels of 2 factors: wool type (2 levels) and tension level (3 levels)
# 

# To clear environment from previous analysis:
rm(list = ls())

# Step 1: load the data

data(warpbreaks)
warpbreaks_tbl <- dplyr::as_data_frame(warpbreaks)

glimpse(warpbreaks_tbl)

# Check the observation counts in each level of the factors
warpbreaks_tbl %>%
  group_by(wool,tension) %>%
  summarize(cnt = n())

# Plot the data. We're looking for
# 1) Main effects, and equality of variance- pairwise boxplots
# 2) Interaction effects- interaction plots
# 3) Normality of the raw data

# 1)
warpbreaks_boxplot_1 <- warpbreaks_tbl %>%
  ggplot(aes(x = wool,y = breaks)) +
  theme_classic() +
  geom_boxplot() +
  labs(title = "Boxplot of breaks by wool type",
       x = "Wool Type",
       y = "Breaks")

warpbreaks_boxplot_2 <- warpbreaks_tbl %>%
  ggplot(aes(x = tension,y = breaks)) +
  theme_classic() +
  geom_boxplot() +
  labs(title = "Boxplot of breaks by level of tension",
       x = "Tension Level",
       y = "Breaks")

cowplot::plot_grid(warpbreaks_boxplot_1,warpbreaks_boxplot_2,nrow=1)

# Comments?

# 2)

warpbreaks_interactionplot_1 <- warpbreaks_tbl %>%
  group_by(wool,tension) %>%
  summarize(group_mean = mean(breaks)) %>%
  ggplot(aes(x = wool,y = group_mean,group = tension)) +
  theme_classic() +
  facet_grid(~tension) +
  geom_point() +
  geom_line() +
  labs(title = "Interaction Plot, wool x tension",
       x = "Wool",
       y = "Mean # of breaks")

warpbreaks_interactionplot_2 <- warpbreaks_tbl %>%
  group_by(wool,tension) %>%
  summarize(group_mean = mean(breaks)) %>%
  ggplot(aes(x = tension,y = group_mean,group = wool)) +
  theme_classic() +
  facet_grid(~wool) +
  geom_point() +
  geom_line() +
  labs(title = "Interaction Plot, tension x wool",
       x = "Tension",
       y = "Mean # of breaks")

cowplot::plot_grid(warpbreaks_interactionplot_1,warpbreaks_interactionplot_2)

# Do you think there is a significant interaction?
# Does the relationship between the mean # of breaks for wool A and wool B differ across different levels of tension?
# Does the relationship between the mean # of breaks between low, medium and high tensions differ for each wool type?

# 3) 
warpbreaks_tbl %>%
  mutate_at("breaks",funs( (. - mean(.)) / sd(.))) %>%
  arrange(breaks) %>%
  mutate(q = qnorm(1:n() / (n() + 1))) %>%
  ggplot(aes(x = q,y = breaks)) +
  theme_classic() +
  geom_point() +
  geom_abline(slope = 1,intercept = 0) +
  labs(title = "Normal QQ-plot",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")

# That looks problematic. What to do?
# We need to find a suitable transformation that makes the response look closer to normally distributed
# Look closer

warpbreaks_tbl %>%
  ggplot(aes(x = breaks)) +
  theme_classic() +
  geom_histogram(aes(y = ..density..),colour="black",fill="#E89A3D",bins = 15) +
  geom_density(colour = "#A200E4") +
  labs(title = "Histogram and density for observed breaks data",
       x = "Breaks",
       y = "Density")

# The problem appears to be a long right tail.
# The response is positive- use a box-cox transformation

library(MASS)

breaks_boxcox <- boxcox(breaks ~ 1,data=warpbreaks_tbl)
breaks_boxcox$x[which(breaks_boxcox$y == max(breaks_boxcox$y))]
# Indicates that lambda = -0.2222... is the best value for a power transformation, but that pretty much any value
# between around -.75 and .5 are acceptable.
# Lambda = 0 corresponds to a log transformation, so let's do that, since the result remains interpretable

warpbreaks_tbl_transform <- warpbreaks_tbl %>%
  mutate(log_breaks = log(breaks))

# If you get stupid ggplot errors like "invalid graphics state", whatever
# that means, try typing dev.off() to reset the internal graphics
# display device

warpbreaks_tbl_transform %>%
  ggplot(aes(x = log_breaks)) +
  theme_classic() +
  geom_histogram(aes(y = ..density..),colour="black",fill="#E89A3D",bins = 15) +
  geom_density(colour = "#A200E4") +
  labs(title = "Histogram and density for log-transformed observed breaks data",
       x = "Breaks",
       y = "Density")

# I mean, it looks better...
warpbreaks_tbl_transform %>%
  mutate_at("log_breaks",funs( (. - mean(.)) / sd(.))) %>%
  arrange(log_breaks) %>%
  mutate(q = qnorm(1:n() / (n() + 1))) %>%
  ggplot(aes(x = q,y = log_breaks)) +
  theme_classic() +
  geom_point() +
  geom_abline(slope = 1,intercept = 0) +
  labs(title = "Normal QQ-plot",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")

# That looks WAY better!
# But now, we have to re-do the original plots, don't we?

# 1)
warpbreaks_transformed_boxplot_1 <- warpbreaks_tbl_transform %>%
  ggplot(aes(x = wool,y = log_breaks)) +
  theme_classic() +
  geom_boxplot() +
  labs(title = "Boxplot of log(breaks) by wool type",
       x = "Wool Type",
       y = "log(Breaks)")

warpbreaks_transformed_boxplot_2 <- warpbreaks_tbl_transform %>%
  ggplot(aes(x = tension,y = log_breaks)) +
  theme_classic() +
  geom_boxplot() +
  labs(title = "Boxplot of log(breaks) by level of tension",
       x = "Tension Level",
       y = "log(Breaks)")

cowplot::plot_grid(warpbreaks_transformed_boxplot_1,warpbreaks_transformed_boxplot_2,nrow=1)

# Comments?
# The outliers disappeared, and now equality of variance looks more reasonable

# 2)

warpbreaks_transform_interactionplot_1 <- warpbreaks_tbl_transform %>%
  group_by(wool,tension) %>%
  summarize(group_mean = mean(log_breaks)) %>%
  ggplot(aes(x = wool,y = group_mean,group = tension)) +
  theme_classic() +
  facet_grid(~tension) +
  geom_point() +
  geom_line() +
  labs(title = "Interaction Plot, wool x tension",
       x = "Wool",
       y = "Mean log(# of breaks)")

warpbreaks_transform_interactionplot_2 <- warpbreaks_tbl_transform %>%
  group_by(wool,tension) %>%
  summarize(group_mean = mean(log_breaks)) %>%
  ggplot(aes(x = tension,y = group_mean,group = wool)) +
  theme_classic() +
  facet_grid(~wool) +
  geom_point() +
  geom_line() +
  labs(title = "Interaction Plot, tension x wool",
       x = "Tension",
       y = "Mean log(# of breaks)")

cowplot::plot_grid(warpbreaks_transform_interactionplot_1,warpbreaks_transform_interactionplot_2)

# Looks similar, but note the scale of the y-axis
# 
# We can now perform an ANOVA in order to test hypotheses of interest.
# There are two possible main effects, and one interaction, that we may wish to test
# We can first test the interaction, then test each of the main effects in question
# 
# First let's look at the ANOVA table

warpbreaks_aov1 <- aov(log_breaks ~ wool*tension,data=warpbreaks_tbl_transform)
summary(warpbreaks_aov1)

# We have p-values on the cusp of the arbitrary .05 significance level, so we need to be more careful than usual

warpbreaks_lm1 <- lm(log_breaks ~ wool*tension,data=warpbreaks_tbl_transform)
summary(warpbreaks_lm1)

warpbreaks_lm2 <- lm(log_breaks ~ wool + tension,data=warpbreaks_tbl_transform)
summary(warpbreaks_lm2)

warpbreaks_lm3 <- lm(log_breaks ~ wool,data=warpbreaks_tbl_transform)
summary(warpbreaks_lm3)

warpbreaks_lm4 <- lm(log_breaks ~ tension,data=warpbreaks_tbl_transform)
summary(warpbreaks_lm4)

# Test the interaction term
anova(warpbreaks_lm2,warpbreaks_lm1)

# It's on the edge. Let's jointly test whether a model with only each single main effect fits the data as well
# as a model with both main effects and their interaction
anova(warpbreaks_lm1,warpbreaks_lm3)
anova(warpbreaks_lm1,warpbreaks_lm4)


# What do we conclude?
