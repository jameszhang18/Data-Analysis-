# Supplementary R code for lectures 10 + 11, STA303 summer 2018
# Alex Stringer
# 

library(tidyverse)
library(faraway)

# Pulp data: experiment to test the paper brightness depending on a shift operator

data(pulp)
pulp <- as_data_frame(pulp)
glimpse(pulp)

# 20 observations, how many operators and obs per operator? Plus brief summary
# 5 times
pulp %>%
  group_by(operator) %>%
  summarize(cnt = n(),
            meanbright = mean(bright),
            sdbright = sd(bright))

# Single continuous response and discrete factor- pairwise boxplots
# one-way anova? => not correct (need independent)
pulp %>%
  ggplot(aes(x = operator,y = bright)) +
  theme_classic() + 
  geom_boxplot() +
  labs(title = "Brightness by Shift Operator",
       subtitle = "Pulp Data",
       x = "Operator",
       y = "Brightness")

# Wow! So the constant variance assumption doesn't look reasonable, mostly
# because of operator A (what are they even up to?)
# 
# Fixed Effects Anova
lmod_fixed <- aov(bright ~ operator,data = pulp)
summary(lmod_fixed)
coef(lmod_fixed) # Different than Faraway
# (Intercept)   operatorb       operatorc              operatord
#(y_a bar) (y_b bar - y_a bar) (y_c bar - y_a bar) (y_d bar - y_a bar)

# Grand mean
gm <- pulp %>% pull(bright) %>% mean()
gm
# Mean of each operator
pulp %>% filter(operator == 'a') %>% pull(bright) %>% mean()
pulp %>% filter(operator == 'b') %>% pull(bright) %>% mean()
pulp %>% filter(operator == 'c') %>% pull(bright) %>% mean()
pulp %>% filter(operator == 'd') %>% pull(bright) %>% mean()

coef(lmod_fixed)[1] #mean 1= beta 0
coef(lmod_fixed)[1] + coef(lmod_fixed)[2] #mean2 = beta 0 + beta1
coef(lmod_fixed)[1] + coef(lmod_fixed)[3] #mean3 = beta 0 + beta2
coef(lmod_fixed)[1] + coef(lmod_fixed)[4] #mean4 = beta 0 + beta3

# Okay then.

# The estimated error variance is 0.106
# We can get the estimated treatment variance corresponding to the random effects model
# from this ANOVA table:
(.4467 - .1062) / 5
# 5 = measurement per operator
# estimated treatment variance = (MST - MSE)/n
# Standard deviation due to operator:
(.4467 - .1062) / 5 %>% sqrt()
# Unbiased


# Now use Maximum Likelihood, fully general (doesn't assume balanced)
# REML = TRUE is restricted ML, the default
# REML = FALSE gives regular ML (biased)
library(lme4)
lmod_mixed_ml <- lmer(bright ~ 1 + (1|operator),data = pulp,REML = FALSE)
summary(lmod_mixed_ml)
# not same as the variance before, because MLE is unbiased


# Example of fitting with fixed effect of operator:
# bright ~ 1 + operator
# The 1 represents the intercept, and is included by default (you don't need
# to type it)
# To make operator a random effect, enclose it in brackets with a |
# bright ~ 1 + (1|operator)
# The first "1" is the fixed effects intercept
# The (1) inside the operator term is the RANDOM intercept

# Variability due to operator?
.04575 / (.04575 + .10625)
# About 30%
# Total variance is sigma^2_b + sigma^2
# So proportion of variance due to operator
# is sigma^2_b / (sigma^2_b + sigma^2)

# The fixed effect intercept is the grand mean
# The Residual variance is the same as that got from the MSE
# What's the deal with the operator variance? Before it was 0.0681, now it is 0.04575
# Bias!


lmod_mixed_reml <- lmer(bright ~ 1 + (1|operator),data = pulp,REML = TRUE)
summary(lmod_mixed_reml)

# There we go. Notice fixed effects estimates unchanged
# Also notice this difference is different than what Faraway gets, because the package
# has changed (dramatically) in the 12 years since the book was published

# Variability due to operator
.0681 / (.0681 + .1062)
# About 39%, much greater than in the ML model

# Rat growth data (SM page 459 example 9.18)
library(SMPracticals)
data(rat.growth)
rat <- as_data_frame(rat.growth)
glimpse(rat)

# 5 weekly measurements on 30 rats (150 obs total)
# y is weight, units unknown
# How to plot these data?

rat %>%
  ggplot(aes(x = week,y = y,group = rat)) +
  theme_classic() + 
  geom_line(alpha = 0.3) +
  geom_line(data = rat %>% group_by(week) %>% summarize(y = mean(y)),
            aes(x = week,y = y,group = 1),
            colour = "red",
            size = 1) +
  labs(title = "Growth by Week, for each rat",
       subtitle = "Rat Growth Data",
       x = "Week",
       y = "Weight (Units unknown)")

# It looks like the actual amount of growth varies across rats (random intercept)
# but the growth rate/pattern does not (no random slope)
# 
# Can we graphically investigate whether the week itself induces variation?
rat %>%
  ggplot(aes(x = rat,y = y,group = week)) +
  theme_classic() + 
  geom_line() +
  labs(title = "Rat Growth by Rat, for each week",
       subtitle = "Rat Growth Data",
       x = "Rat",
       y = "Weight (Units unknown)")

# Better: pairwise boxplots for week
rat %>%
  mutate_at("week",as.factor) %>%
  ggplot(aes(x = week,y = y)) +
  theme_classic() + 
  geom_boxplot() +
  labs(title = "Pairwise boxplots of rat weights by week",
       subtitle = "Rat Growth Data",
       x = "Rat",
       y = "Weight (Units unknown)")

# Be careful when interpreting a graph like this. The x axis is not ordered,
# so the lines don't mean anything
# What DOES mean something is that the lines all have the same shape
# We could do a stacked bar chart, but that would be messy and doesn't scale to more rats

# Model asumptions: normality of data?

rat %>%
  dplyr::select(y) %>%
  arrange(y) %>%
  mutate_at("y",funs( (. - mean(.)) / sd(.))) %>%
  mutate(q = qnorm(seq(1:nrow(rat)) / (1 + nrow(rat)))) %>%
  ggplot(aes(x = q,y = y)) +
  theme_classic() +
  geom_point() +
  geom_abline(slope = 1,intercept = 0,colour = "red") +
  labs(title = "Normal QQ Plot, Rat Growth",
       subtitle = "Evaluating normality of weight, rat growth data",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")

# Definitely some concern in the tails. We will have to proceed anyways.

# Model 1: random intercept for rat, fixed effect for week
ratmod1 <- lmer(y ~ 1+(1|rat)+week,data = rat)
summary(ratmod1)
191.86 / (191.86 + 64.29) # % of variance due to rat
# 75%!

# Sidebar: we should have considered week to be a discrete variable,
# not a continuous variable
# Because the labels don't mean anything
# But for illustration of LME models, we will proceed as-is


# Model 2: add a random slope for rat. Each rat has a different baseline weight,
# and grows as a different linear function of week
ratmod2 <- lmer(y ~ 1+(1+week|rat)+week,data = rat)
summary(ratmod2)
119.53 / (119.53 + 12.49 + 33.84) # 72%
12.49  / (119.53 + 12.49 + 33.84) # 7%

# We didn't gain much by adding in the random slope
# This agrees with what we thought based on looking at the plot
# To formally test this hypothesis is complicated (boundary problem)

# Test another model assumption: normality of the random effects

data_frame(b = ranef(ratmod1)$rat[ ,1]) %>%
  mutate_at("b",funs( (. - mean(.)) / sd(.))) %>%
  arrange(b) %>%
  mutate(q = qnorm(seq(1:nrow(ranef(ratmod1)$rat))/(1 + nrow(ranef(ratmod1)$rat)))) %>%
  ggplot(aes(x = q,y = b)) +
  theme_classic() +
  geom_point() +
  geom_abline(slope = 1,intercept = 0,colour = "red") +
  labs(title = "Normal QQ Plot, Predicted Random Intercepts",
       subtitle = "Evaluating normality of predicted random intercepts, rat growth data",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")

# Not bad!

# Illustration of lmer formula syntax, side-by-side

# Saw already: fixed intercept for week, random intercept for rat
ratmod1 <- lmer(y ~ 1+(1|rat)+as.factor(week),data = rat)
summary(ratmod1)

# Proportion of variance in growth that is attributable to rat?
194.46 / (194.46 + 51.29)
# About 79%

# Saw already: fixed intercept for week, random intercept for rat,
# and random rat x week slope
ratmod2 <- lmer(y ~ 1+(week|rat)+week,data = rat)
summary(ratmod2)

# Same model, but with intercept and slope uncorrelated.
# The previous estimate for the correlation was 18%, not that small but
# not that large
ratmod3 <- lmer(y ~ 1+(week||rat)+week,data = rat)
summary(ratmod3)

# Forcing correlation to zero inflates the variance estimates (why?)

# More repeated measures data

# Yearly income, age, and gender information for American households
# from 1968 - 1990
data(psid)
psid <- as_data_frame(psid)
glimpse(psid)

# How many people, observations per person?
psid %>% pull(person) %>% unique() %>% length() # 85 households
psid %>%
  group_by(person) %>%
  summarize(nrows = n()) %>%
  group_by(nrows) %>%
  summarize(ntimes = n()) %>%
  ggplot(aes(x = nrows,y = ntimes)) +
  theme_classic() +
  geom_bar(stat="identity",colour = "black",fill = "lightblue") +
  labs(title = "Number of times each person appears in the dataset",
       subtitle = "PSID Data",
       x = "# of People",
       y = "# of Times") +
  scale_x_continuous(breaks = 11:23) +
  scale_y_continuous(breaks = seq(3,27,3))

# Let's look at 20 randomly selected people's yearly income for the period
# Make sure we select only people who have complete data (just for the purposes
# of this plot)
psid %>%
  inner_join(psid %>%
              group_by(person) %>%
              summarise(nrows = n()) %>%
              filter(nrows == 23) %>%
              sample_n(20),
            by = "person") %>%
  ggplot(aes(x = year,y = income,colour = sex,group = person)) +
  theme_classic() +
  facet_wrap(~person) +
  geom_line() +
  labs(title = "Income, 1978 - 1990, individual people",
       subtitle = "PSID Data",
       x = "Year",
       y = "Income",
       colour = "Gender") +
  scale_y_continuous(labels = scales::dollar_format())

# Faraway indicates that one goal of analysis is to compare incomes for
# males and females
psid %>%
  ggplot(aes(x = year,y = income,group = person)) +
  theme_classic() +
  facet_grid(~sex) +
  geom_line() +
  labs(title = "Individuals' Income by Year, for each gender",
       subtitle = "PSID Data",
       x = "Year",
       y = "Income") +
  scale_y_continuous(labels = scales::dollar_format())


# Interesting...  next let's do income in 1968 vs income in 1990
# by age, for each gender.
psid %>%
  filter(year == 68 | year == 90) %>%
  mutate_at("year",as.factor) %>%
  ggplot(aes(x = age,y = income,fill = year)) +
  theme_classic() +
  facet_grid(~sex) +
  geom_bar(stat="identity",position = "dodge") +
  labs(title = "Income by age, gender, and year",
       subtitle = "PSID Data",
       x = "Age",
       y = "Income") +
  scale_y_continuous(labels = scales::dollar_format())

# Now build a model
# Person is definitely a random effect
# Gender: fixed (not a random sample from population of all genders)
# Age: also fixed, not a random sample from population of possible ages
# Year: again, fixed
# May want to try year x age, gender x age and year x gender interactions
# First, a basic model
psid_lmer1 <- lmer(income ~ year + (1|person),data = psid)
summary(psid_lmer1)

114750124 / (114750124 + 70863916) # 62% of the variability in income is person-level

# Those variances are huge! Remember the units of the data though
# Faraway suggests a log transformation of income. This seems reasonable since
# income is strictly positive; let's make a quick plot though to see what values
# are actually in the data
psid %>%
  ggplot(aes(x = income)) +
  theme_classic() + 
  geom_histogram(colour = "black",fill = "#ff9933") +
  labs(title = "Histogram of Income",
       subtitle = "PSID Data",
       x = "Income",
       y = "# People") +
  scale_x_continuous(labels = scales::dollar_format())

# Definitely a large right skew, common in strictly positive data, and expected
# based on our rough socio-economic understanding of income dynamics (which for me can
# be summarized as "some people make a lot but most of us don't")
# Log transformation? Check if there are any zeroes
psid %>%
  filter(income == 0)
# No. What is lowest income?
psid %>% pull(income) %>% summary()
# Min income is 3... take a look at the bottom 100
psid %>%
  arrange(income) %>%
  slice(1:100) %>%
  View()

# It doesn't look like an outlier; there are just some really tiny incomes in
# the data
# Try the transformation
psid %>%
  mutate_at("income",log) %>%
  ggplot(aes(x = income)) +
  theme_classic() + 
  geom_histogram(colour = "black",fill = "#ff9933") +
  labs(title = "Histogram of log(Income)",
       subtitle = "PSID Data",
       x = "log(Income)",
       y = "# People")

# Well that seems better. Check the qqplots...
psid_qq1 <- psid %>%
  arrange(income) %>%
  mutate_at("income",funs( (. - mean(.)) / sd(.))) %>%
  mutate(q = qnorm(seq(1:nrow(psid)) / (1 + nrow(psid)))) %>%
  ggplot(aes(x = q,y = income)) +
  theme_classic() +
  geom_point() + 
  geom_abline(slope = 1,intercept = 0,colour = "red") +
  labs(title = "Normal QQ-Plot for un-transformed income",
       subtitle = "PSID Data",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")

psid_qq2 <- psid %>%
  mutate_at("income",log) %>%
  arrange(income) %>%
  mutate_at("income",funs( (. - mean(.)) / sd(.))) %>%
  mutate(q = qnorm(seq(1:nrow(psid)) / (1 + nrow(psid)))) %>%
  ggplot(aes(x = q,y = income)) +
  theme_classic() +
  geom_point() + 
  geom_abline(slope = 1,intercept = 0,colour = "red") +
  labs(title = "Normal QQ-Plot for log-transformed income",
       subtitle = "PSID Data",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")

cowplot::plot_grid(psid_qq1,psid_qq2)

# Well... the second is clearly better than the first, but that long left
# tail is troubling
# Further steps would probably be to discuss whether we could remove the low-incomes
# from the dataset- not because they are values we don't like, but rather we would
# have to figure out whether it made sense to analyze them. E.g. are they self-employed,
# etc? Are they somehow likely to be different/unrepresentative
# of American households, at least those that we want to make inferences about?
# 
# Anyways...

psid_lmer2 <- lmer(log(income) ~ year + (1|person),data = psid)
summary(psid_lmer2)

.6721 / (.6721 + .5561) # 54.7%

# Here's the predicted regression line for each person:
as_data_frame(coef(psid_lmer2)$person)

# Predicted and actual
# Average income over the whole period
pred_avg_income <- as_data_frame(coef(psid_lmer2)$person) %>%
  bind_cols(psid %>% 
              group_by(person) %>% 
              summarize(income = mean(log(income)),
                        minyear = min(year),
                        maxyear = max(year))) %>%
  rename(intercept = `(Intercept)`) %>%
  mutate(predincome = intercept + year * (minyear + (maxyear - minyear)/2))
pred_avg_income

pred_avg_income %>%
  ggplot(aes(x = predincome,y = income)) +
  theme_classic() +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1,intercept = 0,colour = "orange") +
  labs(title = "Predicted vs Actual Log Income",
       subtitle = "PSID Data, Linear Mixed Effects Model 1",
       x = "Predicted",
       y = "Actual")

# This is common: LME are usually pretty fantastic for prediction, because they accurately
# harness trended data at the individual level. People's pasts are generally pretty predictive
# of their futures.

# Yearly Income- harder
pred_yearly_income <- pred_avg_income %>%
  rename(avgpredincome = predincome,
         avgincome = income,
         slopeyear = year) %>%
  left_join(psid,by = "person") %>%
  mutate(income = log(income),
         predincome = intercept + year * slopeyear)
pred_yearly_income

# Randomly select 20 people and plot their predicted vs actual
set.seed(3247)
people_range <- pred_yearly_income %>% pull(person) %>% unique()
pred_yearly_income %>%
  inner_join(data_frame(person = sample(people_range,20))) %>%
  ggplot(aes(x = year,y = income,group = person)) +
  theme_classic() +
  facet_wrap(~person) +
  geom_line(colour = "#ff9933") +
  geom_line(aes(y = predincome),colour = "#33ccff") +
  labs(title = "Predicted vs Actual Yearly Log Income",
       subtitle = "PSID Data, Linear Mixed Effects Model. Orange = Actual, Blue = Predicted",
       x = "Year",
       y = "Log(Income)")


# Try a model with additional fixed effects
psid_lmer3 <- lmer(log(income) ~ year + age + sex + (1|person),data = psid)
summary(psid_lmer3)

exp(1.109) # 3.03
# The data suggests that a man is expected to make 3 times as much as a woman
# if they are they same age and it's the same year
# Does that effect change across ages, or years?

psid_lmer4 <- lmer(log(income) ~ year * sex + age * sex + (1|person),data = psid)
summary(psid_lmer4)

# Older men can be expected to make more than older women
# The effect is decreasing, though very slightly, with time
# Note that the large correlation between the fixed effects estimates suggests that
# the actual point estimates for these coefficents is not very precise. In particular,
# look at the standard error for the sexM point estimate- even though it is 1.36, the standard
# error of 1.11 suggests anything between -0.86 and 3.58 is reasonable given the data.
# So I wouldn't put too much trust in this.
# 
# I'm not about to start making conclusions based on all this.
# The point is that you need to balance all the sources of "information": prior knowledge,
# summary tables and plots (interpreted in the context of the data collected), regression
# model output, the nature of the predictions from such a model... it's not enough to just
# fit a model and look at a p-value. And you should be critical of everything along the way.
# 
# But by now, having taken this course, you should be prepared to do all of the above :)
# 
# Anyways, that's it for STA303- thanks for a great class.

