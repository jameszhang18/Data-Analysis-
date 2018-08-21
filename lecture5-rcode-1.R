# Supplementary R code for lecture 5, STA303 summer 2018
# Alex Stringer
# 

library(tidyverse)
library(faraway)
library(MASS)

# Load the Galapagos Islands data

data(gala)
gala_tbl <- gala %>%
  as_data_frame() %>%
  dplyr::select(-Endemics) # Alternative response variable; remove from analysis

glimpse(gala_tbl)

# Want to model species as a function of the other variables
# Look at its distribution
gala_tbl %>%
  ggplot(aes(x = Species)) +
  theme_classic() +
  geom_histogram(bins = 10,colour = "black",fill = "orange") +
  labs(title = "Histogram of Species",
       subtitle = "Galapagos Islands plant species",
       x = "Species",
       y = "# Islands")

# Definitely not normal
# Some high counts, but a lot of low ones
# What happens if we do a normal linear regression?

gala_lm1 <- lm(Species ~ .,data=gala_tbl)
summary(gala_lm1)

# The fit isn't even that bad
# But the diagnostics...

data_frame(x = residuals(gala_lm1)) %>%
  mutate_at("x",funs( (. - mean(.)) / sd(.))) %>%
  arrange(x) %>%
  mutate(q = qnorm(seq(1,length(residuals(gala_lm1))) / (1+length(residuals(gala_lm1))))) %>%
  ggplot(aes(x = q,y = x)) +
  theme_classic() +
  geom_point() +
  geom_abline(slope = 1,intercept = 0) +
  labs(title = "Normal QQ-Plot of Residuals",
       subtitle = "Normal linear model for Galapagos Islands plant species data",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")

# Not actually that bad. A lot of the counts are pretty large

data_frame(resid = residuals(gala_lm1),
           fitted = fitted(gala_lm1)) %>%
  mutate_at("resid",funs( (. - mean(.)) / sd(.))) %>%
  ggplot(aes(x = fitted,y = resid)) +
  theme_classic() +
  geom_point() +
  #geom_smooth() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -2,linetype="dashed",colour="red") +
  geom_hline(yintercept = 2,linetype="dashed",colour="red") +
  labs(title = "Residuals vs Fitted Values",
       subtitle = "Normal linear model for Galapagos Islands plant species data",
       x = "Fitted Values",
       y = "Residuals")

# Do you see the problem? This is classic
# Transformation? We have the variance stabilizing transformation. Let's try the boxcox one too and compare them
# 
gala_boxcox <- boxcox(gala_lm1)

# 0.5 is within the range of "good" values. So, let's try a square root transformation
# Boxcox actually should be (y^lambda - 1)/lambda
# Try "variance-stabilizing transformation":
gala_lm2 <- lm(sqrt(Species) ~ .,data=gala_tbl)
summary(gala_lm2)

# Fit is reasonable. Diagnostics?

data_frame(x = residuals(gala_lm2)) %>%
  mutate_at("x",funs( (. - mean(.)) / sd(.))) %>%
  arrange(x) %>%
  mutate(q = qnorm(seq(1,length(residuals(gala_lm2))) / (1+length(residuals(gala_lm2))))) %>%
  ggplot(aes(x = q,y = x)) +
  theme_classic() +
  geom_point() +
  geom_abline(slope = 1,intercept = 0) +
  labs(title = "Normal QQ-Plot of Residuals",
       subtitle = "Normal linear model for Galapagos Islands plant species data with square-root transformation",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")


data_frame(resid = residuals(gala_lm2),
           fitted = fitted(gala_lm2)) %>%
  mutate_at("resid",funs( (. - mean(.)) / sd(.))) %>%
  ggplot(aes(x = fitted,y = resid)) +
  theme_classic() +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -2,linetype="dashed",colour="red") +
  geom_hline(yintercept = 2,linetype="dashed",colour="red") +
  labs(title = "Residuals vs Fitted Values",
       subtitle = "Normal linear model for Galapagos Islands plant species data with square-root transformation",
       x = "Fitted Values",
       y = "Residuals")

# Much better.
# 
# Before we do any modelling, really, we should at least take a look at the relationship
# between the response and each predictor

plot_gala <- function(x) {
  gala_tbl %>%
    mutate(Species = log(Species)) %>%
    ggplot(aes_string(x = x,y = "Species")) +
    theme_classic() +
    geom_point() +
    labs(subtitle = "Galapagos Island Species Data",
         title = stringr::str_c("Species vs ",x),
         x = x,
         y = "log(Species)") +
    scale_x_continuous(labels = scales::comma_format())
  
}

cowplot::plot_grid(
  plot_gala("Area"),
  plot_gala("Elevation"),
  plot_gala("Nearest"),
  plot_gala("Scruz"),
  plot_gala("Adjacent"),
  ncol = 3
)

# There are several variables that clearly have nonlinear relationships with log(Species)
# Let's investigate the Poisson GLM

### Poisson GLM ###

gala_glm1 <- glm(Species ~ .,data=gala_tbl,family=poisson)
summary(gala_glm1)

# Coefficients are all small, but remember their effect is multiplicative now, not additive
# That residual deviance is a bit ridiculous. Is the fit really that bad?
# It could be caused by an outlier. But note that the null deviance is also huge.
# What proportion of the null deviance is explained by the model?
(gala_glm1$null.deviance - deviance(gala_glm1)) / gala_glm1$null.deviance

# So the model does explain much of the deviance
# 
# While in the normal linear model we assumed constant variance, 
# here was are assuming that Var(Y) = E(Y)
# This is a model assumption, and needs to be checked
# Note that the graph in Faraway is labelled incorrectly- it's a log-log plot

data_frame(muhat = predict(gala_glm1,type="response"),
           varhat = (gala_tbl$Species - muhat)^2) %>%
  mutate_all(log) %>%
  ggplot(aes(x = muhat,y = varhat)) +
  theme_classic() +
  geom_point() +
  geom_abline(slope = 1,intercept = 0,colour = "red") +
  labs(title = "Mean-Variance Relationship",
       subtitle = "Poisson GLM, Galapagos Islands species data",
       x = "Predicted Mean",
       y = "Estimated Variance at Predicted Mean")

# It looks like we have overdispersion
# 
# Estimate dispersion
n <- nrow(gala_tbl)
p <- length(coef(gala_glm1))
phi <- data_frame(muhat = predict(gala_glm1,type="response"),
                  varhat = (gala_tbl$Species - muhat)^2) %>%
  summarize(phi = sum(varhat / muhat) / (n - p)) %>%
  pull(phi)

phi

# Yes, lots of overdispersion
# Compute the new standard errors

summary(gala_glm1,dispersion = phi)
summary(gala_glm1)

sqrt(diag(vcov(gala_glm1)))
sqrt(diag(vcov(gala_glm1)) * phi)

# Offset. Area is potentially an offset
# First, try a model with log(Area)
# 
gala_glm2 <- glm(Species ~ log(Area) + Elevation + Nearest + Scruz + Adjacent,
                 data = gala_tbl,
                 family = poisson)
summary(gala_glm2)

# Coefficient for log(Area) is about 0.4
# This is probably too small to consider to be "close to 1"
# Especially relative to its own standard error of 0.015
# Probably in practice, would NOT fix it to be 1
# But, let's see what happens when we do

gala_glm3 <- glm(Species ~ offset(log(Area)) + Elevation + Nearest + Scruz + Adjacent,
                 data = gala_tbl,
                 family = poisson)
summary(gala_glm3)

# Exercise: can you try transformations (e.g. log) on other predictors
# to get a better model?
# Or interactions?
# If you're super keen, try to get a significantly better model than me,
# and show me!

### Example ###

# Use the salmonella data, exercise 2 Faraway chapter 3 page 73
# Data represents number of colonies (a count) for levels of dose of quinoline
# 
# Load the data and do some initial data analysis

data(salmonella)
salmonella_tbl <- as_data_frame(salmonella)

glimpse(salmonella_tbl)

# Plot it. Technically dose is continuous, but there are only 6 unique values,
# so it would be better to do boxplots. Except for the additional fact that there are only
# 3 datapoints per level.
# 
# So let's just plot them.
salmonella_tbl %>%
  ggplot(aes(x = dose,y = colonies)) +
  theme_classic() +
  geom_point() +
  labs(title = "Salmonella Data",
       subtitle = "# of colonies for each observed dose of quinoline",
       x = "Dose",
       y = "# colonies")

# That doesn't tell us much. Dose increases exponentially; let's log it
salmonella_tbl %>%
  mutate_at("dose",funs(log(. + 1))) %>%
  ggplot(aes(x = dose,y = colonies)) +
  theme_classic() +
  geom_point() +
  labs(title = "Salmonella Data",
       subtitle = "# of colonies for each observed dose of quinoline",
       x = "log(Dose+1)",
       y = "# colonies")

# Hmmm... maybe the boxplots were a better idea
# We could treat dose as categorical in the model

salmonella_tbl %>%
  mutate_at("dose",as.factor) %>%
  ggplot(aes(x = dose,y = colonies)) +
  theme_classic() +
  geom_boxplot() +
  labs(title = "Salmonella Data",
       subtitle = "Boxplot of # of colonies for each observed dose of quinoline",
       x = "Dose",
       y = "# colonies")

# That's better!
# 
# Fit an initial glm
dose_glm1 <- glm(colonies ~ dose,data=salmonella_tbl,family = poisson)
summary(dose_glm1)

# Not a good fit. Why? There could be reasons other than overdispersion
# Let's look at the fitted model
curve_data <- data_frame(x = seq(1,1000,by=0.1),
                         y = exp(coef(dose_glm1)[1] + x*coef(dose_glm1)[2]))

salmonella_tbl %>%
  ggplot(aes(x = dose,y = colonies)) +
  theme_classic() +
  geom_point() +
  geom_line(data = curve_data,aes(x = x,y = y)) +
  labs(title = "Salmonella Data - Fitted Poisson GLM",
       subtitle = "# of colonies for each observed dose of quinoline",
       x = "Dose",
       y = "# colonies")

# It's more likely here that we got the functional form of the model wrong
# We already said that dose was increasing exponentially, so let's log it in the model
dose_glm2 <- glm(colonies ~ log(dose+1),data=salmonella_tbl,family = poisson)
summary(dose_glm2)

# Check out the residual deviance- it's not great, but it's better than it was
curve_data <- data_frame(x = seq(1,1000,by=0.1),
                         y = exp(coef(dose_glm2)[1] + log(x+1)*coef(dose_glm2)[2]))

salmonella_tbl %>%
  ggplot(aes(x = dose,y = colonies)) +
  theme_classic() +
  geom_point() +
  geom_line(data = curve_data,aes(x = x,y = y)) +
  labs(title = "Salmonella Data - Fitted Poisson GLM",
       subtitle = "# of colonies for each observed dose of quinoline",
       x = "Dose",
       y = "# colonies")

# Better. Can we do any better than this?
# From before, it looked like there might be some sort of a quadratic relationship
# between log dose and colonies
dose_glm3 <- glm(colonies ~ log(dose+1) + I(log(dose+1)^2),data=salmonella_tbl,family = poisson)
summary(dose_glm3)

# Nah, the residual deviance didn't really change
# 
# Now estimate the overdispersion from the better model

n <- nrow(salmonella_tbl)
p <- length(coef(dose_glm2))

phi <- data_frame(muhat = predict(dose_glm2,type="response"),
                  varhat = (salmonella_tbl$colonies - muhat)^2) %>%
  summarize(phi = sum(varhat / muhat) / (n - p)) %>%
  pull(phi)
phi

summary(dose_glm2,dispersion = phi)

