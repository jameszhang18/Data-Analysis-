# Supplementary R code for lecture 9, STA303 summer 2018
# Alex Stringer
# 

library(tidyverse)
library(faraway)

# An experiment measuring death rates for insects, with 30 insects at each of five treatment levels.
# BE CAREFUL: make sure you are loading the right dataset from right package
# Specify package if not sure, and always check the docs.
data(bliss,package = "faraway")
bliss

# Function to implement IRWLS for binomial data
irwls_binomial_glm <- function(y,n,dat) {
  # y: observed vector of PROPORTIONS
  # n: integer, denominator of y (i.e. how many binomial trials)
  # dat: dataframe containing the covariates as columns
  
  # Set initial value for mu
  mu <- y
  
  # Do 30 iterations
  # Because that's what glm() does.
  for (i in 1:30) {
    if (i == 1) {
      # eta = g(mu)
      eta <- logit(mu)
    }
    else {
      # Update eta, then update mu
      eta <- lmod$fit
      mu <- ilogit(eta)
    }
    
    # Linearized response
    # d eta / d mu = (mu*(1-mu))^-1
    z <- eta + (y - mu) / (mu * (1-mu))
    # Get the weights
    w <- n * mu * (1-mu)
    # Regress z on X, with weights w
    lmod <- lm(z ~ .,data = dat,weights = w)
    
    cat("Iteration:",i," Estimated Coefficients: ",coef(lmod),"\n")
  }
  
  return(coef(lmod))
}

irwls_binomial_glm(y = bliss$dead/30,n = 30,dat = dplyr::select(bliss,conc))

# Compare with glm()
glm_mod <- glm(cbind(dead,alive) ~ conc,data = bliss,family = binomial)
summary(glm_mod)

glm_coef <- coef(glm_mod)
glm_coef

# What if we just used newton's method on the regular log-likelihood?
# Remember at the maximum, must satisfy normal equations
library(rootSolve)

normal_equations <- function(beta,y,X) {
  # Add an intercept to X
  X <- cbind(rep(1,nrow(X)),X)
  # The predicted values for a binomial regression:
  ypred <- apply(X,1,function(x) 1 / (1 + exp(-sum(x * beta))))
  # This part is the same for every GLM
  t(X) %*% (cbind(y - ypred))
}

# Check the glm() answer was good
normal_equations(beta = glm_coef,y = bliss$dead/30,X = cbind(bliss$conc))

# Yep!
# Now find the zero(es) of this equation numerically
rootSolve::multiroot(normal_equations,start = c(0,0),y = bliss$dead/30,X = cbind(bliss$conc))

# Though we motivated it from a regression stance, it can be shown that IRWLS is just NR for GLM
# (that's a lot of TLAs!) (TLA == "Three Letter Acronym" #meta)

# Diagnostics
# Re-do the Galapagos Islands analysis
data(gala)
gala_tbl <- gala %>%
  as_data_frame() %>%
  dplyr::select(-Endemics) # Alternative response variable; remove from analysis

glimpse(gala_tbl)

gala_glm1 <- glm(Species ~ .,data=gala_tbl,family=poisson)
summary(gala_glm1)

# Based on residual deviance of 716.85 on 24 df,
# doesn't look like it fits as well as saturated.

# Plot the deviance residuals against the fitted linear predictor
data_frame(devres = residuals(gala_glm1,type="deviance"),
           fittedeta = predict(gala_glm1,type="link")) %>%
  ggplot(aes(x = fittedeta,y = devres)) +
  theme_classic() +
  geom_point() +
  geom_hline(yintercept = 0,colour="red") +
  labs(title = "Deviance Residuals vs Fitted Linear Predictor",
       subtitle = "Galapagos Islands Species Data",
       x = "Fitted Linear Predictor",
       y = "Deviance Residuals")

# Neither of the things we're looking for - misspecification of the structural
# form of the model (nonlinearity) and misspecification of the mean/variance
# relationship - are clearly present here
# 
# 
# Investigate the relationship between the response and each predictor
# Compute the "linearized" response from the fitted model
mu <- predict(gala_glm1,type = "response")
eta <- predict(gala_glm1,type = "link")
z <- eta + (pull(gala_tbl,Species) - mu) / mu

# Recreate the plots from Faraway investigating the relationship with Area
plt1 <- gala_tbl %>%
  ggplot(aes(x = Area,y = Species)) +
  theme_classic() +
  geom_point() +
  labs(title = "Species vs Area",
       subtitle = "Galapagos Islands Species Data",
       x = "Area",
       y = "Species") +
  scale_x_continuous(labels = scales::comma_format())

plt2 <- gala_tbl %>%
  mutate(Area = log(Area)) %>%
  ggplot(aes(x = Area,y = Species)) +
  theme_classic() +
  geom_point() +
  labs(title = "Species vs log(Area)",
       subtitle = "Galapagos Islands Species Data",
       x = "log(Area)",
       y = "Species") +
  scale_x_continuous(labels = scales::comma_format())

plt3 <- data_frame(z = z,Area = gala_tbl %>% pull(Area)) %>%
  mutate(Area = log(Area)) %>%
  ggplot(aes(x = Area,y = z)) +
  theme_classic() +
  geom_point() +
  labs(title = "Linearized Response vs log(Area)",
       subtitle = "Galapagos Islands Species Data",
       x = "log(Area)",
       y = "Linearized Response") +
  scale_x_continuous(labels = scales::comma_format())

cowplot::plot_grid(plt1,plt2,plt3,nrow=1)

# This illustrates that a log transformation on Area is probably reasonable
# Notice the huge difference in plotting the linearized response instead of the regular response
# Use partial residual plots to further assess whether the relationship between the (linked) response
# and the covariates is linear, and whether there are any outlying points
beta <- coef(gala_glm1)
beta

partial_residuals <- gala_tbl %>%
  mutate(z = z,
         eta = eta) %>%
  mutate(Area_pres = z - (eta - Area * beta['Area']),
         Elevation_pres = z - eta + Elevation * beta['Elevation'],
         Nearest_pres = z - eta + Nearest * beta['Nearest'],
         Scruz_pres = z - eta + Scruz * beta['Scruz'],
         Adjacent_pres = z - eta + Adjacent * beta['Adjacent'])

# Plot them
plot_pres <- function(var) {
  partial_residuals %>%
    ggplot(aes_string(x = var,y = stringr::str_c(var,"_pres"))) +
    theme_classic() +
    geom_point() +
    labs(title = stringr::str_c("Partial Residual Plot for ",var),
         subtitle = "Galapagos Islands Species Data",
         x = var,
         y = "Partial Residuals") +
    scale_x_continuous(labels = scales::comma_format())
}

cowplot::plot_grid(
  plot_pres("Area"),
  plot_pres("Elevation"),
  plot_pres("Nearest"),
  plot_pres("Scruz"),
  plot_pres("Adjacent"),
  nrow = 2
)
