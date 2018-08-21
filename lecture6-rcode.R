# Supplementary R code for lecture 6, STA303 summer 2018
# Alex Stringer
# 

library(tidyverse)
library(faraway)

# Create the particle data
wafer <- data_frame(
  y = c(320,14,80,36),
  particle = rep(c("no","yes"),2),
  quality = c(rep("good",2),rep("bad",2))
)

wafer
xtabs(y ~ quality + particle,data = wafer)

# Full data, including marginal totals
wafer_with_totals <- wafer %>%
  left_join(wafer %>%
              group_by(particle) %>%
              summarize(col_total = sum(y)),
            by = "particle") %>%
  left_join(wafer %>%
              group_by(quality) %>%
              summarize(row_total = sum(y)),
            by = "quality")
  
  

# Poisson model
pois1 <- glm(y ~ particle+quality,data = wafer,family = poisson)
summary(pois1)

# Null deviance high
# Residual deviance also high, but much better
# Corresponds to hypothesis test of independence... why?
# So are particle and quality independent?
# 
# Deviance?

2 * sum(wafer$y * log(wafer$y / predict(pois1,type="response")))

# How do the model coefficients relate to the marginals?

model.matrix(pois1)
t(model.matrix(pois1)) %*% cbind(wafer %>% pull(y))

t(model.matrix(pois1)) %*% cbind(predict(pois1,type="response"))

# What happens if we put in an interaction?
pois2 <- glm(y ~ particle*quality,data = wafer,family = poisson)
summary(pois2)

t(model.matrix(pois2)) %*% cbind(predict(pois2,type="response"))
# and:
t(model.matrix(pois2)) %*% cbind(wafer %>% pull(y))
# ...but that's just because:
cbind(predict(pois2,type="response"))
cbind(wafer %>% pull(y))
# i.e:
qr(model.matrix(pois2))$rank
# Rank == dimension, i.e "full rank", i.e. has no null-space. So X^t(y - mu) = 0 implies that y = mu.

# Multinomial probabilities

mult_prob_wafer <- wafer_with_totals %>%
  mutate(unrestricted_mle = y / sum(y),
         restricted_mle = (col_total / sum(y)) * (row_total / sum(y)),
         expected_y = restricted_mle * sum(y)
  )
  
mult_prob_wafer

# LRT
mult_prob_wafer %>%
  summarize(lrt = 2*sum(y * log(sum(y)*y / (col_total * row_total))),
            lrt_a_different_way = 2*sum(y*log(unrestricted_mle / restricted_mle)))

# Compare to a chisq on (2 - 1) * (2 - 1) df- very large
# Exactly the same as testing the deviance in the poisson model
# 

# Larger two-way table example

data(haireye)
haireye_tbl <- as_data_frame(haireye)

haireye_tbl

# The inference methods are all the same as when there were only two levels per class
# The difference is from a data-analytic perspective. We can't just look at the table and get ideas about independence

xtabs(y ~ eye + hair,data = haireye_tbl)

# Let's look at some visualizations

# Basic bar chart
# May have ripped some of the below from Stackoverflow :0
haireye_tbl %>%
  ggplot(aes(x = hair, y = y, fill = hair)) + 
  theme_classic() +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = y), vjust = -0.25) +
  facet_grid(~eye, switch = "x", scales = "free_x", space = "free_x") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank()) +
  labs(title = "Hair vs Eye Colour",
       x = "Hair (lower) and Eye (upper) Colours",
       y = "# People") +
  guides(fill = FALSE) # Remove legend

# That's a start. But there are a different number of people in each group.
# Replace with proportions?

eyepropn_plt <- haireye_tbl %>%
  left_join(haireye_tbl %>%
              group_by(eye) %>%
              summarize(eye_count = sum(y)),
            by = "eye") %>%
  mutate(propn = y / eye_count) %>%
  ggplot(aes(x = hair, y = propn, fill = hair)) + 
  theme_classic() +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = scales::percent(propn)), vjust = -0.25) +
  facet_grid(~eye, switch = "x", scales = "free_x", space = "free_x") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank()) +
  labs(title = "Hair vs Eye Colour",
       subtitle = "Proportion of people with each hair colour out of people with each eye colour",
       x = "Hair (lower) and Eye (upper) Colours",
       y = "% People with each eye colour") +
  guides(fill = FALSE) +
  scale_y_continuous(labels = scales::percent_format())

eyepropn_plt
# ...but now we should do it using the other conditional probability
hairpropn_plt <- haireye_tbl %>%
  left_join(haireye_tbl %>%
              group_by(hair) %>%
              summarize(hair_count = sum(y)),
            by = "hair") %>%
  mutate(propn = y / hair_count) %>%
  ggplot(aes(x = eye, y = propn, fill = eye)) + 
  theme_classic() +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = scales::percent(propn)), vjust = -0.25) +
  facet_grid(~hair, switch = "x", scales = "free_x", space = "free_x") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank()) +
  labs(title = "Hair vs Eye Colour",
       subtitle = "Proportion of people with each eye colour out of people with each hair colour",
       x = "Hair (upper) and Eye (lower) Colours",
       y = "% People with each hair colour") +
  guides(fill = FALSE) +
  scale_y_continuous(labels = scales::percent_format())

cowplot::plot_grid(eyepropn_plt,hairpropn_plt,nrow = 1)

# It is good to look at both
# The dot plot in Faraway is just the left bar chart with a point geom instead of a bar geom,
# and the axes flipped

haireye_tbl %>%
  ggplot(aes(x = hair, y = y)) + 
  theme_classic() +
  #geom_bar(stat = "identity", width = 1) +
  geom_point(pch = 21) +
  facet_grid(eye~., switch = "x", scales = "free_x", space = "free_x") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank()) +
  labs(title = "Hair vs Eye Colour",
       subtitle = "Proportion of people with each hair colour out of people with each eye colour",
       x = "Hair (lower) and Eye (upper) Colours",
       y = "# People with each eye colour") +
  coord_flip()

# ... but to me, this plot is impossible to read.
# What about the mosaic plot? That one is more annoying.
# We can make something similar though with circles that don't touch
haireye_tbl %>%
  group_by(hair,eye) %>%
  summarize(cnt = sum(y)) %>%
  ggplot(aes(x = hair,y = eye,size = cnt)) +
  theme_classic() +
  geom_point() +
  labs(title = "Attempted size plot for hair/eye data",
       subtitle = "Size of the circles are proportional to counts for each category",
       x = "Hair Colour",
       y = "Eye Colour")

# Test independence between hair and eye colour using two methods
# Do we expect these variables to be independent, based on our initial examination
# of the data?

hairpois1 <- glm(y ~ eye + hair,data = haireye_tbl,family = poisson)
summary(hairpois1)

# You can look at the output and test the hypothesis of independence
# Reminder: what do the coefficients represent?

haireye_tbl %>%
  dplyr::select(observed = y) %>%
  mutate(expected = exp(model.matrix(hairpois1) %*% cbind(coef(hairpois1))),
         expected_another_way = predict(hairpois1,type = "response"))

# Compute the multinomial probabilities
haireye_with_totals <- haireye_tbl %>%
  left_join(haireye_tbl %>%
              group_by(hair) %>%
              summarize(hair_total = sum(y)),
            by = "hair") %>%
  left_join(haireye_tbl %>%
              group_by(eye) %>%
              summarize(eye_total = sum(y)),
            by = "eye")

haireye_with_totals

haireye_with_mle <- haireye_with_totals %>%
  mutate(unrestricted_mle = y / sum(y),
         restricted_mle = (hair_total / sum(y)) * (eye_total / sum(y)))

haireye_with_mle

haireye_with_mle %>%
  summarize(lrt = 2 * sum(y * log(unrestricted_mle / restricted_mle))) %>%
  pull(lrt)





# Three-way data

data(femsmoke)
femsmoke_tbl <- as_data_frame(femsmoke)

glimpse(femsmoke_tbl)

xtabs(y ~ smoker + age + dead,data=femsmoke_tbl)
xtabs(y ~ age + dead + smoker,data=femsmoke_tbl)
xtabs(y ~ dead + age + smoker,data=femsmoke_tbl)

# Harder to visualize 3-way data. What about proportions?
xtabs(y ~ smoker + age + dead,data=femsmoke_tbl) %>% prop.table() %>% round(2)
xtabs(y ~ age + dead + smoker,data=femsmoke_tbl) %>% prop.table() %>% round(2)
xtabs(y ~ dead + age + smoker,data=femsmoke_tbl) %>% prop.table() %>% round(2)

# Still hard. What if we marginalize over age?

xtabs(y ~ smoker + dead,data = femsmoke_tbl) %>% prop.table(1)
xtabs(y ~ smoker + dead,data = femsmoke_tbl) %>% prop.table(2)

# So smoking correlates with... lower chance of death?
# Be careful when marginalizing. Look at smoker vs age now:

xtabs(y ~ smoker + age,data = femsmoke_tbl) %>% prop.table(2)


# Way more of the non-smokers are old than the smokers (maybe smokers pass away younger...?)
# Older folks are more likely to pass away. Especially when you note that the response here is
# actually "death in the 20 years after initial contact"- of course the 75+ age group is
# going to have the highest incidence. Most people don't live to over 95, smoker or not.
# But, many smokers would pass away before 75, and hence there aren't many smokers over 75
# included in the study.

# Let's first examine whether all three factors are mutually independent

smoke_glm1 <- glm(y ~ smoker + age + dead,data = femsmoke_tbl,family = poisson)
summary(smoke_glm1)

# Doesn't look like it.
# 
# Joint independence: are smoking and death jointly independent of age?
smoke_glm2 <- glm(y ~ smoker * dead + age,data = femsmoke_tbl,family = poisson)
summary(smoke_glm2)

anova(smoke_glm1,smoke_glm2)

# So it's a bit better fit than the mutual independence model, but still not a good fit
# 
# Conditional independence
# Given age, are smoking and death independent?
# 
smoke_glm3 <- glm(y ~ smoker * age + dead * age,data = femsmoke_tbl,family = poisson)
summary(smoke_glm3)

anova(smoke_glm1,smoke_glm3)

# WAY better! The conditional independence model seems to fit as well as the saturated model
