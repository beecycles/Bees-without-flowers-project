set.seed(1)
library(lme4)
library(mvtnorm)
library(ggplot2)

d = read.csv("~/Downloads/Meiners_BeeHoneydew_data.csv")

formula = as.formula(
  "Bee_Count ~ Mold * Insecticide + 
    Sugar * Paint + 
    scale(min_day) + 
    Site +
    (1|Plant_Code) + 
    (1|julDate)"
)

Honeydew <- glmer.nb(
  formula,
  data=d,
  control = glmerControl(optimizer = "bobyqa")
)


Honeydew_poisson <- glmer(
  formula,
  data=d,
  family = poisson,
  control = glmerControl(optimizer = "bobyqa")
)

anova(Honeydew, Honeydew_poisson)

summary(Honeydew, correlation = FALSE)

# Treatments in the order listed below, plus Control in treatment 3
Mold =        c(1, 1, 0, 0, 0, 0, 0)
Insecticide = c(0, 1, 0, 1, 0, 0, 0)
Sugar =       c(0, 0, 0, 0, 0, 1, 1)
Paint =       c(0, 0, 0, 0, 1, 0, 1)

newdata = cbind(
  `(Intercept)` = 1,
  Mold = Mold,
  Insecticide = Insecticide,
  Sugar = Sugar,
  Paint = Paint,
  `scale(min_day)` = 0,
  SiteB = 0,
  SiteC = 0,
  `Mold:Insecticide` = Mold * Insecticide,
  `Sugar:Paint` = Sugar * Paint
)

treat_names = c("Natural Mold", "Natural Mold + Insecticide", "Insecticide", "Black Paint", "Sugar", "Sugar + Black Paint")
treat_names = factor(treat_names, levels = rev(treat_names))

# Grab the point estimates and uncertainty values from the model
mu = fixef(Honeydew)
Sigma = as.matrix(vcov(Honeydew))

# Double-check that my columns match up
stopifnot(all(names(mu) == colnames(newdata)))

# Sample a million coefficient combinations that are consistent with the data
coef_samples = rmvnorm(1E6, mu, Sigma)

# Expected log-visits for each of the samples
linpred = t(newdata %*% t(coef_samples))

# Mean and +/- 1SE
threenum = function(x){
  structure(
    c(mean = mean(x), quantile(x, c(.025, .975))),
    names = c("mean", "lower", "upper")
  )
}

# Subtract off the log-expectation from the control for each treatment,
# then calculate means +/- 1SE
plot_data = cbind(
  treat_names = treat_names,
  as.data.frame(t(apply(exp(linpred[,-3] - linpred[,3]), 2, threenum)))
)

library(dplyr)

d %>% group_by(Treatment_Code) %>% summarize(y = mean(Bee_Count))

plot_data

f = function(x){sum(d$Bee_Count[x])}

y = c(
  `Natural Mold` = f(d$Mold & !d$Insecticide),
  `Natural Mold + Insecticide` = f(d$Mold & d$Insecticide),
  `Insecticide` = f(!d$Mold & d$Insecticide),
  `Black Paint` = f(!d$Sugar & d$Paint),
  `Sugar` = f(d$Sugar & !d$Paint),
  `Sugar + Black Paint` = f(d$Sugar & d$Paint)
)

stopifnot(all(names(y) == plot_data$treat_names))

plot_data$y = y


ggplot(plot_data, aes(y = y))+
  geom_vline(xintercept = 1, color = "#FF000070") + 
  geom_errorbarh(aes(x = mean, xmin = lower, xmax = upper, height = 5)) +
  geom_point(aes(x = mean), color = c("orange", "darkgray", "darkgray", "darkgray", "red", "red"), size = 3) + 
  theme_bw() + 
  xlab("Multiplicative factor of experimental treatment over control treatment") +
  ylab("Total bee count per treatment") +
  scale_x_continuous(
    limits = c(0, 33), 
    expand = c(0, 0),
    breaks = c(0, 1, 10, 20, 30),
    labels = c(0, 1, 10, 20, 30),
    minor_breaks = c(5, 15, 25)
  ) +
  scale_y_continuous(limits = c(0, 125), expand = c(0, 0)) + 
  geom_text(aes(x = upper, label = treat_names), hjust = -0.1, color = "#0000AA") +
  cowplot::theme_cowplot()


# P-value for "Sugar effect > Mold effect"
mean(coef_samples[,"Mold"] > coef_samples[,"Sugar"])

# P-value for "Sugar+Paint effect < Sugar effect"
mean(coef_samples[,"Sugar"] + coef_samples[,"Paint"] + coef_samples[,"Sugar:Paint"] > coef_samples[,"Sugar"])

