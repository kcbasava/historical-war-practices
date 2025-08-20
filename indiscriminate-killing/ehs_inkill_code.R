#libraries----
library(brms)
library(cmdstanr)
library(dagitty)
library(dplyr)
library(ggdag)
library(sjPlot)
library(mice)
library(here)
library(ggplot2)
options(scipen=999) #turn off scientific notation

here()


#DAGs----
enemy_dag <- dagify(K ~ M + E + C,
                    E ~ C,
                    M ~ C)

#backdoor paths
adjustmentSets(enemy_dag, exposure="C", outcome="K") 
#empty set, so kill ~ PC

adjustmentSets(enemy_dag, exposure="M", outcome="K") 
#social complexity, so kill ~ military + PC

adjustmentSets(enemy_dag, exposure="E", outcome="K") 
#social complexity, so kill ~ expansion + PC

impliedConditionalIndependencies(enemy_dag) 
#expansion and military organization independent of each other conditional on political complexity

####**start here to do imputation on original dataset**#### 
# read in pre-imputed data
inkill_preimp <- read.csv(here("historical-war-practices", "inkill_preimputation.csv"))
inkill_preimp<- inkill_preimp[c(1,3,4,7,8,13,14,17:23)]

#impute (default should be predictive mean matching)
inkill_mice <- mice(inkill_preimp, method="pmm")
inkill_imp <- complete(inkill_mice)
inkill_imp[c(1,7,8)]

inkill_mice$loggedEvents

# linear model sanity checks----
summary(lm(inkill_imp$pop_res1~as.factor(inkill_imp$ter_res))) #societies with higher populations more likely to have political expansion
summary(lm(pop_res1~ as.factor(mil_res), data=inkill_imp)) #and more likely to have formal militaries 
summary(lm(pop_res1 ~ pol_cent, data=inkill_imp)) #population and centralization correlated
summary(lm(pop_res1 ~ ext_hier, data=inkill_imp)) #population and hierarchy correlated
summary(lm(pol_cent ~ ext_hier, data=inkill_imp)) #hierarchy and centralization correlated
summary(glm(ter_res ~ PC1.st, family=binomial, data=inkill_imp25))
summary(glm(mil_res ~ PC1.st, family=binomial, data=inkill_imp25))

# #PCA for political complexity----
###define standardize function from rethinking----
# standardize and normalize
standardize <- function(x) {
  x <- scale(x)
  z <- as.numeric(x)
  attr(z,"scaled:center") <- attr(x,"scaled:center")
  attr(z,"scaled:scale") <- attr(x,"scaled:scale")
  return(z)
}

ikpc.list <- list(cent=standardize(inkill_imp$pol_cent), pop=standardize(inkill_imp$pop_res1), hier=standardize(inkill_imp$ext_hier))
ikpc_df <- data.frame(ikpc.list$cent, ikpc.list$pop, ikpc.list$hier, row.names=ikpc.list$society)
ikpc <- prcomp(ikpc_df)
screeplot(ikpc, type="lines")
ikpc$sdev^2 #lambda for PC1 > 1
pve.ik <- (ikpc$sdev^2)/sum(ikpc$sdev^2) #proportion of variance explained by each PC
summary(ikpc) #PC1 explains 0.88 of variance

ikdf.pcs <- data.frame(ikpc$x) #save PCAs in a dataframe
inkill_imp$PC1 <- ikdf.pcs$PC1
inkill_imp$PC1.st <- standardize(inkill_imp$PC1)

#can't have sd = 0 for the brm, so adding making se for societies w/ no error se = 0.00000001

se_01 <- ifelse (inkill_imp$kill_se == 0.0, inkill_imp$kill_se + 0.00000000000001, inkill_imp$kill_se)
inkill_imp$k_sd <- se_01 / sd(inkill_imp$kill_mean)

#add mean dates and locations----
inkill_imp$meandate <- (inkill_imp$start_date + inkill_imp$end_date)/2


#######**start here with imputed dataset**######### 
#read in imputed and prepped data file----
inkill_imp25 <- read.csv(here("historical-war-practices","inkill_imp25.csv"))

#make categorical variables
inkill_imp25$ter_res <- as.factor(inkill_imp25$ter_res) 
inkill_imp25$mil_res <- as.factor(inkill_imp25$mil_res)

#basic regresions----
summary(lm(kill_mean ~ PC1.st + as.factor(ter_res) + as.factor(mil_res), data=inkill_imp25))
summary(lm(kill_mean ~ as.factor(ter_res) + as.factor(mil_res), data=inkill_imp25))
summary(lm(kill_mean ~ as.factor(ter_res), data=inkill_imp25))
summary(lm(kill_mean ~ PC1.st, data=inkill_imp25))


#prior only model----
ik_prioronly <-  brm(data = inkill_imp25, family = gaussian,
                     kill_mean | se(k_sd) ~ PC1.st + ter_res + mil_res,
                     prior = c(prior(normal(4, 1.3), class = Intercept),
                               prior(normal(0, 1), class=b)),
                     sample_prior = "only", 
                     iter = 2000, warmup = 1000, cores = 2, chains = 2,
                     backend="cmdstanr")
summary(ik_prioronly)
ik_prioronly$prior

#intercept only model----
ik_intercept <-  brm(data = inkill_imp25, family = gaussian,
                     kill_mean | se(k_sd, sigma=TRUE) ~ 1,
                     prior = c(prior(normal(4, 1.3), class = Intercept),
                               prior(exponential(1), class=sigma)), 
                     iter = 4000, warmup = 2000, cores = 4, chains = 4,
                     backend="cmdstanr")
summary(ik_intercept)
plot(ik_intercept)

#complexity only model----
m_p1 <- 
  brm(data = inkill_imp25, family = gaussian,
      kill_mean | se(k_sd, sigma=TRUE) ~ PC1.st,
      prior = c(prior(normal(4, 1.3), class = Intercept),
                prior(normal(0, 1), class=b),
                prior(exponential(1), class=sigma)), 
      iter = 4000, warmup = 2000, cores = 4, chains = 4,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 13),
      backend="cmdstanr")
summary(m_p1) 

mcmc_plot(m2_p, variable="b_PC1.st", type="intervals")


#complexity and military----
#effect of military after accounting for complexity and complexity -> exp + effect of comp including effects through exp, independent of effects through military
m_pm <- brm(data = inkill_imp25, family = gaussian,
            kill_mean | se(k_sd, sigma=TRUE) ~ PC1.st + mil_res,
            prior = c(prior(normal(4,1.3), class = Intercept),
                      prior(normal(0,1), class=b),
                      prior(exponential(1), class=sigma)),
            iter = 4000, warmup = 2000, cores = 4, chains = 4,
            control = list( adapt_delta = 0.99,
                            max_treedepth = 13),
            backend="cmdstanr")
summary(m_pm)
mcmc_plot(m2_pm, variable=c("b_PC1.st", "b_mil_res"), type="intervals")

#expansion and complexity----
#expansion, controlling for complexity
m_ep <- brm(data = inkill_imp25, family = gaussian,
            kill_mean | se(k_sd, sigma=TRUE) ~ ter_res + PC1.st,
            prior = c(prior(normal(4,1.3), class = Intercept),
                      prior(normal(0,1), class=b),
                      prior(exponential(1), class=sigma)),
            iter = 4000, warmup = 2000, cores = 4, chains = 4,
            control = list( adapt_delta = 0.99,
                            max_treedepth = 13),
            backend="cmdstanr")
summary(m_ep) #expansion positive but CI crosses 0
mcmc_plot(m2_ep, variable=c("b_PC1.st", "b_ter_res"), type="intervals")

#all vars----
m_pem <- brm(data = inkill_imp25, family = gaussian,
             kill_mean | se(k_sd, sigma=TRUE) ~ 1 + PC1.st + ter_res + mil_res,
             prior = c(prior(normal(4,1.3), class = Intercept),
                       prior(normal(0,1), class=b),
                       prior(exponential(1), class = sigma)),
             iter = 4000, warmup = 2000, cores = 4, chains = 4,
             control = list(adapt_delta = 0.99,
                            max_treedepth = 13),
             backend="cmdstanr") 
summary(m_pem)
mcmc_plot(m2_all1, variable=c("b_PC1.st", "b_ter_res", "b_mil_res"), type="intervals")


#model comparison----
model_weights(m_p1, m_ep, m_pm, m_pem, weights="loo")
#m_pem 0.42955619 ; m_ep 0.36150256 ; m_p1 0.10912573 ; m_pm 0.09981552

loo(m_p1) #elpd_loo   -163.4 3.4
loo(m_pm) #elpd_loo   -163.5 3.5
loo(m_ep) #elpd_loo   -162.3 3.3
loo(m_pem) #elpd_loo   -162.1 3.4
loo(ik_intercept) #elpd_loo   -162.3 3.3

summary(ik_intercept)
loos <- loo(m_p1, m_pm, m_ep, m_pem, ik_intercept)
View(loos$diffs)

ch2_coefs <- tab_model(m_pem, m_ep, m_p1, m_pm, transform=NULL, pred.labels =c("Intercept", "Complexity", "Expansion", "Military"), dv.labels=c("All predictors", "Complexity + expansion", "Complexity", "Complexity + military"), file="ch2_coefs.html")
ch2_coefs

#regional autocorrelation----
m_icc <- brm(data = inkill_imp25, family = gaussian,
             kill_mean ~ PC1.st + (1|Region),
             prior = c(prior(normal(4,1.3), class = Intercept),
                       prior(normal(0,1), class=b)),
             iter = 3000, warmup = 1000, cores = 3, chains = 3,
             backend="cmdstanr")
summary(m_icc)
#PC1.st = 0.04, 95% CI = -0.42, 0.48

m_simp <- brm(data = inkill_imp25, family = gaussian,
              kill_mean ~ PC1.st,
              prior = c(prior(normal(4,1.3), class = Intercept),
                        prior(normal(0,1), class=b)),
              iter = 3000, warmup = 1000, cores = 3, chains = 3,
              backend="cmdstanr")
summary(m_simp)
#PC1.st = 0.03, 95% CI = -0.41, 0.46

tab_model(m_icc, m_simp, pred.labels =c("Intercept", "Complexity"), dv.labels=c("Model with region", "Model without region"))

#geographic distance ----
library(geosphere)
dist_mat <- as.matrix(inkill_imp25[c(6,5)])
dist_mat <- distm(dist_mat, fun=distHaversine)
rownames(dist_mat) <- inkill_imp25$society
colnames(dist_mat) <- inkill_imp25$society

diag(dist_mat) <- NA

inkill_imp25$nndist <- apply(dist_mat, 1, min, na.rm = TRUE)
inkill_imp25$nndist <- standardize(inkill_imp25$nndist)

lm_geo <- lm(kill_mean ~ PC1.st + as.factor(ter_res) + as.factor(mil_res) + nndist, data=inkill_imp25)
summary(lm_geo)

#time autocorrelation----
time_dat <- inkill_imp25
time_dat$date.st <- standardize(time_dat$meandate)

m_time <- brm(data = time_dat, family = gaussian,
              kill_mean | se(k_sd, sigma=TRUE) ~ PC1.st + ter_res + mil_res + date.st,
              prior = c(prior(normal(4,1.3), class = Intercept),
                        prior(normal(0,1), class=b),
                        prior(exponential(1), class = sigma)),
              iter = 4000, warmup = 2000, cores = 3, chains = 3,
              control = list(adapt_delta = 0.99,
                             max_treedepth = 13),
              backend="cmdstanr") 
summary(m_time)
tab_model(m_time, m_pem, pred.labels =c("Intercept", "Complexity", "Expansion", "Military", "Date"), dv.labels=c("Model with dates", "Model without dates"))

loo(m_time)


#plots----
library(ggplot2)

#ranges by society, ordered top to bottom by complexity
ordsoc <- inkill_imp25[order(inkill_imp25$PC1.st),]
ordsoc_plot <- ggplot(data=ordsoc, aes(y=1:73, x=kill_mean, xmin=(kill_mean-kill_se), xmax=(kill_mean+kill_se))) +
  geom_errorbarh(height=.5) +
  scale_y_continuous(breaks=1:nrow(ordsoc), labels=ordsoc$society_name) +
  scale_x_continuous(breaks = scales::breaks_width(1) ) +
  labs(x='Ranges of enemies killed', y = 'Societies ordered by ascending political complexity') +
  theme_minimal(base_size=10)
ordsoc_plot

#histogram of killing scores
library(tidyr)
df_long <- inkill_imp25[c(13,14)] %>%
  pivot_longer(cols = c(kill_min, kill_max),
               names_to = "variable",
               values_to = "k_value")

ggplot(inkill_imp25, aes(x=as.factor(kill_mean), fill=as.factor(kill_mean) )) + 
  geom_bar(fill="hotpink4") +
  labs(x = "kill_mean", y = "Frequency") +
  theme_minimal() +
  theme(legend.position="none")

ggplot(df_long, aes(x = k_value, fill = variable)) +
  geom_bar(position = "dodge") +
  scale_x_continuous(breaks = c(1:8) ) +
  scale_fill_manual(values=c("hotpink4","black")) +
  labs(x = "Indiscriminate killing scores", y = "Frequency", fill = "Variable") +
  theme_minimal()
