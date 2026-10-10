#AN phylogenetic analyses of headhunting 

setwd('~/trophy-taking/AN-hh')

options(scipen=999) #if you want, turn off scientific notation
library(dplyr)
library(brms)
library(cmdstanr)

#standardize function from rethinking package (McElreath)
standardize <- function(x) {
  x <- scale(x)
  z <- as.numeric(x)
  attr(z,"scaled:center") <- attr(x,"scaled:center")
  attr(z,"scaled:scale") <- attr(x,"scaled:scale")
  return(z)
}

# #data prep and imputation----
#library(mice) 
# AN_preimp <- read.csv("AN_preimp.csv")
# AN_toimp <- AN_preimp[c(2,4:7)] 
# AN_toimp$v88.polauth <- as.ordered(AN_toimp$v88.polauth)
# AN_toimp$v61.agr <- as.ordered(AN_toimp$v61.agr)
# AN_toimp$v52.commpop <- as.ordered(AN_toimp$v52.commpop)
# AN_toimp$v48.polhier <- as.ordered(AN_toimp$v48.polhier)
# 
# AN_imp <- complete(mice(AN_toimp, method='polr'))
# 
# AN_imp <- full_join(AN_preimp[c(1:3)], AN_imp, by=join_by('taxon'))
# 
# AN_imp$v88.polauth <- as.numeric(AN_imp$v88.polauth)
# AN_imp$v61.agr <- as.numeric(AN_imp$v61.agr)
# AN_imp$v52.commpop <- as.numeric(AN_imp$v52.commpop)
# AN_imp$v48.polhier <- as.numeric(AN_imp$v48.polhier)
# 
# ##PCA----
# hh.pca1 <- prcomp(data.frame(A=standardize(AN_imp$v88.polauth), C=standardize(AN_imp$v52.commpop), H=standardize(AN_imp$v48.polhier), row.names=an.imp$taxon))
# summary(hh.pca1) #1st component = 0.75
# hhpca1df <- as.data.frame(hh.pca1$x)
# AN_imp$PC1 <- -1*hhpca1df$PC1 
# 
# #phylogeny----
#library(ape)
# an.tree <- read.nexus(file="gray2009.trees")
# 
# ## comparative data----
#library(caper) 
# an.comp <- comparative.data(an.tree, AN_imp[c(2,3,4,6,8)], names.col="taxon", force.root=TRUE)
# an.compdat <- an.comp$data
# an.compdat$society <- rownames(an.compdat)
# an.compphy <- an.comp$phy
# 
# ## covariance matrix----
# covmat <- vcv(an.compphy) 
# library(Matrix)
# cormat.pd <- nearPD(covmat, corr=TRUE, base.matrix=TRUE, ensureSymmetry=TRUE) #make positive definite
# cormat.pd <- cormat.pd$mat
# rownames(cormat.pd) <- AN_imp$taxon
# colnames(cormat.pd) <- AN_imp$taxon

#start here for prepped correlation matrix and data----
## read in prepped correlation matrix----
cormat <- readRDS('cormatpd.rds')

#phylogenetic logistic regressions----
#read in final dataset for analyses
AN_imp <- read.csv('AN_hh_dataset.csv')

## agriculture----
m.ANagr <- brm(bf(headhunting ~ as.ordered(v61.agr) + (1 | gr(society, cov = cormat))),
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,1), class = b)),
               data = an.compdat,
               data2 = list(cormat = cormat),
               family = bernoulli(),
               iter = 4000, chains = 4,
               control = list(adapt_delta = 0.89),
               backend = "cmdstanr",
               cores = 4)
summary(m.ANagr)

## political levels----
m.ANpol <- brm(bf(headhunting ~ as.ordered(v48.polhier) + (1 | gr(society, cov = cormat))),
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,1), class = b)),
               data = an.compdat,
               data2 = list(cormat = cormat),
               family = bernoulli(),
               iter = 4000, chains = 4,
               control = list(adapt_delta = 0.89),
               backend = "cmdstanr",
               cores = 4)
summary(m.ANpol)

## PC1----
m.ANpc <- brm(bf(headhunting ~ PC1 + (1 | gr(society, cov = cormat))),
              prior = c(prior(normal(0,1), class = Intercept),
                        prior(normal(0,1), class = b)),
              data = an.compdat,
              data2 = list(cormat = cormat),
              family = bernoulli(),
              iter = 4000, chains = 4,
              control = list(adapt_delta = 0.89),
              backend = "cmdstanr",
              cores = 4)
summary(m.ANpc)

## PC1 + agriculture----
m.ANpcagr <- brm(bf(headhunting ~ PC1 + as.ordered(v61.agr) + (1 | gr(taxon, cov = cormat))),
                 prior = c(prior(normal(0,1), class = Intercept),
                           prior(normal(0,1), class = b)),
                 data = an.imp,
                 data2 = list(cormat = cormat),
                 family = bernoulli(),
                 iter = 4000, chains = 4,
                 control = list(adapt_delta = 0.89),
                 backend = "cmdstanr",
                 cores = 4)
summary(m.ANpcagr)

#signal ----
compphy_rooted <- root(an.compphy, node=115) #root tree
compphy_rooted1 <- multi2di(compphy_rooted) #resolve polytomies
compphy_rooted1$node.label <- NULL #remove node labels because new nodes are unnamed and therefore duplicates
compphy_rooted1$edge.length[which(compphy_rooted1$edge.length==0)]=0.000000000001 #0 length branches
hhsig <- phylo.d(an.imp, compphy_rooted1, names.col=taxon, binvar=headhunting) #and finally signal 
hhsig$Pval0 #probability of resulting from Brownian phylogenetic structure = 0.462
hhsig$Pval1 #probability of resulting from no (random) phylogenetic structure = 0

#plotting----
library(viridis)
x <- an.compdat$headhunting
y <- ifelse(an.compdat$v48.polhier >1, 1, 0)
z <-  ifelse(an.compdat$v61.agr >2, 1, 0)
names(x) <- names(y) <- names(z) <- rownames(an.compdat)
co <- c("gray87","#fde725")
co1 <- c("gray87","#21918c")
co2 <- c("gray87","#453581")

plot(an.compphy, label.offset = 0.15, cex=0.5, font=0, no.margin=TRUE)
tiplabels(pch = 22, bg = co[as.numeric(x)+1], cex = 1)
tiplabels(pch = 22, bg = co1[as.numeric(z)+1], cex = 1, offset= 0.05)
tiplabels(pch = 22, bg = co2[as.numeric(y)+1], cex = 1, offset = 0.1)

