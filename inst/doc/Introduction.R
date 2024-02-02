## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dpi = 300
)

## ----setup--------------------------------------------------------------------
# devtools::install_github("HaoWang47/PCGII")
library(PCGII)
library(corpcor)
library(glmnet)
library(igraph)
library(Matrix)
library(mvtnorm)
library(tidyverse)
set.seed(010120000)

## ----simulate omega-----------------------------------------------------------
## Simulate precision matrix as your true network structure
N_nodes = 100
N_samples = 30
Omega = make_sf_precision_mat(e = 1, p = N_nodes)

## ----display the true network, fig.align='center', fig.cap="True Network Connections, auto-layout", fig.height=8, fig.width=8, out.width="80%"----
## Display the true network
nodenames = 1:N_nodes
links = which(lower.tri(Omega) & Omega!=0, arr.ind = TRUE)
dim(links) # display number of connections, two columns correspond to the connected nodes
my_net <- graph_from_data_frame(d=links, vertices=nodenames, directed=F) 
Ecolrs=c("gray50")
Vcolrs=c("gold")
plot(my_net, 
     edge.arrow.size=.5, 
     edge.color=Ecolrs,
     vertex.frame.color="#ffffff",
     vertex.label.color="black",
     vertex.size=3,
     layout=layout.auto(my_net)) 

## ----fig.align='center', fig.cap="True Network Connections, circle-layout", fig.height=8, fig.width=8, out.width="80%"----
plot(my_net, 
     edge.arrow.size=.5, 
     edge.color=Ecolrs,
     vertex.frame.color="#ffffff",
     vertex.label.color="black",
     vertex.size=3,
     layout=layout.circle(my_net)) 

## ----data---------------------------------------------------------------------
## Simulate Normal data
Sigma = solve(Omega)
mu = exp(qnorm(seq(from = 0.01, to = 0.99, length.out = N_nodes), mean = 2, sd=1))
norm_data = mvtnorm::rmvnorm(n = N_samples, mean = mu, sigma = Sigma)
## Convert simulated normal data to expression count data
norm_data[norm_data<0] = 0
Exp_data = round(norm_data)
head(Exp_data[,1:10])

## ----prior, fig.align='center', fig.cap="Prior Network Connections, circle-layout", fig.height=8, fig.width=8, out.width="80%"----
prior_links = links[sample(1:nrow(links), .2*nrow(links)),]
types = rep(1, 99)
types[prior_links] = 2 # 1 for unknown true connections, 2 for known prior connections
undirected_prior_links = undirected_prior(prior_links)
prior_net = graph_from_data_frame(d=cbind(links,types), vertices=nodenames, directed=F) 
Ecolrs=c("grey90", "grey10") # dark grey shows known prior connections
E(prior_net)$color = Ecolrs[E(prior_net)$types]
Vcolrs=c("gold")
plot(prior_net, 
     edge.arrow.size=.5, 
     vertex.frame.color="#ffffff",
     vertex.label.color="black",
     vertex.size=3,
     layout=layout.circle(my_net)) 

## ----pcgii--------------------------------------------------------------------
fixed_lamdba = sqrt(2*log(N_nodes/sqrt(N_samples))/N_samples)
PCGII.out = PCGII(df = Exp_data, prior = undirected_prior_links, lambda = fixed_lamdba)
PCGII.inf = inference(PCGII.out, alpha = .1)

## ----PCGII_net, fig.align='center', fig.cap="Network Recovered by PCGII, circle-layout", fig.height=8, fig.width=8, out.width="80%"----
PCGII_links = which(lower.tri(PCGII.inf) & PCGII.inf!=0, arr.ind = TRUE)
dim(PCGII_links) # display number of connections detected by PCGII
PCGII_net <- graph_from_data_frame(d=PCGII_links, vertices=nodenames, directed=F) 
Ecolrs=c("blue")
Vcolrs=c("gold")
plot(PCGII_net, 
     edge.arrow.size=.5, 
     edge.color=Ecolrs,
     vertex.frame.color="#ffffff",
     vertex.label.color="black",
     vertex.size=3,
     layout=layout.circle(PCGII_net)) 

## ----clevel-------------------------------------------------------------------
CLEVEL.out = clevel(df = Exp_data, lambda = fixed_lamdba)
CLEVEL.inf = inference(CLEVEL.out, alpha = .1)

## ----CLEVEL_net, fig.align='center', fig.cap="Network Recovered by CLEVEL, circle-layout", fig.height=8, fig.width=8, out.width="80%"----
CLEVEL_links = which(lower.tri(CLEVEL.inf) & CLEVEL.inf!=0, arr.ind = TRUE)
dim(CLEVEL_links) # display number of connections detected by PCGII
CLEVEL_net <- graph_from_data_frame(d=CLEVEL_links, vertices=nodenames, directed=F) 
Ecolrs=c("blue")
Vcolrs=c("gold")
plot(CLEVEL_net, 
     edge.arrow.size=.5, 
     edge.color=Ecolrs,
     vertex.frame.color="#ffffff",
     vertex.label.color="black",
     vertex.size=3,
     layout=layout.circle(CLEVEL_net)) 

