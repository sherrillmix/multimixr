library('rstan')
nThreads<-50
if(parallel::detectCores()>10){
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
}

#sample type 1 for control
#sample type 2 for healthy
#sample type 3 for disease
stanCode<-"
  data {
    int<lower=0> nSamples;
    int<lower=0> nSpecies;
    int<lower=0> nPair;
    int<lower=0> nControl;
    int<lower=0> nHealthy;
    int<lower=0> nSarcoid;
    int<lower=0> counts[nSamples,nSpecies];
    int<lower=1,upper=3> sampleType[nSamples];
    int<lower=1> pairIds[nSamples];
  }
  parameters {
    real metaOtuMu[3,nSpecies-1];
    real<lower=0> metaOtuSigma[3,nSpecies];
    real<lower=0,upper=1> otuSigmas[nPair,nSpecies,3];
    real<lower=0,upper=1> sampleProps[nSamples,3];
  }
  transformed parameters{
    real<lower=0,upper=1> pairProps[nPair,nOtu,3];
    for(ii in 1:nPair){
      for(jj in 1:3){
        pairProps[ii,,jj]<-exp(otuSigmas[ii,,jj])/sum(exp(otuSigmas[ii,,jj]))
      }
    }
  }
  model {
    for (ii in 1:nPair){
      for(jj in 1:3){
        otuSigmas[ii,,jj] ~ normal(c(metaOtuMu[jj,],0),metaOtuSigma[jj,]); //last one 0 for identifiability
      }
    }
    for(ii in 1:nSamples){
      if(sampleType[ii]==1)sampleProps[ii,]~dirichlet(c(1,0,0))
      if(sampleType[ii]==2)sampleProps[ii,]~dirichlet(c(1,1,0))
      if(sampleType[ii]==3)sampleProps[ii,]~dirichlet(c(1,1,1))
      counts ~ multinomial(otuProps[ii, pairProps[pairIds[ii],] %*% sampleProps[ii,]);
    }
  }
"

tmp<-new.env()
load('data/bal_bact_counts.RData',tmp)
load('data/bal_bact_metadata.RData',tmp)
counts<-with(tmp,o)[['counts']]
taxa<-with(tmp,o)[['metadata']]
info<-with(tmp,s)
info$SampleID<-as.character(info$SampleID)
info$DerivingSampleID<-as.character(info$DerivingSampleID)
counts<-counts[,info$SampleID]


assignGroups<-function(x,selector=rep(TRUE,length(x)),outId=99999){
  cladeBs<-unique(x[selector])
  notCladeBs<-unique(x[!selector])
  cladeBIds<-structure(c(1:length(cladeBs),rep(outId,length(notCladeBs))),.Names=c(cladeBs,notCladeBs))
  return(cladeBIds)
}
pairID<-assignGroups(sort(info$DerivingSampleID))


dat<-list(
  nSamples=nrow(counts),  
  nSpecies=ncol(counts),
  nPair=length(pairID),
  counts=counts,
  sampleType=,
  pairIds=,
)
#,control=list(adapt_delta=.99,stepsize=.01)
fit <- cacheOperation('work/stanFit.Rdat',stan,model_code = stanCode, data = dat, iter=100000, chains=nThreads,thin=25)

