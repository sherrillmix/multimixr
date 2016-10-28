library('rstan')
if(parallel::detectCores()>10){
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
}

#sample type 1 for control
#sample type 2 for healthy
stanCode<-"
  data {
    int<lower=0> nSamples;
    int<lower=0> nSpecies;
    int<lower=0> nPair;
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
    real<lower=0,upper=1> pairProps[nPair,nOtu];
    for(ii in 1:nPair){
      for(jj in 1:3){
        pairProps[ii,]<-exp(otuSigmas[ii,,jj])/sum(exp(otuSigmas[ii,,jj]))
      }
    }
  }
  model {
    for (ii in 1:nPair){
      for(jj in 1:3){
        otuSigmas[ii,,jj] ~ normal(c(metaOtuMu[jj,],0),metaOtuSigma[jj,]);
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
