library('rstan')

nChain<-10
nIter<-10
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
    int<lower=0> nNegativePair;
    int<lower=0> nHealthyPair;
    int<lower=0> nSickPair;
    //dont actually need this yet
    int<lower=0> nControlSample;
    int<lower=0> nTissueSample;
    int<lower=0> nDiseaseSample;
    //could combine these three into one vector
    int<lower=1> negativeID[nSample];
    int<lower=1> healthyID[nSample];
    int<lower=1> sickID[nSample];
    int<lower=1,upper=3> sampleType[nSamples];
    int<lower=1,upper=3> samplePairType[nSamples];
    int<lower=0> counts[nSamples,nSpecies];
  }
  parameters {
    real metaOtuMu[3,nSpecies-1]; //last one fixed at 0 for identifiability
    real<lower=0> metaOtuSigma[3,nSpecies];
    real<lower=0,upper=1> otuSigmasNegative[nNegativePair,nSpecies,1];
    real<lower=0,upper=1> otuSigmasHealthy[nHealthyPair,nSpecies,2];
    real<lower=0,upper=1> otuSigmasSick[nSickPair,nSpecies,3];
    real<lower=0,upper=1> samplePropsTissue[nTissueSample,2];
    real<lower=0,upper=1> samplePropsDisease[nDiseaseSample,3];
  }
  transformed parameters{
    real<lower=0,upper=1> otuPropNegative[nNegativePair,nOtu,1];
    real<lower=0,upper=1> otuPropHealthy[nHealthyPair,nOtu,2];
    real<lower=0,upper=1> otuPropSick[nSickPair,nOtu,3];

    //convert raw proportion values to proportions
    for(ii in 1:nNegativePair) otuPropControl[ii,,1] <- exp(otuSigmasNegative[ii,,1])/sum(exp(otuSigmasNegative[ii,,1]))
    for(ii in 1:nHealthyPair){
      for(jj in 1:2) otuPropHealthy[ii,,jj] <- exp(otuSigmasHealthy[ii,,jj])/sum(exp(otuSigmasHealthy[ii,,jj]))
    }
    for(ii in 1:nSickPair){
      for(jj in 1:3) otuPropSick[ii,,jj] <- exp(otuSigmasSick[ii,,jj])/sum(exp(otuSigmasSick[ii,,jj]))
    }
  }
  model {
    #get the raw OTU proportion values for each pair
    for(ii in 1:nNegativePair) otuSigmasNegative[ii,,1] ~ normal(c(metaOtuMu[1,],0),metaOtuSigma[1,]);
    for(ii in 1:nHealthyPair){
      for(jj in 1:2) otuSigmasHealthy[ii,,jj] ~ normal(c(metaOtuMu[jj,],0),metaOtuSigma[jj,]);
    }
    for(ii in 1:nSickPair){
      for(jj in 1:3) otuSigmasSick[ii,,jj] ~ normal(c(metaOtuMu[jj,],0),metaOtuSigma[jj,]);
    }

    //set up mixing proportions
    for(ii in 1:nTissueSample)samplePropsTissue[ii,]~dirichlet(c(1,1))
    for(ii in 1:nDiseaseSample)samplePropsDisease[ii,]~dirichlet(c(1,1,1))

    for(ii in 1:nSamples){
      if(sampleType[ii]==1){
        if(samplePairType[ii]==1)counts[ii,] ~ multinomial(otuPropControl[negativeID[ii],,1]);
        if(samplePairType[ii]==2)counts[ii,] ~ multinomial(otuPropHealthy[healthyID[ii],,1]);
        if(samplePairType[ii]==3)counts[ii,] ~ multinomial(otuPropDisease[sickID[ii],,1]);
      }
      if(sampleType[ii]==2){
        //samplePairType 1 would be an error tissue sample from negative control
        if(samplePairType[ii]==2)counts[ii,] ~ multinomial(samplePropsTissue[ii,] %*% otuPropHealthy[healthyID[ii],,])
        //doesn't occur currently
        if(samplePairType[ii]==3)counts[ii,] ~ multinomial(samplePropsDisease[ii,1:2] %*% otuPropDisease[sickID[ii],,1:2]);
      }
      if(sampleType[ii]==3){
        //samplePairType 1 or 2 would be an error disease sample from negative control or healthy control
        if(samplePairType[ii]==3)counts[ii,] ~ multinomial(samplePropsDisease[ii,1:3] %*% otuPropDisease[sickID[ii],,1:3]);
      }
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
info$StudyGroup<-as.character(info$StudyGroup)
counts<-counts[,info$SampleID]


assignGroups<-function(x,selector=rep(TRUE,length(x)),outId=99999){
  cladeBs<-unique(x[selector])
  notCladeBs<-unique(x[!selector])
  cladeBIds<-structure(c(1:length(cladeBs),rep(outId,length(notCladeBs))),.Names=c(cladeBs,notCladeBs))
  return(cladeBIds)
}
negativeID<-assignGroups(info$DerivingSampleID,info$StudyGroup=='ext_ctrl')
healthyID<-assignGroups(info$DerivingSampleID,info$StudyGroup=='healthy')
sickID<-assignGroups(info$DerivingSampleID,info$StudyGroup=='sarcoid')
samplePairType<-ifelse(info$StudyGroup=='ext_ctrl',1,ifelse(info$StudyGroup=='healthy',2,3))
sampleType<-ifelse(info$SampleType!='BAL',1,ifelse(info$StudyGroup=='healthy',2,3))

dat<-list(
  nSamples=nrow(counts),  
  nSpecies=ncol(counts),
  nNegativePair=max(negativeID[negativeID<99999]),
  nHealthyPair=max(healthyID[healthyID<99999]),
  nSickPair=max(sickID[sickID<99999]),
  nControlSample=sum(sampleType==1),
  nTissueSample=sum(sampleType==2),
  nDiseaseSample=sum(sampleType==3),
  negativeID=negativeID[info$DerivingSampleID],
  healthyID=healthyID[info$DerivingSampleID],
  sickID=sickID[info$DerivingSampleID],
  sampleType=sampleType,
  samplePairType=samplePairType,
  counts=counts
)
# ,control=list(adapt_delta=.99,stepsize=.01)
# cacheOperation('work/stanFit.Rdat',
fit <- stan(model_code = stanCode, data = dat, iter=nIter, chains=nChain,thin=25)

