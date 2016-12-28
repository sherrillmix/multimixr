library('rstan')

nChain<-50
nIter<-1000
thin<-10
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
    int<lower=1> negativeID[nSamples];
    int<lower=1> healthyID[nSamples];
    int<lower=1> sickID[nSamples];
    int<lower=1,upper=3> sampleType[nSamples];
    int<lower=1,upper=3> samplePairType[nSamples];
    int<lower=0> counts[nSamples,nSpecies];
    vector[2] tissuePrior;
    vector[3] diseasePrior;
  }
  parameters {
    matrix[3,nSpecies-1] metaOtuMu; //last one fixed at 0 for identifiability
    real<lower=0> metaOtuSigma[3,nSpecies];
    matrix[nSpecies,1] otuSigmasNegative[nNegativePair];
    matrix[nSpecies,2] otuSigmasHealthy[nHealthyPair];
    matrix[nSpecies,3] otuSigmasSick[nSickPair];
    simplex[2] samplePropsTissue[nTissueSample];
    simplex[3] samplePropsDisease[nDiseaseSample];
    real sickLocation;
    real<lower=0> sickScale;
  }
  transformed parameters{
    simplex[nSpecies] otuPropNegative[nNegativePair];
    simplex[nSpecies] otuPropHealthy[nHealthyPair,2];
    simplex[nSpecies] otuPropSick[nSickPair,3];

    //convert raw proportion values to proportions
    for(ii in 1:nNegativePair) otuPropNegative[ii,] = exp(otuSigmasNegative[ii][,1]) / sum(exp(otuSigmasNegative[ii][,1]));
    for(ii in 1:nHealthyPair){
      for(jj in 1:2) otuPropHealthy[ii,jj] = exp(otuSigmasHealthy[ii][,jj]) / sum(exp(otuSigmasHealthy[ii][,jj]));
    }
    for(ii in 1:nSickPair){
      for(jj in 1:3) otuPropSick[ii,jj] = exp(otuSigmasSick[ii][,jj]) / sum(exp(otuSigmasSick[ii][,jj]));
    }
  }
  model {
    metaOtuMu[3,]~double_exponential(sickLocation,sickScale);
    for(ii in 1:3)metaOtuSigma[ii,]~gamma(1.0,0.01);
    #get the raw OTU proportion values for each pair
    for(ii in 1:nNegativePair) otuSigmasNegative[ii][,1] ~ normal(append_col(metaOtuMu[1,],0.0),metaOtuSigma[1,]);
    for(ii in 1:nHealthyPair){
      for(jj in 1:2) otuSigmasHealthy[ii][,jj] ~ normal(append_col(metaOtuMu[jj,],0.0),metaOtuSigma[jj,]);
    }
    for(ii in 1:nSickPair){
      for(jj in 1:3) otuSigmasSick[ii][,jj] ~ normal(append_col(metaOtuMu[jj,],0.0),metaOtuSigma[jj,]);
    }

    //set up mixing proportions
    for(ii in 1:nTissueSample)samplePropsTissue[ii]~dirichlet(tissuePrior);
    for(ii in 1:nDiseaseSample)samplePropsDisease[ii]~dirichlet(diseasePrior);

    for(ii in 1:nSamples){
      if(sampleType[ii]==1){
        if(samplePairType[ii]==1)counts[ii,] ~ multinomial(otuPropNegative[negativeID[ii]]);
        if(samplePairType[ii]==2)counts[ii,] ~ multinomial(otuPropHealthy[healthyID[ii],1]);
        if(samplePairType[ii]==3)counts[ii,] ~ multinomial(otuPropSick[sickID[ii],1]);
      }
      if(sampleType[ii]==2){
        //samplePairType 1 would be an error tissue sample from negative control
        if(samplePairType[ii]==2)counts[ii,] ~ multinomial(samplePropsTissue[healthyID[ii]][1] * otuPropHealthy[healthyID[ii],1] + samplePropsTissue[healthyID[ii]][2] * otuPropHealthy[healthyID[ii],2]);
        //doesn't occur currently
        if(samplePairType[ii]==3)counts[ii,] ~ multinomial(samplePropsDisease[sickID[ii]][1] * otuPropSick[sickID[ii],1] + samplePropsDisease[sickID[ii]][2] * otuPropSick[sickID[ii],2]);
      }
      if(sampleType[ii]==3){
        //samplePairType 1 or 2 would be an error disease sample from negative control or healthy control
        if(samplePairType[ii]==3)counts[ii,] ~ multinomial(samplePropsDisease[sickID[ii]][1] * otuPropSick[sickID[ii],1] + samplePropsDisease[sickID[ii]][2] * otuPropSick[sickID[ii],2] + samplePropsDisease[sickID[ii]][3] * otuPropSick[sickID[ii],3]);
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
taxaSplit<-strsplit(taxa,'; ')
taxaSplit<-do.call(rbind,lapply(taxaSplit,function(x,n)c(x,rep(NA,n-length(x))),max(sapply(taxaSplit,length))))
taxaSplit[grep('^[a-z]__$',taxaSplit)]<-NA
otuMaxProp<-apply(apply(counts,2,function(x)x/sum(x)),1,max)
minorSelect<-otuMaxProp<.01
bak<-counts
counts<-rbind(counts[!minorSelect,],'MinorOTUs'=apply(counts[minorSelect,],2,sum))
counts<-counts[order(apply(counts,1,sum),decreasing=TRUE),]


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
  nSamples=ncol(counts),  
  nSpecies=nrow(counts),
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
  counts=t(counts),
  tissuePrior=c(1,1),
  diseasePrior=c(1,1,1)
)

# ,control=list(adapt_delta=.99,stepsize=.01)
# cacheOperation('work/stanFit.Rdat',
fit <- stan(model_code = stanCode, data = dat, iter=nIter, chains=nChain,thin=thin)
#save(fit,file='tmp.Rdat')
sims<-as.array(fit)
dim(sims)<-c(prod(dim(sims)[c(1,2)]),dim(sims)[3])
colnames(sims)<-dimnames(as.array(fit))[[3]]

apply(sims[,grep('metaOtuMu\\[3,[0-9]+\\]',colnames(sims))],2,quantile,c(.025,.975))
apply(sims[,grep('metaOtuSigma\\[3,[0-9]+\\]',colnames(sims))],2,mean)
apply(sims[,grep('samplePropsDisease\\[[0-9]+,3',colnames(sims))],2,mean)
apply(sims[,grep('samplePropsTissue\\[[0-9]+,2',colnames(sims))],2,mean)
apply(sims[,grep('otuPropNegative\\[[0-9]+,118',colnames(sims))],2,mean)
otuPred<-apply(sims[,grep('metaOtuMu',colnames(sims))],2,mean)
pdf('otuPred.pdf',width=10)
plot(1,1,type='n',xlim=c(0,nrow(counts)+1),ylim=c(1e-20,1),log='y',xaxs='i',xlab='OTU',ylab='Abundance')
for(ii in 1:3){
  exps<-exp(c(otuPred[grep(sprintf('\\[%d,',ii),names(otuPred))],0))
  print(tail(exps/sum(exps)))
  points(1:(nrow(counts)),exps/sum(exps),col=ii,cex=.7)
}
abline(v=1:nrow(counts),col='#00000033')
image(1:nrow(counts),1:ncol(counts),apply(counts[,order(sampleType,samplePairType)],2,function(x)x/sum(x)),col=rev(heat.colors(100)),xlab='OTU',ylab='',yaxt='n')
axis(2,1:ncol(counts),paste(substring(info$SampleType,1,3),substring(info$StudyGroup,1,3))[order(sampleType,samplePairType)],las=1,cex=.8)
dev.off()

#print(fit,pars='otuSigmasSick')
pdf('test.pdf',width=10)
print(plot(fit)+theme_minimal(base_family='Helvetica'))
print(traceplot(fit)+theme_minimal(base_family='Helvetica'))
print(traceplot(fit,inc_warmup=TRUE)+theme_minimal(base_family='Helvetica'))
dev.off()
