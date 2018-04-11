
multiSitePrebas_v2 <- function(nYearsMS,
                    pCROBAS = pCROB,
                    pPRELES = pPREL,
                    PREBASversion = 0,
                    etmodel = 0,
                    pYASSO =pYAS,
                    pAWEN = parsAWEN,
                    siteInfo = NA,
                    multiInitVar = NA,
                    multiThin = NA,
                    multiNthin = NA,
                    multiInitClearCut = NA,
                    nLayers = NA,
                    climIDs=NA,
                    nSp=NA,
                    PAR,
                    TAir,
                    VPD,
                    Precip,
                    CO2,
                    multiP0=NA,
                    soilC = NA,
                    weatherYasso = NA,
                    litterSize = NA,
                    soilCtot = NA,
                    defaultThin = 1.,
                    ClCut = 1.,
                    HarvLim = NA,
                    minDharv = 15,
                    inDclct = NA,
                    inAclct = NA,
                    yassoRun = 0){

  maxYears <- max(nYearsMS)
  nSites <- length(nYearsMS)
  if(all(is.na(nLayers))) nLayers <- rep(3,nSites)
  if(is.na(nSp)) nSp <- rep(3,nSites)
  allSp = ncol(pCROBAS)
  if(is.na(siteInfo)){
    siteInfo = matrix(c(1,1,3,160,0,0,20),nSites,7,byrow = T) ###default values for nspecies and site type = 3
    siteInfo[,1] <- 1:nSites
  }
  if(length(climIDs)==1) climIDs <- rep(climIDs,nSites)
  siteInfo[,2] <- climIDs
  if(length(HarvLim)==1) HarvLim <- rep(HarvLim,maxYears)
  if(any(is.na(HarvLim))) HarvLim[which(is.na(HarvLim))] <- 0.

  varNam <- getVarNam()
  nVar <- length(varNam)

  # nClimID <- dim(PAR)[1]
  if(all(!is.na(climIDs))){
    nClimID <- length(unique(climIDs))
    if(!all((1:nClimID) %in% climIDs) | length(climIDs) != nSites) return("check consistency between weather inputs and climIDs")
  } else{
    nClimID <- nSites
    climIDs <- 1:nSites}


  maxNlayers <- max(nLayers)
  layerNam <- paste("layer",1:maxNlayers)
  multiOut <- array(0, dim=c(nSites,(maxYears),nVar,maxNlayers,2),
                    dimnames = list(NULL,NULL,varNam,layerNam,
                                    c("stand","thinned")))
  initClearcut = c(1.5,0.5,0.0431969,0.,NA)
  if (is.na(multiInitClearCut)) multiInitClearCut <- matrix(initClearcut,nSites,5,byrow = T)

  ###process yasso inputs if missing
  if(is.na(soilC)) soilC <- array(0,dim=c(nSites,maxYears,5,3,maxNlayers))
  if(is.na(weatherYasso)) weatherYasso <- array(0,dim=c(nClimID,maxYears,3))
  if(is.na(litterSize)) litterSize <- array(0,dim=c(nSites,3,maxNlayers))
  if(is.na(soilCtot)) soilCtot <- matrix(0,nSites,maxYears)

  if (length(defaultThin) == 1) defaultThin=as.double(rep(defaultThin,nSites))
  if (length(ClCut) == 1) ClCut=as.double(rep(ClCut,nSites))
  if (length(inDclct) == 1) inDclct=matrix(inDclct,nSites,allSp)
  if (length(inAclct) == 1) inAclct=matrix(inAclct,nSites,allSp)
  if (length(inDclct) == nSites) inDclct=matrix(inDclct,nSites,allSp)
  if (length(inAclct) == nSites) inAclct=matrix(inAclct,nSites,allSp)
  if (length(yassoRun) == 1) yassoRun=as.double(rep(yassoRun,nSites))
  if (length(PREBASversion) == 1) PREBASversion=as.double(rep(PREBASversion,nSites))

  ###process ETS
  multiETS <- matrix(NA,nClimID,maxYears)
  for(climID in 1:nClimID){
    nYearsX <- max(nYearsMS[which(climIDs==climID)])
    Temp <- TAir[climID,1:(365*nYearsX)]-5
    ETS <- pmax(0,Temp,na.rm=T)
    ETS <- matrix(ETS,365,nYearsX); ETS <- colSums(ETS)
    multiETS[climID,(1:nYearsX)] <- ETS

    xx <- min(10,nYearsX)
    Ainit = 6 + 2*3.5 - 0.005*mean(ETS[1:xx]) + 2.25 ## this is not dependent to site type? check with Annikki
    sitesClimID <- which(climIDs==climID)
    multiInitClearCut[sitesClimID,5] <- replace(multiInitClearCut[sitesClimID,5],
            which(is.na(multiInitClearCut[sitesClimID,5])),round(Ainit))
  }
    ETSthres <- 1000; ETSmean <- rowMeans(multiETS,na.rm=T)

    ####process clearcut
  for(i in 1: nSites){
    if(ClCut[i]==1 & all(is.na(inDclct[i,]))) inDclct[i,] <-
        c(ClCutD_Pine(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          ClCutD_Spruce(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          ClCutD_Birch(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]))
    if(ClCut[i]==1 & all(is.na(inAclct[i,]))) inAclct[i,] <-
        c(ClCutA_Pine(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          ClCutA_Spruce(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]),
          ClCutA_Birch(ETSmean[climIDs[i]],ETSthres,siteInfo[i,3]))
    if(any(!is.na(inDclct[i,]))) inDclct[i,is.na(inDclct[i,])] <- max(inDclct[i,],na.rm=T)
    if(all(is.na(inDclct[i,]))) inDclct[i,] <- 9999999.99
    if(any(!is.na(inAclct[i,]))) inAclct[i,is.na(inAclct[i,])] <- max(inAclct[i,],na.rm=T)
    if(all(is.na(inAclct[i,]))) inAclct[i,] <- 9999999.99
  }

  maxThin <- max(multiNthin)
  ###thinning if missing.  To improve
  if(all(is.na(multiThin))){
    multiNthin <- rep(0,nSites)
    maxThin <- 2
    multiThin <- array(0, dim=c(nSites,maxThin,8))
  }
  multiThin[is.na(multiThin)] <- -999

###PROCESS weather inputs for prebas
  multiweather <- array(-999,dim=c(nClimID,maxYears,365,5))

  ##extract weather inputs
  for(i in 1:nClimID){
    nYearsX <- max(nYearsMS[which(climIDs==i)])
    weatherPreles <- array(c(PAR[i,1:(365*nYearsX)],TAir[i,1:(365*nYearsX)],
                             VPD[i,1:(365*nYearsX)],Precip[i,1:(365*nYearsX)],
                             CO2[i,1:(365*nYearsX)]),dim=c(365,nYearsX,5))

    weatherPreles <- aperm(weatherPreles, c(2,1,3))

    multiweather[i,(1:nYearsMS[i]),,] <- weatherPreles
  }

  ### compute P0
  ###if P0 is not provided use preles to compute P0
  if(is.na(multiP0)){
    multiP0 <- matrix(NA,nClimID,maxYears)
    for(climID in 1:nClimID){
      nYearsX <- max(nYearsMS[which(climIDs==climID)])
      P0 <- PRELES(DOY=rep(1:365,nYearsX),PAR=PAR[climID,1:(365*nYearsX)],
                   TAir=TAir[climID,1:(365*nYearsX)],VPD=VPD[climID,1:(365*nYearsX)],
                   Precip=Precip[climID,1:(365*nYearsX)],CO2=rep(380,(365*nYearsX)),
                   fAPAR=rep(1,(365*nYearsX)),LOGFLAG=0,p=pPRELES)$GPP
      P0 <- matrix(P0,365,nYearsX)
      multiP0[climID,(1:nYearsX)] <- colSums(P0)
    }}

  if (all(is.na(multiInitVar))){
    multiInitVar <- array(NA,dim=c(nSites,6,maxNlayers))
      multiInitVar[,1,] <- rep(1:maxNlayers,each=nSites)
      multiInitVar[,3,] <- initClearcut[1]; multiInitVar[,4,] <- initClearcut[2]
      multiInitVar[,5,] <- initClearcut[3]/maxNlayers; multiInitVar[,6,] <- initClearcut[4]
      multiInitVar[,2,] <- matrix(multiInitClearCut[,5],nSites,maxNlayers)
  }

  siteOrder <- matrix(1:nSites,nSites,maxYears)
  siteOrder <- apply(siteOrder,2,sample,nSites)

  prebas <- .Fortran("multiPrebas_v2",
                     siteOrder = as.matrix(siteOrder),
                     HarvLim = as.double(HarvLim),
                     minDharv = as.double(minDharv),
                     multiOut = as.array(multiOut),
                     nSites = as.integer(nSites),
                     nClimID = as.integer(nClimID),
                     nLayers = as.integer(nLayers),######
                     nSp = as.integer(nSp),######
                     maxYears = as.integer(maxYears),
                     maxThin = as.integer(maxThin),
                     nYears = as.integer(nYearsMS),
                     thinning=as.array(multiThin),
                     pCROBAS = as.matrix(pCROBAS),    ####
                     allSp = as.integer(allSp),       ####
                     siteInfo = as.matrix(siteInfo),  ####
                     maxNlayers = as.integer(maxNlayers), ####
                     nThinning=as.integer(multiNthin),
                     fAPAR=matrix(0.7,nSites,maxYears),
                     initClearcut=as.matrix(multiInitClearCut),
                     ETSy=as.matrix(multiETS),
                     P0y=as.matrix(multiP0),
                     multiInitVar=as.array(multiInitVar),
                     weather=as.array(multiweather),
                     DOY= as.integer(1:365),
                     pPRELES=as.double(pPRELES),
                     etmodel=as.integer(etmodel),
                     soilC = as.array(soilC),
                     pYASSO=as.double(pYASSO),
                     pAWEN = as.matrix(pAWEN),
                     weatherYasso = as.array(weatherYasso),
                     litterSize = as.array(litterSize),
                     soilCtot = as.matrix(soilCtot),
                     defaultThin=as.double(defaultThin),
                     ClCut=as.double(ClCut),
                     inDclct=as.matrix(inDclct),
                     inAclct=as.matrix(inAclct),
                     dailyPRELES = array(-999,dim=c(nSites,(maxYears*365),3)),
                     yassoRun=as.double(yassoRun),
                     PREBASversion=as.double(PREBASversion))
  class(prebas) <- "multiPrebas"
  prebas$totHarv <- apply(prebas$multiOut[,,37,,1],2,sum)
  return(prebas)
}



