 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!subroutine bridging  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine regionPrebas(siteOrder,HarvLim,minDharv,multiOut,nSites,nClimID,nLayers,nSp,maxYears,maxThin, &
		nYears,thinning,pCrobas,allSP,siteInfo, maxNlayers, &
		nThinning,fAPAR,initClearcut, ETSy,P0y, initVar,&
		weatherPRELES,DOY,pPRELES,etmodel, soilCinOut,pYasso,&
		pAWEN,weatherYasso,litterSize,soilCtotInOut, &
		defaultThin,ClCut,inDclct,inAclct,dailyPRELES,yassoRun,prebasVersion)

implicit none

integer, parameter :: nVar=46,npar=27!, nSp=3
integer, intent(in) :: nYears(nSites),nLayers(nSites),nSp(nSites),allSP
integer :: i,climID,ij,iz,ijj,ki,n
integer, intent(in) :: nSites, maxYears, maxThin,nClimID,maxNlayers,siteOrder(nSites,maxYears)
real (kind=8), intent(in) :: weatherPRELES(nClimID,maxYears,365,5),HarvLim(maxYears),minDharv
 integer, intent(in) :: DOY(365),etmodel
 real (kind=8), intent(in) :: pPRELES(30),pCrobas(npar,allSP)
 real (kind=8), intent(inout) :: siteInfo(nSites,7)
 real (kind=8), intent(in) :: thinning(nSites,maxThin,8),pAWEN(12,allSP)
 real (kind=8), intent(inout) :: dailyPRELES(nSites,(maxYears*365),3)
 real (kind=8), intent(inout) :: initClearcut(nSites,5)	!initial stand conditions after clear cut. (H,D,totBA,Hc,Ainit)
! real (kind=8), intent(in) :: pSp1(npar),pSp2(npar),pSp3(npar)!,par_common
 real (kind=8), intent(in) :: defaultThin(nSites),ClCut(nSites),yassoRun(nSites),prebasVersion(nSites)
 real (kind=8), intent(in) :: inDclct(nSites,allSP),inAclct(nSites,allSP)
! integer, intent(in) :: siteThinning(nSites)
 integer, intent(inout) :: nThinning(nSites)
 real (kind=8), intent(out) :: fAPAR(nSites,maxYears)
 real (kind=8), intent(inout) :: initVar(nSites,6,maxNlayers),P0y(nClimID,maxYears),ETSy(nClimID,maxYears)!,par_common
 real (kind=8), intent(inout) :: multiOut(nSites,maxYears,nVar,maxNlayers,2)
 real (kind=8), intent(inout) :: soilCinOut(nSites,maxYears,5,3,maxNlayers),soilCtotInOut(nSites,maxYears) !dimensions = nyears,AWENH,treeOrgans(woody,fineWoody,Foliage),species
 real (kind=8), intent(in) :: pYasso(35), weatherYasso(nClimID,maxYears,3),litterSize(nSites,3,maxNlayers) !litterSize dimensions: treeOrgans,species
 real (kind=8) :: output(1,nVar,maxNlayers,2),totBA(nSites), relBA(nSites,maxNlayers)
 real (kind=8) :: ClCutX, HarvArea,defaultThinX,maxState(nSites),check(maxYears)
 integer :: maxYearSite = 300,yearX(nSites),Ainit,sitex,ops(1)

!initialize run
yearX(:) = 0
!do i = 1,nSites
!! totBA(i) = sum(initVar(i,5,:))
! relBA(i,1:nLayers(i)) = initVar(i,5,1:nLayers(i))/sum(initVar(i,5,1:nLayers(i)))
!enddo

do ij = 1,maxYears
 HarvArea = 0.
! do i = 1,nSites
 do iz = 1,nSites
	i=siteOrder(iz,ij)
	ClCutX = ClCut(i)
	defaultThinX = defaultThin(i)

!!!check if the limit has been exceeded if yes no havest (thinning or clearcut will be performed)
	if (HarvLim(ij) > 0. .and. HarvArea >= HarvLim(ij)) then
	 ClCutX = 0.
	 defaultThinX = 0.
	endif
!!!
	climID = siteInfo(i,2)
	if(ij==yearX(i))then
	 yearX(i) = 0
	 
	 do ijj = 1,nLayers(i)
	  initVar(i,1,ijj) = multiOut(i,1,4,ijj,1)
	  initVar(i,2,ijj) = initClearcut(i,5)
	  initVar(i,3,ijj) = initClearcut(i,1)
	  initVar(i,4,ijj) = initClearcut(i,2)
	  initVar(i,5,ijj) = initClearcut(i,3) * relBA(i,ijj)!multiOut(i,int(ij-initClearcut(i,5)-2),13,ijj,1)/ totBA(i)
	  initVar(i,6,ijj) = initClearcut(i,4)
	  do ki = 1,int(initClearcut(i,5)+1)
	   multiOut(i,int(ij-initClearcut(i,5)+ki-1),7,ijj,1) = ki !#!#
	  enddo !ki
	 enddo !ijj
	endif

	if(prebasVersion(i)==0.) then
	  call prebas_v0(1,nLayers(i),nSp(i),siteInfo(i,:),pCrobas,initVar(i,:,1:nLayers(i)),&
		thinning(i,:,:),output(1,:,:,:),nThinning(i),maxYearSite,fAPAR(i,ij),initClearcut(i,:),&
		ETSy(climID,ij),P0y(climID,ij),weatherPRELES(climID,ij,:,:),DOY,pPRELES,etmodel, &
		soilCinOut(i,ij,:,:,1:nLayers(i)),pYasso,pAWEN,weatherYasso(climID,ij,:),&
		litterSize(i,:,1:nLayers(i)),soilCtotInOut(i,ij),&
		defaultThinX,ClCutX,inDclct(i,:),inAclct(i,:),dailyPRELES(i,(((ij-1)*365)+1):(ij*365),:),yassoRun(i))
	elseif(prebasVersion(i)==1.) then
	  call prebas_v1(1,nLayers(i),nSp(i),siteInfo(i,:),pCrobas,initVar(i,:,1:nLayers(i)),&
		thinning(i,:,:),output(1,:,:,:),nThinning(i),maxYearSite,fAPAR(i,ij),initClearcut(i,:),&
		ETSy(climID,ij),P0y(climID,ij),weatherPRELES(climID,ij,:,:),DOY,pPRELES,etmodel, &
		soilCinOut(i,ij,:,:,1:nLayers(i)),pYasso,pAWEN,weatherYasso(climID,ij,:),&
		litterSize(i,:,1:nLayers(i)),soilCtotInOut(i,ij),&
		defaultThinX,ClCutX,inDclct(i,:),inAclct(i,:),dailyPRELES(i,(((ij-1)*365)+1):(ij*365),:),yassoRun(i))
	endif

	if(sum(output(1,11,:,1))==0 .and. yearX(i) == 0) then
	 if((maxYears-ij)<10) then
	  Ainit = nint(6 + 2*3.5 - 0.005*ETSy(climID,ij) + 2.25)
	 else
	  Ainit = nint(6 + 2*3.5 - 0.005*(sum(ETSy(climID,(ij+1):(ij+10)))/10) + 2.25)
	 endif
	 yearX(i) = Ainit + ij + 1
	 initClearcut(i,5) = Ainit
	 if(ij==1) then
	  relBA(i,1:nLayers(i)) = initVar(i,5,1:nLayers(i))/sum(initVar(i,5,1:nLayers(i)))
	 else
	  relBA(i,1:nLayers(i)) = multiOut(i,(ij-1),13,1:nLayers(i),1)/sum(multiOut(i,(ij-1),13,1:nLayers(i),1))
	 endif
	endif

	multiOut(i,ij,:,:,:) = output(1,:,:,:)
	do ijj = 1,nLayers(i)
	  multiOut(i,ij,38,ijj,1) = sum(multiOut(i,1:ij,30,ijj,2)) + &
		sum(multiOut(i,1:ij,42,ijj,1)) + multiOut(i,ij,30,ijj,1)
	enddo !ijj

	initVar(i,1,:) = output(1,4,:,1)
	initVar(i,2,:) = output(1,7,:,1)
	initVar(i,3:6,:) = output(1,11:14,:,1)
	HarvArea = HarvArea + sum(output(1,37,:,1))
 end do !iz i

 !!! check if the haverst limit of the area has been reached otherwise clearcut the stands sorted by basal area 
 if (HarvArea < HarvLim(ij) .and. HarvLim(ij) /= 0.) then 
  n = 0
  do while(n < nSites .and. HarvArea < HarvLim(ij))
   n = n + 1
   do i = 1, nSites
	maxState(i) = maxval(multiOut(i,ij,12,:,1))
   enddo ! i
   ops = maxloc(maxState)
   siteX = int(ops(1))
   climID = siteInfo(siteX,2)

if(maxState(siteX)>minDharv) then
   !!clearcut!!
   HarvArea = HarvArea + sum(multiOut(siteX,ij,30,:,1))
   multiOut(siteX,ij,37,:,1) = multiOut(siteX,ij,37,:,1) + multiOut(siteX,ij,30,:,1)
   do ijj = 1, nLayers(siteX)
    multiOut(siteX,ij,6:nVar,ijj,2) = multiOut(siteX,ij,6:nVar,ijj,1) 
    multiOut(siteX,ij,26,ijj,1) = multiOut(siteX,ij,33,ijj,1) + multiOut(siteX,ij,26,ijj,1)
    multiOut(siteX,ij,27,ijj,1) = multiOut(siteX,ij,25,ijj,1) + multiOut(siteX,ij,27,ijj,1)
    multiOut(siteX,ij,28,ijj,1) = multiOut(siteX,ij,24,ijj,1) + multiOut(siteX,ij,28,ijj,1)
    multiOut(siteX,ij,29,ijj,1) = multiOut(siteX,ij,31,ijj,1)* 0.1 + & 
	multiOut(siteX,ij,32,ijj,1) + multiOut(siteX,ij,29,ijj,1) !0.1 takes into account of the stem residuals after clearcuts
    multiOut(siteX,ij,8:21,ijj,1) = 0.
    multiOut(siteX,ij,23:36,ijj,1) = 0 !#!#
    multiOut(siteX,ij,43:44,ijj,1) = 0
    multiOut(siteX,ij,38,ijj,1) = sum(multiOut(siteX,1:ij,30,ijj,2)) + &
		sum(multiOut(siteX,1:ij,42,ijj,1)) + multiOut(siteX,ij,30,ijj,1)
   enddo
	 if((maxYears-ij)<10) then
	  Ainit = nint(6 + 2*3.5 - 0.005*ETSy(climID,ij) + 2.25)
	 else
	  Ainit = nint(6 + 2*3.5 - 0.005*(sum(ETSy(climID,(ij+1):(ij+10)))/10) + 2.25)
	 endif
	 yearX(siteX) = Ainit + ij + 1
	 initClearcut(siteX,5) = Ainit
	 if(ij==1) then
	  relBA(siteX,1:nLayers(siteX)) = initVar(siteX,5,1:nLayers(siteX))/ & 
		sum(initVar(siteX,5,1:nLayers(siteX)))
	 else
	  relBA(siteX,1:nLayers(siteX)) = multiOut(siteX,(ij-1),13,1:nLayers(siteX),1)/ & 
		sum(multiOut(siteX,(ij-1),13,1:nLayers(siteX),1))
	 endif

  initVar(siteX,1,:) = 0. !output(1,4,:,1)
  initVar(siteX,2,:) = 0.!output(1,7,:,1)
  initVar(siteX,3:6,:) = 0.!output(1,11:14,:,1)
endif !(maxState(i)>minDharv)
  enddo
 endif !HarvArea < HarvLim .and. HarvLim /= 0.


end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

