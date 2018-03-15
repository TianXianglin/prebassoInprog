 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!subroutine bridging  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multiPrebas(multiOut,nSites,nClimID,nLayers,nSp,maxYears,maxThin, &
		nYears,thinning,pCrobas,allSP,siteInfo, maxNlayers, &
		nThinning,fAPAR,initClearcut, ETSy,P0y, initVar,&
		weatherPRELES,DOY,pPRELES,etmodel, soilCinOut,pYasso,weatherYasso,litterSize,soilCtotInOut, &
		defaultThin,ClCut,inDclct,inAclct,dailyPRELES,yassoRun)

implicit none

integer, parameter :: nVar=46,npar=27!, nSp=3
integer, intent(in) :: nYears(nSites),nLayers(nSites),nSp(nSites),allSP
integer :: i,climID
integer, intent(in) :: nSites, maxYears, maxThin,nClimID,maxNlayers
real (kind=8), intent(in) :: weatherPRELES(nClimID,maxYears,365,5)
 integer, intent(in) :: DOY(365),etmodel
 real (kind=8), intent(in) :: pPRELES(30),pCrobas(npar,allSP),siteInfo(nSites,3)
 real (kind=8), intent(in) :: thinning(nSites,maxThin,8)
 real (kind=8), intent(inout) :: dailyPRELES(nSites,(maxYears*365),3)
 real (kind=8), intent(inout) :: initClearcut(nSites,5)	!initial stand conditions after clear cut. (H,D,totBA,Hc,Ainit)
! real (kind=8), intent(in) :: pSp1(npar),pSp2(npar),pSp3(npar)!,par_common
 real (kind=8), intent(in) :: defaultThin(nSites),ClCut(nSites),yassoRun(nSites)
 real (kind=8), intent(in) :: inDclct(nSites,allSP),inAclct(nSites,allSP)
! integer, intent(in) :: siteThinning(nSites)
 integer, intent(inout) :: nThinning(nSites)
 real (kind=8), intent(out) :: fAPAR(nSites,maxYears)
 real (kind=8), intent(in) :: initVar(nSites,6,maxNlayers),P0y(nSites,maxYears),ETSy(nSites,maxYears)!,par_common
 real (kind=8), intent(inout) :: multiOut(nSites,maxYears,nVar,maxNlayers,2)
 real (kind=8), intent(inout) :: soilCinOut(nSites,maxYears,5,3,maxNlayers),soilCtotInOut(nSites,maxYears) !dimensions = nyears,AWENH,treeOrgans(woody,fineWoody,Foliage),species
 real (kind=8), intent(in) :: pYasso(35), weatherYasso(nClimID,maxYears,3),litterSize(nSites,3,maxNlayers) !litterSize dimensions: treeOrgans,species
 real (kind=8) :: output(maxYears,nVar,maxNlayers,2)
 integer :: maxYearSite = 300


do i = 1,nSites
	climID = siteInfo(i,2)
	call prebas(nYears(i),nLayers(i),nSp(i),siteInfo(i,:),pCrobas,initVar(i,:,1:nLayers(i)),&
		thinning(i,:,:),output,nThinning(i),maxYearSite,fAPAR(i,:),initClearcut(i,:),&
		ETSy(i,:),P0y(i,:),weatherPRELES(climID,:,:,:),DOY,pPRELES,etmodel, &
		soilCinOut(i,:,:,:,1:nLayers(i)),pYasso,weatherYasso(climID,:,:),&
		litterSize(i,:,1:nLayers(i)),soilCtotInOut(i,:),&
		defaultThin(i),ClCut(i),inDclct(i,:),inAclct(i,:),dailyPRELES(i,:,:),yassoRun(i))
	multiOut(i,:,:,:,:) = output
end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

