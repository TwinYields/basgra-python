module soil

use parameters_site
use parameters_plant
implicit none

real*8:: FO2, fPerm, Tsurf, WCL

real*8:: DRAIN, RUNOFF
real*8:: dCLITT, rCLITT, rCSOMF, Rsoil
real*8:: dCLITTrsoil, dCLITTsomf, dCSOMF, dCSOMFrsoil, dCSOMFsoms, dCSOMS
real*8:: Nemission, NemissionN2O, NemissionNO, Nfixation, Nleaching
real*8:: NLITTnmin, NLITTsomf, Nmineralisation
real*8:: dNLITT, dNSOMF, dNSOMS, NSOMFnmin, NSOMFsoms, rNLITT, rNSOMF
real*8:: fTsoil

Contains

Subroutine SoilWaterContent(Fdepth,ROOTD,WAL)
  real*8:: Fdepth,ROOTD,WAL
  if (Fdepth < ROOTD) then
    WCL = WAL*0.001 / (ROOTD-Fdepth)
  else
    WCL = 0
  end if
end Subroutine SoilWaterContent

Subroutine Physics(DAVTMP,Fdepth,ROOTD,Sdepth,WAS, Frate)
  real*8:: DAVTMP,Fdepth,ROOTD,Sdepth,WAS
  real*8:: Frate
  if (Fdepth > 0.) then
    Tsurf = DAVTMP / (1. + 10. * (Sdepth / Fdepth) )
    fPerm = 0.
  else
    Tsurf = DAVTMP * exp(-gamma*Sdepth)
    fPerm = 1.
  end if
  call Frozensoil(Fdepth,ROOTD,WAS, Frate)
end Subroutine Physics

   Subroutine FrozenSoil(Fdepth,ROOTD,WAS, Frate)
     real*8:: Fdepth,ROOTD,WAS
     real*8:: Frate
     real*8:: alpha, PFrate, WCeff
     ! Determining the amount of solid water that contributes in transportation of heat to surface 'WCeff'
     if (Fdepth > ROOTD) then
       WCeff = WCFC
     else if (Fdepth > 0.) then
       WCeff = (0.001*WAS) / Fdepth
     else
       WCeff = WCL
     end if
     ! Calculating potential frost rate 'PFrate'
     if (((Fdepth == 0.).and.(Tsurf>0.)).or.(WCeff == 0.)) then ! No soil frost present AND no frost starting
       PFrate = 0.
     else
       alpha  = LAMBDAsoil / ( RHOwater * WCeff * LatentHeat )
       PFrate = Sqrt( max(0.,Fdepth**2 - 2.*alpha*Tsurf) ) - Fdepth
     end if
     if ((PFrate >= 0.).and.(Fdepth > 0.).and.(Fdepth < ROOTD)) then
       Frate = PFrate * (0.001*WAS/Fdepth) / WCFC ! Soil frost increasing
     else if ((PFrate+Fdepth/DELT) < 0.) then
       Frate = -Fdepth / DELT                     ! Remaining soil frost thaws away
     else
       Frate = PFrate
     end if
   end Subroutine FrozenSoil

Subroutine FRDRUNIR(EVAP,Fdepth,Frate,INFIL,poolDRAIN,ROOTD,TRAN,WAL,WAS, &
                                               FREEZEL,IRRIG,THAWS)
  real*8:: EVAP,Fdepth,Frate,INFIL,poolDRAIN,ROOTD,TRAN,WAL,WAS
  real*8:: FREEZEL,IRRIG,THAWS
  real*8:: INFILTOT,WAFC,WAST
  WAFC   = 1000. * WCFC * max(0.,(ROOTD-Fdepth))                      ! (mm)
  WAST   = 1000. * WCST * max(0.,(ROOTD-Fdepth))                      ! (mm)
  INFILTOT = INFIL + poolDrain
  if (Fdepth < ROOTD) then
    FREEZEL = max(0., min( WAL/DELT + (INFILTOT - EVAP - TRAN), &
                         (Frate/(ROOTD-Fdepth))*WAL))                 ! (mm d-1)
  else
    FREEZEL = 0.                                                      ! (mm d-1)
  end if
  if ((Fdepth > 0.) .and. (Fdepth <= ROOTD)) then
    THAWS   = max(0.,min( WAS/DELT, -Frate*WAS/Fdepth ))              ! (mm d-1)
  else
    THAWS   = 0.                                                      ! (mm d-1)
  end if
  DRAIN  = max(0.,min( DRATE, (WAL-WAFC)/DELT + &
         (INFILTOT - EVAP - TRAN - FREEZEL + THAWS) ))                ! (mm d-1)
  RUNOFF = max(0.,            (WAL-WAST)/DELT + &
         (INFILTOT - EVAP - TRAN - FREEZEL + THAWS - DRAIN) )         ! (mm d-1)
  IRRIG  = IRRIGF *  (        (WAFC-WAL)/DELT - &
         (INFILTOT - EVAP - TRAN - FREEZEL + THAWS - DRAIN - RUNOFF)) ! (mm d-1)
end Subroutine FRDRUNIR

Subroutine O2status(O2,ROOTD)
  real*8:: O2,ROOTD
  FO2 = O2 / (ROOTD * FGAS * 1000./22.4)
end Subroutine O2status

Subroutine O2fluxes(O2,PERMgas,ROOTD,RplantAer, O2IN,O2OUT)
  real*8:: O2,PERMgas,ROOTD,RplantAer
  real*8:: O2IN,O2OUT
  real*8:: O2MX
  O2OUT = RplantAer * KRTOTAER * 1./12. * 1.
  O2MX  = FO2MX * ROOTD * FGAS * 1000./22.4
  O2IN  = PERMgas * ( (O2MX-O2) + O2OUT*DELT )
end Subroutine O2fluxes

Subroutine N_fert(year,doy,DAYS_FERT,NFERTV, Nfert)
  integer                  :: year,doy,i
  integer,dimension(100,2) :: DAYS_FERT
  real*8  ,dimension(100  ) :: NFERTV
  real*8                    :: Nfert
  Nfert   = 0
  do i=1,100
    if ( (year==DAYS_FERT (i,1)) .and. (doy==DAYS_FERT (i,2)) ) then
      Nfert   = NFERTV (i)
    end if
  end do
end Subroutine N_fert

Subroutine N_dep(year,doy,DAYS_NDEP,NDEPV, Ndep)
  integer                  :: year,doy,j
  integer,dimension(100,2) :: DAYS_NDEP
  real*8  ,dimension(100  ) :: NDEPV
  integer                  :: idep
  real*8                    :: NDEPV_interval,t
  real*8  ,dimension(100)   :: tNdep
  real*8                    :: Ndep
  t     = year           + (doy           -0.5)/366
  tNdep = DAYS_NDEP(:,1) + (DAYS_NDEP(:,2)-0.5)/366
  do j = 2,100
   if ( (tNdep(j-1)<t) .and. (tNdep(j)>=t) ) idep = j-1
  end do
  NDEPV_interval = NDEPV(idep+1) - NDEPV(idep)
  Ndep           = NDEPV(idep) + NDEPV_interval * (t            -tNdep(idep)) / &
                                                  (tNdep(idep+1)-tNdep(idep))
end Subroutine N_dep

Subroutine CNsoil(ROOTD,RWA,WFPS,WAL,GCR,CLITT,CSOMF,NLITT,NSOMF,NSOMS,NMIN,CSOMS)
  real*8:: CLITT, CSOMF, CSOMS, fN2O, GCR, NLITT, NMIN, NSOMF, NSOMS
  real*8:: ROOTD, RWA, WAL, WFPS
  ! Soil temperature effect
  fTsoil = exp((Tsurf-10.)*(2.*TMAXF-Tsurf-10.)/(2.*TSIGMAF**2.))
  ! C Litter
  rCLITT      = ((CLITT*RUNOFF) / ROOTD) * RRUNBULK * 0.001
  dCLITT      =  (CLITT*fTsoil) / TCLITT
  dCLITTsomf  = FLITTSOMF * dCLITT
  dCLITTrsoil = dCLITT - dCLITTsomf
  ! C SOM fast
  rCSOMF      = ((CSOMF*RUNOFF) / ROOTD) * RRUNBULK * 0.001
  dCSOMF      =  (CSOMF*fTsoil) / TCSOMF
  dCSOMFsoms  = FSOMFSOMS * dCSOMF
  dCSOMFrsoil = dCSOMF - dCSOMFSOMS
  ! C SOM slow
  dCSOMS      = (CSOMS*fTsoil) / TCSOMS
  ! Respiration
  Rsoil       = dCLITTrsoil + dCSOMFrsoil + dCSOMS
  ! N Litter
  rNLITT      = ((NLITT*RUNOFF) / ROOTD) * RRUNBULK * 0.001
  dNLITT      =  (NLITT*dCLITT) / CLITT
  NLITTsomf   = dNLITT * FLITTSOMF
  NLITTnmin   = dNLITT - NLITTsomf
  ! N SOM fast
  rNSOMF      = ((NSOMF*RUNOFF) / ROOTD) * RRUNBULK * 0.001
  dNSOMF      =  (NSOMF*dCSOMF) / CSOMF
  NSOMFsoms   = dNSOMF * FSOMFSOMS
  NSOMFnmin   = dNSOMF - NSOMFsoms
  ! N SOM slow
  dNSOMS      = (NSOMS*dCSOMS) / CSOMS
  ! N mineralisation, fixation, leaching, emission
  Nmineralisation = NLITTnmin + NSOMFnmin + dNSOMS
  Nfixation       = gCR * KNFIX
  ! Nleaching       = (NMIN*RNLEACH*DRAIN) / WAL
  if ((WAL > 0.) .and. (NMIN > 0.)) then
    Nleaching       = (NMIN*RNLEACH*DRAIN) / WAL
  else
    Nleaching       = 0.
  end if
  Nemission       =  NMIN * KNEMIT * RWA
  fN2O            = 1. / (1. + exp(-RFN2O*(WFPS-WFPS50N2O)))
  NemissionN2O    = Nemission *     fN2O
  NemissionNO     = Nemission * (1.-fN2O)
end Subroutine CNsoil

end module soil
