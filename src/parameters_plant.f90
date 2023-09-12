module parameters_plant

implicit none

! Initial constants
  real*8:: LOG10CLVI, LOG10CRESI, LOG10CRTI, CSTI, LOG10LAII
  real*8::      CLVI,      CRESI,      CRTI,            LAII
  real*8:: PHENI, TILTOTI, FRTILGI, FRTILGG1I

! Initial constants, continued
  real, parameter :: CLVDI  = 0.
  real, parameter :: YIELDI = 0.
  real, parameter :: CSTUBI = 0.
  real*8           :: LT50I

! Process parameters
  real*8:: CLAIV   , COCRESMX, CSTAVM, DAYLB   , DAYLG1G2, DAYLP  , DLMXGE, FSLAMIN
  real*8:: FSMAX   , HAGERE  , K     , KLUETILG, LAICR   , LAIEFT , LAITIL, LFWIDG
  real*8:: LFWIDV  , NELLVM  , PHENCR, PHY     , RDRSCO  , RDRSMX , RDRTEM, RGENMX
  real*8:: RGRTG1G2, ROOTDM  , RRDMAX, RUBISC  , SHAPE   , SIMAX1T, SLAMAX, SLAMIN
  real*8:: TBASE   , TCRES   , TOPTGE, TRANCO  , YG
  real*8:: RDRTMIN , TVERN
  real*8:: NCSHMAX , NCR
  real*8:: RDRROOT , RDRSTUB
  real*8:: FNCGSHMIN, TCNSHMOB, TCNUPT

  real*8:: F_WALL_LV_FMIN, F_WALL_LV_MAX, F_WALL_ST_FMIN, F_WALL_ST_MAX
  real*8:: F_DIGEST_WALL_FMIN, F_DIGEST_WALL_MAX

! Process parameters, continued
  real*8           :: Dparam, Hparam, KRDRANAER, KRESPHARD, KRSR3H
  real*8           :: LDT50A, LDT50B, LT50MN, LT50MX, RATEDMX
  real*8           :: reHardRedDay
  real, parameter :: reHardRedEnd = 91. ! If LAT<0, this is changed to 91+183 in plant.f90
  real*8           :: THARDMX, TsurfDiff

end module parameters_plant
