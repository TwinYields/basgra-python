module parameters_site

! Simulation period and time step
  real, parameter       :: DELT   =   1.0

! Geography
  real*8                 :: LAT

! Atmospheric conditions
  real, parameter       :: CO2A   = 350

! Soil
  real, parameter       :: DRATE  =  50
  real*8                 :: WCI
  real*8                 :: FWCAD, FWCWP, FWCFC, FWCWET, WCST
  real*8                 ::  WCAD,  WCWP,  WCFC,  WCWET
  real, parameter       :: KNFIX = 0, RRUNBULK = 0.05

! Soil - WINTER PARAMETERS
  real*8                 :: FGAS, FO2MX, gamma, KRTOTAER, KSNOW
  real, parameter       :: LAMBDAice      = 1.9354e+005
  real*8                 :: LAMBDAsoil
  real, parameter       :: LatentHeat     = 335000.
  real, parameter       :: poolInfilLimit =      0.2
  real*8                 :: RHOnewSnow, RHOpack
  real, parameter       :: RHOwater       =   1000.
  real*8                 :: SWret, SWrf, TmeltFreeze, TrainSnow
  real*8                 :: WpoolMax

! Soil initial values (parameters)
real*8:: CLITT0, CSOM0, CNLITT0, CNSOMF0, CNSOMS0, FCSOMF0, NMIN0
real*8:: FLITTSOMF, FSOMFSOMS, RNLEACH, KNEMIT
real*8:: TCLITT, TCSOMF, TCSOMS, TMAXF, TSIGMAF, RFN2O
real*8:: WFPS50N2O

! Soil initial constants
  real, parameter       :: DRYSTORI = 0.
  real, parameter       :: FdepthI  = 0.
  real, parameter       :: SDEPTHI  = 0.
  real, parameter       :: TANAERI  = 0.
  real, parameter       :: WAPLI    = 0.
  real, parameter       :: WAPSI    = 0.
  real, parameter       :: WASI     = 0.
  real, parameter       :: WETSTORI = 0.

! Management: harvest dates and irrigation
  integer, dimension(3) :: doyHA
  real, parameter       :: IRRIGF = 0.

! Mathematical constants
  real, parameter       :: pi   = 3.141592653589793
!  real, parameter       :: Freq = 2.*pi / 365.
  real, parameter       :: Kmin = 4.
  real, parameter       :: Ampl = 0.625
  real, parameter       :: Bias = Kmin + Ampl

! SA parameters
  real*8                 :: NFERTMULT

end module parameters_site

