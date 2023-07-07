!#NOTPUBLIC  email:rubinmar@umich.edu  expires:12/31/2099
module ModVarIndexes

  use ModExtraVariables, Redefine1 => Pe_

  implicit none

  save

  ! This equation module contains the standard MHD equations 
  ! with electron pressure
  character (len=*), parameter :: NameEquation= &
       'Multi-fluid Pe MHD for Europa'

  integer, parameter :: nVar = 19

  integer, parameter :: nFluid    = 3
  integer, parameter :: IonFirst_ = 2        ! First individual ion fluid
  integer, parameter :: IonLast_  = 3        ! Last individual ion fluid
  logical, parameter :: IsMhd     = .true.   ! First total ion fluid obeys MHD
  real               :: MassFluid_I(2:nFluid) = (/ 16.0, 32.0 /)

  character (len=3), parameter :: NameFluid_I(nFluid) = &
       (/ 'All', 'Op ', 'O2p' /)


  ! Named indexes for State_VGB and other variables
  ! These indexes should go subsequently, from 1 to nVar+nFluid.
  ! The energy is handled as an extra variable, so that we can use
  ! both conservative and non-conservative scheme and switch between them.
  integer, parameter :: &
       Rho_       =  1, &
       RhoUx_     =  2, Ux_ = 2, &
       RhoUy_     =  3, Uy_ = 3, &
       RhoUz_     =  4, Uz_ = 4, &
       Bx_        =  5, &
       By_        =  6, &
       Bz_        =  7, &
       Pe_        =  8, &
       p_         =  9, &
       OpRho_     = 10, &
       OpRhoUx_   = 11, &
       OpRhoUy_   = 12, &
       OpRhoUz_   = 13, &
       OpP_       = 14, &
       O2pRho_    = 15, &
       O2pRhoUx_  = 16, &
       O2pRhoUy_  = 17, &
       O2pRhoUz_  = 18, &
       O2pP_      = 19, &
       Energy_    = nVar+1, &
       OpEnergy_  = nVar+2, &
       O2pEnergy_ = nVar+3
 
  ! This allows to calculate RhoUx_ as RhoU_+x_ and so on.
  integer, parameter :: U_ = Ux_ - 1, RhoU_ = RhoUx_-1, B_ = Bx_-1

  ! These arrays are useful for multifluid
  integer, parameter :: &
       iRho_I(nFluid)  =(/Rho_,   OpRho_,   O2pRho_ /) ,&
       iRhoUx_I(nFluid)=(/RhoUx_, OpRhoUx_, O2pRhoUx_ /),&
       iRhoUy_I(nFluid)=(/RhoUy_, OpRhoUy_, O2pRhoUy_ /),&
       iRhoUz_I(nFluid)=(/RhoUz_, OpRhoUz_, O2pRhoUz_ /),&
       iP_I(nFluid)    =(/p_,     OpP_,     O2pP_ /)

  ! The default values for the state variables:
  ! Variables which are physically positive should be set to 1,
  ! variables that can be positive or negative should be set to 0:
  real, parameter :: DefaultState_V(nVar+nFluid) = (/ & 
       1.0, & ! Rho_
       0.0, & ! RhoUx_
       0.0, & ! RhoUy_
       0.0, & ! RhoUz_
       0.0, & ! Bx_
       0.0, & ! By_
       0.0, & ! Bz_
       1.0, & ! Pe_
       1.0, & ! p_
       1.0, & ! OpRho_
       0.0, & ! OpRhoUx_
       0.0, & ! OpRhoUy_
       0.0, & ! OpRhoUz_
       1.0, & ! OpP_
       1.0, & ! O2pRho_
       0.0, & ! O2pRhoUx_
       0.0, & ! O2pRhoUy_
       0.0, & ! O2pRhoUz_
       1.0, & ! O2pP_
       1.0, & ! Energy_
       1.0, & ! OpEnergy_
       1.0 /) ! O2pEnergy_
  
  ! The names of the variables used in i/o
  character(len=7) :: NameVar_V(nVar+nFluid) = (/ &
       'Rho    ', & ! Rho_
       'Mx     ', & ! RhoUx_
       'My     ', & ! RhoUy_
       'Mz     ', & ! RhoUz_
       'Bx     ', & ! Bx_
       'By     ', & ! By_
       'Bz     ', & ! Bz_
       'Pe     ', & ! Pe_
       'p      ', & ! p_
       'OpRho  ', & ! OpRho_
       'OpMx   ', & ! OpRhoUx_
       'OpMy   ', & ! OpRhoUy_
       'OpMz   ', & ! OpRhoUz_
       'OpP    ', & ! OpP_
       'O2pRho ', & ! O2pRho_
       'O2pMx  ', & ! O2pRhoUx_
       'O2pMy  ', & ! O2pRhoUy_
       'O2pMz  ', & ! O2pRhoUz_
       'O2pP   ', & ! O2pP_
       'E      ', & ! Energy_
       'OpE    ', & ! OpEnergy_
       'O2pE   ' /) ! O2pEnergy_

  ! The space separated list of nVar conservative variables for plotting
  character(len=*), parameter :: NameConservativeVar = &
       'Rho Mx My Mz Bx By Bz E '// &
       'OpRho OpMx OpMy OpMz OpE '// &
       'O2pRho O2pMx O2pMy O2pMz O2pE '

  ! The space separated list of nVar primitive variables for plotting
  character(len=*), parameter :: NamePrimitiveVar = &
       'Rho Ux Uy Uz Bx By Bz Pe P ' // &
       'OpRho OpUx OpUy OpUz OpP '// &
       'O2pRho O2pUx O2pUy O2pUz O2pP '

  ! The space separated list of nVar primitive variables for TECplot output
  character(len=*), parameter :: NamePrimitiveVarTec = &
       '"`r", "U_x", "U_y", "U_z", "B_x", "B_y", "B_z", "p_e", "p",' // &
       '"`r^O^+", "U_x^O^+", "U_y^O^+", "U_z^O^+", "P^O^+"'// &
       '"`r^O2^+", "U_x^O2^+", "U_y^O2^+", "U_z^O2^+", "P^O2^+"'

  ! Names of the user units for IDL and TECPlot output
  character(len=20) :: &
       NameUnitUserIdl_V(nVar+nFluid) = '', NameUnitUserTec_V(nVar+nFluid) = ''

  ! The user defined units for the variables
  real :: UnitUser_V(nVar+nFluid) = 1.0

  ! There are no extra scalars (Pe has its own flux)
  integer, parameter :: ScalarFirst_ = 2, ScalarLast_ = 1

  ! There are no multi-species
  logical, parameter :: UseMultiSpecies = .false.

  ! Declare the following variables to satisfy the compiler
  integer, parameter :: SpeciesFirst_ = 1, SpeciesLast_ = 1
  real               :: MassSpecies_V(SpeciesFirst_:SpeciesLast_)

contains

  subroutine init_mod_equation

    ! Initialize user units and names for the MHD variables
    call init_mhd_variables

  end subroutine init_mod_equation

end module ModVarIndexes

