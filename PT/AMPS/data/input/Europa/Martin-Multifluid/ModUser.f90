!#NOTPUBLIC  email:rubinmar@umich.edu  expires:12/31/2099
!This code is a copyright protected software (c) 2002- University of Michigan
!========================================================================

module ModUser

  use ModMultiFluid
  use ModSize,       ONLY: nI, nJ, nK, nBlk
  use ModVarIndexes, ONLY: nVar
  use ModAdvance,    ONLY: Pe_, UseElectronPressure
  use ModUserEmpty,                               &
       IMPLEMENTED1  => user_read_inputs,         &
       IMPLEMENTED2  => user_calc_sources,        &
       IMPLEMENTED3  => user_update_states,       &
       IMPLEMENTED4  => user_set_face_boundary,   &
       IMPLEMENTED5  => user_set_resistivity,     &
       IMPLEMENTED6  => user_material_properties, &
       IMPLEMENTED7  => user_init_point_implicit, &
       IMPLEMENTED8  => user_init_session,        &
       IMPLEMENTED9  => user_set_plot_var,        &
       IMPLEMENTED10 => user_set_ICs,             &
       IMPLEMENTED11 => user_get_log_var,         & 
       IMPLEMENTED12 => user_set_boundary_cells,  &
       IMPLEMENTED13 => user_set_cell_boundary

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 0.1
  character (len=*), parameter :: NameUserModule = &
       'Rubin, 2-fluid Europa MHD module, Jul 2013'

  integer, parameter, public :: nNeutral = 1 !! number of neutral species
  !! Neutral species names
  character (len=6), parameter, public :: NameNeutral_I(nNeutral) = &
       (/ 'O2 ' /)
  integer, parameter, public :: O2_  =  1
  !! Ion species names
  integer, parameter, public :: Op_  =  1
  integer, parameter, public :: O2p_ =  2

  real, dimension(nNeutral,MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1) :: NnNeutral_IG, &
       UnxNeutral_IG, UnyNeutral_IG, UnzNeutral_IG, TnNeutral_IG
  real, dimension(1:nNeutral) :: NeutralMass_I(nNeutral)
  real, dimension(1:nIonFluid) :: BodyNDimUp_I, BodyTDimUp_I, &
       BodyNDimDo_I, BodyTDimDo_I
  real, dimension(3) :: Dir_I, SunPos_D

  integer :: iNeutralBlockLast = -1

  real :: nH1, nH2, H1, H2, rcos, mi, Tmin, vO2toOp, vO2toO2p

  logical :: UseResistivePlanet = .false.
  logical, public :: UseFixedBoundary = .false., UseShadow = .false.
  real :: PlanetDensity=-1., PlanetPressure=-1., PlanetRadius=-1., &
       PlanetDensitySi=-1., PlanetPressureSi=-1., PlanetRadiusSi=-1.
  integer :: iLayer, nLayer =0 ! Number of points in planet resistivity profile
  real, allocatable :: PlanetRadiusSi_I(:),PlanetRadius_I(:),&
       ResistivitySi_I(:),Resistivity_I(:)
  real, allocatable :: ResistivityRate(:)

  !! Plotting array to be used for testing
  ! real, dimension(8,MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1,nBLK) :: TestArray = 0.
  real, dimension(8,MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1,nBLK) :: MassLoading = 0.

  !! To plot the cells where the heat flux is capped by the limiter
  ! real, public, dimension(MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1,nBLK) :: FluxLimited_GB = 0.
  !! Add the following lines to get_impl_heat_cond_state 
  ! use ModUser, ONLY: FluxLimited_GB
  !! and after "The threshold heat flux limiter model"
  ! FluxLimited_GB(i,j,k,iBlock) = 0.
  ! if(HeatCoef*GradTe > FreeStreamFlux) &
  !    FluxLimited_GB(i,j,k,iBlock) = 1.
  !! before limiter is applied


contains

  !========================================================================

  subroutine user_read_inputs

    use ModMain,      ONLY: lVerbose
    use ModProcMH,    ONLY: iProc
    use ModReadParam, ONLY: read_var, read_line, read_command
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut

    character (len=100) :: NameCommand
    real :: SunDist

    !-----------------------------------------------------------------------

    if(iProc==0.and.lVerbose > 0)then
       call write_prefix; write(iUnitOut,*)'User read_input Europa starts'
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)

       case("#SUNPOS")
          call read_var('UseShadow', UseShadow)  !! If yes the photoionization is suppressed in Europa's shadow
          if(UseShadow) then
             call read_var('SUNx', SunPos_D(1))  !! x-compoentn of the Europa-Sun vector in GSE coords
             call read_var('SUNy', SunPos_D(2))  !! y-compoentn of the Europa-Sun vector in GSE coords
             call read_var('SUNz', SunPos_D(3))  !! z-compoentn of the Europa-Sun vector in GSE coords
             SunDist=sqrt(sum(SunPos_D(:)**2))   !! Distance for normalization
             SunPos_D(:)=SunPos_D(:)/SunDist     !! Normalization
          end if
       case("#EUROPA")
          call read_var('nH1', nH1)              !! Neutral surface density with H1 scale length [1/cm^3]
          call read_var('H1', H1)                !! Neutral scale height of population 1 [km]
          call read_var('nH2', nH2)              !! Neutral surface density with H2 scale length [1/cm^3]
          call read_var('H2', H2)                !! Neutral scale height of population 2 [km]
          call read_var('rcos', rcos)            !! max fraction in cosine relative to uniform distr. [0..infinity]
          call read_var('vOp', vO2toOp)          !! Ionization rate producing Op from O2 [1/s]
          call read_var('vO2p', vO2toO2p)        !! Ionization rate producing O2p from O2 [1/s]
          call read_var('Tmin', Tmin)            !! Minimum ion temperature (Europa's nightside surface temperature)
          H1=H1*1E3                              !! conversion to SI  
          H2=H2*1E3                              !! conversion to SI  
          nH1=nH1*1E6                            !! conversion to SI
          nH2=nH2*1E6                            !! conversion to SI
       case("#USEFIXEDBOUNDARY")
          call read_var('UseFixedBoundary', UseFixedBoundary) !! UseUpstreamBounday
          if(UseFixedBoundary)then
             do iFluid = 1, nIonFluid
                call read_var('BodyNDimUp', BodyNDimUp_I(iFluid))   !! Upstream particle density
                call read_var('BodyTDimUp', BodyTDimUp_I(iFluid))   !! Upstream plasma temperature
             end do
             do iFluid = 1, nIonFluid
                call read_var('BodyNDimDo', BodyNDimDo_I(iFluid))   !! Downstream particle density
                call read_var('BodyTDimDo', BodyTDimDo_I(iFluid))   !! Downstream plasma temperature
             end do
          end if
       case("#RESISTIVEPLANET")
          call read_var('UseResistivePlanet', UseResistivePlanet)
          call read_var('PlanetDensitySi', PlanetDensitySi)
          call read_var('PlanetPressureSi', PlanetPressureSi)
          call read_var('PlanetRadius', PlanetRadius)
          call read_var('nResistivPoints', nLayer)
          if(nLayer < 2) then
             write(*,*) ' We need minimum 2 points for including resistivity profile'
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
          if(nLayer > 1) then
             allocate(ResistivityRate(nLayer-1),&
                  PlanetRadiusSi_I(nLayer),&
                  PlanetRadius_I(nLayer), &
                  ResistivitySi_I(nLayer),&
                  Resistivity_I(nLayer))
             do iLayer=1,nLayer
                call read_var('Radius',PlanetRadius_I(iLayer))
                call read_var('Resistivity', ResistivitySi_I(iLayer))
             end do
             !Check values
             do iLayer=2,nLayer
                if(PlanetRadius_I(iLayer-1) < &
                     PlanetRadius_I(iLayer)) then
                   write(*,*) 'ERROR: Shoud be decreasing Radius.'
                   call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
                end if
             end do
          end if
       case('#USERINPUTEND')
          if(iProc==0.and.lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input EUROPA ends'
          endif
          EXIT
       case default
          if(iProc==0) then
             call write_myname; write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) '  *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do

  end subroutine user_read_inputs

  !========================================================================
  subroutine user_neutral_atmosphere(iBlock)
    use ModBlockData,  ONLY: get_block_data, set_block_data, put_block_data, &
         use_block_data, MaxBlockData
    use ModPhysics,    ONLY: rPlanetSi, cProtonMass, cPi, rBody, SW_Ux
    use ModMain,       ONLY: nI, nJ, nK, iTest, jTest, kTest, &
         BlkTest, PROCtest, iteration_number, Body1, nIJK, nBlockMax 
    use ModProcMH,     ONLY: iProc 
    use ModGeometry,   ONLY: R_BLK, Xyz_DGB
    use ModMultiFluid, ONLY: MassIon_I
    use ModIO,         ONLY: iUnitOut

    integer,intent(in) :: iBlock
    logical :: DoTest, DoTestMe=.true., init=.true., FirstCall=.true.
    real :: theta, cos_theta, term, Dir_length

    integer :: i, j, k, iNeutral

    ! Dayside neutral distribution centered at vector direction
    Dir_I(1) =  -SW_Ux ! for direction of undisturbed plasma inflow
    Dir_I(2) =  0.0 !! -SW_Uz
    Dir_I(3) =  0.0 !! -SW_Uz

    !----------------------------------------------------------------------

    if(iProc==PROCtest .and. iBlock == BlkTest .and. FirstCall) then
       call set_oktest('user_neutral_atmosphere',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    ! Max number of reals saved in ModBlockData (n,ux,uy,uz,T)
    ! MaxBlockData = nNeutral*5*nIJK*nBlockMax

    if (.not.use_block_data(iBlock)) then
       ! ! calculate and print out total mass loading (integral from rE to infinity)
       ! if(iProc==0.and.init) then
       !    term = 2*nH2*H2**2*rPlanetSi+nH1*H1*rPlanetSi**2+2*nH2*H2**3+& 
       !         2*nH1*H1**2*rPlanetSi+2*nH1*H1**3+nH2*H2*rPlanetSi**2 
       !    write(iUnitOut,*)'Total Mass Loading = ',&
       !         term*cPi*(4+rcos)*(vO2toOp*MassIon_I(Op_)+vO2toO2p*MassIon_I(O2p_))&
       !         *cProtonMass,' [kg/s]'
       !    init=.false.
       ! end if

       Dir_length = sqrt(dot_product(Dir_I,Dir_I))

       do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
          ! Initialize neutral atmosphere
          NnNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1)  = 0.
          UnxNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1) = 0.
          UnyNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1) = 0.
          UnzNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1) = 0.
          TnNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1)  = 0.

          if((Body1).and.(R_BLK(i,j,k,iBlock) < rBody)) CYCLE
          if((R_BLK(i,j,k,iBlock) < PlanetRadius).and.(UseResistivePlanet)) CYCLE
          ! angle of cell position relative to ram direction
          cos_theta=(Dir_I(1)*Xyz_DGB(x_,i,j,k,iBlock)+Dir_I(2)*Xyz_DGB(y_,i,j,k,iBlock)+&
               Dir_I(3)*Xyz_DGB(z_,i,j,k,iBlock))/(R_BLK(i,j,k,iBlock)*Dir_length)

          ! two symmetric distributions w/ different scale heights
          NnNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1) = &
               nH1*exp(-(R_BLK(i,j,k,iBlock)-1.)/(H1/rPlanetSi))+ & ! H1 scale length contribution
               nH2*exp(-(R_BLK(i,j,k,iBlock)-1.)/(H2/rPlanetSi))    ! H2 scale length contribution

          !! ??? test w/o body on
          !if (NnNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1) > nH1+nH2) then
          !   NnNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1) = nH1+nH2
          !end if

          ! relative increase in relation to symmetric part, cosine distribution
          if(cos_theta>=0.) then
             NnNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1) = NnNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1)*(1.+rcos*cos_theta) 
          end if
          UnxNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1) = 0.   !! set bulk flow speed of neutral gas to 0 m/s
          UnyNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1) = 0.
          UnzNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1) = 0.
          TnNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1)  = 600.  !! Kliore et al., Science, 277, 355-358, 1997
       end do;  end do;  end do

       NeutralMass_I(O2_) = 32.*cProtonMass

       call put_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, NnNeutral_IG, &
            DoAllowReplace=.true.)
       call put_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, UnxNeutral_IG, &
            DoAllowReplace=.true.)
       call put_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, UnyNeutral_IG, &
            DoAllowReplace=.true.)
       call put_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, UnzNeutral_IG, &
            DoAllowReplace=.true.)
       call put_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, TnNeutral_IG, &
            DoAllowReplace=.true.)

       call set_block_data(iBlock); !! has to be set in case data is accessed before first iteration is finished

       if(DoTestMe) then
          write(*,*)'user_neutral_atmosphere:'
          theta=acos((Dir_I(1)*Xyz_DGB(x_,iTest,jTest,kTest,BlkTest)+Dir_I(2)*Xyz_DGB(y_,iTest,jTest,kTest,BlkTest)+&
               Dir_I(3)*Xyz_DGB(z_,iTest,jTest,kTest,BlkTest))/(R_BLK(iTest,jTest,kTest,BlkTest)*Dir_length))
          write(*,*)'x      = ',Xyz_DGB(x_,iTest,jTest,kTest,BlkTest)," [rPlanet]"
          write(*,*)'y      = ',Xyz_DGB(y_,iTest,jTest,kTest,BlkTest)," [rPlanet]"
          write(*,*)'z      = ',Xyz_DGB(z_,iTest,jTest,kTest,BlkTest)," [rPlanet]"
          write(*,*)'r      = ',R_BLK(iTest,jTest,kTest,BlkTest)," [rPlanet]"
          write(*,*)'theta  = ',theta," [radians]"
          write(*,*)'costhe = ',cos_theta," []"
          write(*,*)'rcos   = ',rcos," [ ]"
          write(*,*)'nH1    = ',nH1," [m^-3]"
          write(*,*)'H1     = ',H1," [m]"
          write(*,*)'nH2    = ',nH2," [m^-3]"
          write(*,*)'H2     = ',H2," [m]"
          do iNeutral=1,nNeutral
             write(*,*)'Neutral species # ',iNeutral,': ', NameNeutral_I(iNeutral)
             write(*,*)'n_n    = ',NnNeutral_IG(iNeutral,iTest-MinI+1,jTest-MinJ+1,kTest-MinK+1)," [m^-3]"
             write(*,*)'m_n    = ',NeutralMass_I(iNeutral)," [kg]"
             write(*,*)'unx    = ',UnxNeutral_IG(iNeutral,iTest-MinI+1,jTest-MinJ+1,kTest-MinK+1)," [m/s]"
             write(*,*)'uny    = ',UnyNeutral_IG(iNeutral,iTest-MinI+1,jTest-MinJ+1,kTest-MinK+1)," [m/s]"
             write(*,*)'unz    = ',UnzNeutral_IG(iNeutral,iTest-MinI+1,jTest-MinJ+1,kTest-MinK+1)," [m/s]"
             write(*,*)'Tn     = ',TnNeutral_IG(iNeutral,iTest-MinI+1,jTest-MinJ+1,kTest-MinK+1)," [K]"
          end do
          write(*,*)''
          FirstCall = .false.
       end if

    else
       call get_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, NnNeutral_IG)
       call get_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, UnxNeutral_IG)
       call get_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, UnyNeutral_IG)
       call get_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, UnzNeutral_IG)
       call get_block_data(iBlock, nNeutral, MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1, TnNeutral_IG)
    end if

  end subroutine user_neutral_atmosphere

  !========================================================================

  subroutine calc_electron_collision_rates(Te,nElec,i,j,k,iBlock,fen_I,fei_I)

    ! calculate all collision rates for electrons (fen, fei)
    ! (used for sources & resistivity)

    use ModAdvance,    ONLY: State_VGB, Rho_
    use ModPhysics,    ONLY: No2SI_V, UnitN_
    use ModMultiFluid, ONLY: MassIon_I, ChargeIon_I

    integer,intent(in) :: i,j,k,iBlock   
    real,intent(in)    :: Te
    real,intent(in)    :: nElec
    real,intent(out)   :: fen_I(nNeutral)
    real,intent(out)   :: fei_I(nIonFluid)

    real :: sqrtTe
    !----------------------------------------------------------------------

    !! electron-neutral and electron-ion collision rates from Schunk and Nagy, Ionospheres,
    !! Cambridge University Press, 2000
    sqrtTe = sqrt(Te)

    !! initialize all collision rates with zero
    fei_I = 0. ; fen_I = 0.

    !! e - O2, neutral density in [cm^-3]
    fen_I(O2_) = 1.82E-10*NnNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1)/1E6*(1.0+3.6E-2*sqrtTe)*sqrtTe

    !! e - Sp, ion density in [cm^-3]
    fei_I(Op_) = 54.5*ChargeIon_I(Op_)**2*State_VGB(OpRho_,i,j,k,iBlock)/MassIon_I(Op_)* &
         No2SI_V(UnitN_)/1E6/(Te*sqrtTe)                                                   
    !! e - O2p, ion density in [cm^-3]
    fei_I(O2p_) = 54.5*ChargeIon_I(O2p_)**2*State_VGB(O2pRho_,i,j,k,iBlock)/MassIon_I(O2p_)* &
         No2SI_V(UnitN_)/1E6/(Te*sqrtTe)                                                   

  end subroutine calc_electron_collision_rates

  !========================================================================
  subroutine user_calc_rates(Ti_I,Te,i,j,k,iBlock,nElec,nIon_I,fin_II,fii_II,fie_I,alpha_I,kin_IIII,v_II,&
       ve_II,uElec_D,uIon_DI,Qexc_II,Qion_II)

    ! calculate all rates not involving electron collisions

    use ModPhysics,  ONLY: SI2No_V, UnitN_, rPlanetSI, rBody
    use ModConst,    ONLY: cElectronCharge, cBoltzmann, cElectronMass, cProtonMass
    use ModMain,     ONLY: Body1, iTest, jTest, kTest, BlkTest
    use ModNumConst, ONLY: cPi, cHalfPi
    use ModGeometry, ONLY: R_BLK, Xyz_DGB

    integer,intent(in) :: i,j,k,iBlock
    real,intent(in)    :: Ti_I(nIonFluid)
    real,intent(in)    :: Te
    real,intent(in)    :: nElec
    real,intent(in)    :: nIon_I(nIonFluid)
    real,intent(in)    :: uElec_D(3)
    real,intent(in)    :: uIon_DI(3,nIonFluid)
    real,intent(out)   :: fin_II(nIonFluid,nNeutral)
    real,intent(out)   :: fii_II(nIonFluid,nIonFluid)
    real,intent(out)   :: fie_I(nIonFluid)
    real,intent(out)   :: alpha_I(nIonFluid)
    real,intent(out)   :: kin_IIII(nIonFluid,nNeutral,nNeutral,nIonFluid)
    real,intent(out)   :: v_II(nNeutral,nIonFluid)
    real,intent(out)   :: ve_II(nNeutral,nIonFluid)
    real,intent(out)   :: Qexc_II(nNeutral,nIonFluid)
    real,intent(out)   :: Qion_II(nNeutral,nIonFluid)

    real :: Tr, sqrtTe, phi
    !----------------------------------------------------------------------

    !! ion-neutral collision rate from Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    !! densities are in 1/ccm

    !! initialize all collision rates with zero
    fin_II = 0. ; fii_II = 0. ; kin_IIII = 0. ; v_II = 0. ; ve_II = 0. ; alpha_I = 0. ; Qexc_II = 0. ; Qion_II = 0.

    !! Electron excess energies from photoionization (increases electron pressure)
    Qexc_II(O2_,Op_)  = 3.84e-18 ! 24.0 eV Huebner
    Qexc_II(O2_,O2p_) = 2.56e-18 ! 16.0 eV Huebner

    !! Ionization potential for electron impact ionization (needs to be delivered by the ionizing electron)
    Qion_II(O2_,Op_)  = 3.00e-18 ! 18.8 eV (estimate) 5.15 eV bond strength O=O and 13.62 eV O -> O+
    Qion_II(O2_,O2p_) = 1.92e-18 ! 12.0 eV (Samson and Gardener, On the ionization potential of molecular oxygen, 1974)

    !! ********** Ion-neutral collision/charge exchange rates ********** 
    !! Example(s)
    ! resonant H+ & O -> O+ & H  subtracts H+ and adds O+
    ! kin_IIII(Hp_,O_,Op_,H_) = 6.61E-11/1E6*sqrt(Ti_I(Hp_))*(1.0-0.047*log10(Ti_I(Hp)))**2    !! rate in [m^3/s]
    ! resonant O+ & H -> H+ & O  subtracts O+ and adds H+
    ! kin_IIII(Op_,H_,Hp_,O_) = 4.63E-12/1E6*sqrt(TnNeutral(H_,i-MinI+1,j-MinJ+1,k-MinK+1)+TOp_/16.)  !! rate in [m^3/s]


    !! Op & O2 -> Op & O2    ! non-resonant
    fin_II(Op_,O2_) = 6.64E-10*NnNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1)/1E6          !! rate in [s^-1]
    !! Op & O2 -> O & O2p    ! resonant
    !kin_IIII(Op_,O2_,O_,O2p_) = 6.64E-10/1E6                                     !! rate in [m^3/s]

    !! O2p & O2 -> O2 & O2p  ! resonant
    Tr=(Ti_I(O2p_)+TnNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1))/2.

    kin_IIII(O2p_,O2_,O2_,O2p_) = 2.59E-11/1E6*sqrt(Tr)*(1.0-0.073*log10(Tr))**2 !! rate in [m^3/s]
    !! O2p & O2 -> O2p & O2  ! non-resonant
    !fin_II(O2p_,O2_) = 2.59E-11*NnNeutral_IG(O2_,i-MinI+1,j-MinJ+1,k-MinK+1)/1E6* &     !! rate in [s^-1]
    !     sqrt(Tr)*(1.0-0.073*log10(Tr))**2                         


    !! ********** Ion-ion collision rates ********** 
    ! Op - Op is left zero because the they do not result in a change in the source terms anyways
    ! Op - O2p
    fii_II(Op_,O2p_) = 0.26*nIon_I(O2p_)/1E6/(Ti_I(O2p_)*sqrt(Ti_I(O2p_)))
    ! O2p - O2p is left zero because the they do not result in a change in the source terms anyways
    ! O2p - Op
    fii_II(O2p_,Op_) = 0.13*nIon_I(Op_)/1E6/(Ti_I(Op_)*sqrt(Ti_I(Op_)))


    !! Ion - electron collision rates, reduced mass=~me and reduced temperature=~Te
    sqrtTe = sqrt(Te)
    !! Op - e, Schunk and Nagy, Ionospheres, Cambridge University Press, 2000
    fie_I(Op_) = 1.27*sqrt(cElectronMass/cProtonMass)/MassIon_I(Op_)*ChargeIon_I(Op_)**2*nElec/ &
         1E6/(Te*sqrtTe)      !! rate in [1/s]
    !! O2p - e, Schunk and Nagy, Ionospheres, Cambridge University Press, 2000
    fie_I(O2p_) = 1.27*sqrt(cElectronMass/cProtonMass)/MassIon_I(O2p_)*ChargeIon_I(O2p_)**2*nElec/ &
         1E6/(Te*sqrtTe)      !! rate in [1/s]

    !! ********** Ion-electron recombination rates ********** 
    alpha_I(Op_)   = 1E-6*3.7E-12*(250/Te)**0.7                                !! rate in [m^3/s]
    alpha_I(O2p_)  = 1E-6*2.4E-7*(300/Te)**0.7                                 !! rate in [m^3/s]

    !! Photoionization rates [Huebner et al. (1996) at 5.2 AU]
    v_II(O2_,O2p_) = 4.6e-7/(5.2**2)                                           !! rate in [s^-1]
    v_II(O2_,Op_)  = 1.1e-7/(5.2**2)                                           !! rate in [s^-1]

    !! Set photo-rates to zero in Europa's wake
    if(UseShadow) then
       !! angle phi between anti-Sun direction and cell
       phi = acos(sum(-SunPos_D(:)*Xyz_DGB(:,i,j,k,iBlock))/R_BLK(i,j,k,iBlock))
       if(phi < cHalfPi) then
          if (R_BLK(i,j,k,iBlock)*sin(phi) < 1.0) then
             v_II = 0.0
          end if
       end if
    end if

    !! Electron impact ionization rate from PARAM.in
    !ve_II(O2_,Op_)  = vO2toOp                                                  !! rate in [s^-1]
    !ve_II(O2_,O2p_) = vO2toO2p                                                 !! rate in [s^-1]

    !! Electron impact ionization rate using Arrhenius' equation [Schreier et al. (1993)]
    !! ve_II(O2_,O2p_) = 1E-6*9.43e-10*sqrtTe/17.32*exp(-59417./Te)*nElec         !! rate in [s^-1]

    !! Fit to Schilling PhD thesis multiplied by sqrtPi for correction
    if (Te > 50000) then
       ve_II(O2_,O2p_) = 1.77855e-13*((Te-27737.6)/300.)**0.082823* &
            exp(-379245.5/(Te+27737.6))*nElec                                  !! rate in [s^-1]
       !write(*,*)'Te = ', Te, 've = ', ve_II(O2_,O2p_)/nElec
    end if
    ve_II(O2_,Op_)  = 0.1*ve_II(O2_,O2p_) ! estimate                           !! rate in [s^-1]

    !! Use values from input file
    ! ve_II(O2_,O2p_) = vO2toO2p
    ! ve_II(O2_,Op_)  = vO2toOp

    !! Complete separation of magnetospheric and pick-up ions
    !! ve_II(O2_,Op_)   = 1e-9*ve_II(O2_,O2p_) ! no massloading of O+               !! rate in [s^-1]
    !! v_II(O2_,Op_)    = 1e-9*v_II(O2_,Op_)   ! no massloading of O+               !! rate in [s^-1]
    !! ve_II(O2_,O2p_)  = 1.1*ve_II(O2_,O2p_)  ! incl. 10% O2->Op magnetospheric Op !! rate in [s^-1]

    !! Total ionization rate (photo & electron impact)
    v_II(O2_,Op_)  = v_II(O2_,Op_) + ve_II(O2_,Op_)
    v_II(O2_,O2p_) = v_II(O2_,O2p_) + ve_II(O2_,O2p_)

  end subroutine user_calc_rates

  !========================================================================

  subroutine user_calc_sources(iBlock)

    use ModMain,       ONLY: nI, nJ, nK, iTest, jTest, kTest, &
         BlkTest, PROCtest, iteration_number, Dt_BLK
    use ModAdvance,    ONLY: State_VGB, Source_VC, Rho_, RhoUx_, RhoUy_, RhoUz_, &
         Bx_,By_,Bz_, P_, Energy_
    use ModConst,      ONLY: cBoltzmann, cElectronMass, cElectronCharge, cProtonMass
    use ModGeometry,   ONLY: Rmin_BLK, r_BLK, Xyz_DGB
    use ModCurrent,    ONLY: get_current
    use ModProcMH,     ONLY: iProc
    use ModPhysics,    ONLY: SW_Ux, SW_Uy, SW_Uz, UnitN_, UnitRho_, UnitU_, UnitP_, UnitT_, UnitB_, &
         ElectronPressureRatio, ElectronCharge, Si2No_V, No2Si_V, UnitEnergyDens_, UnitJ_, UnitRhoU_, UnitX_
    use ModPointImplicit, ONLY: UsePointImplicit_B, UsePointImplicit, IsPointImplSource

    integer, intent(in) :: iBlock

    real, dimension(1:nI,1:nJ,1:nK) :: nElec_C, Te_C, SBx_C, SBy_C, SBz_C, SPe_C
    real, dimension(4,1:nIonFluid,1:nI,1:nJ,1:nK) :: SRhoTerm_IIC
    real, dimension(5,1:nIonFluid,1:nI,1:nJ,1:nK) :: SRhoUxTerm_IIC, SRhoUyTerm_IIC, SRhoUzTerm_IIC
    real, dimension(8,1:nIonFluid,1:nI,1:nJ,1:nK) :: SPTerm_IIC
    real, dimension(7,1:nI,1:nJ,1:nK) :: SPeTerm_IC

    real, dimension(1:3,1:nI,1:nJ,1:nK) :: Current_DC, uIonMean_DC, uElec_DC
    real, dimension(1:3,1:nIonFluid,1:nI,1:nJ,1:nK) :: uIon_DIC
    real, dimension(1:nIonFluid,1:nNeutral,1:nI,1:nJ,1:nK) :: fin_IIC, uIonNeu2_IIC
    real, dimension(1:nNeutral,1:nI,1:nJ,1:nK) :: fen_IC, uNeuElec2_IC
    real, dimension(1:nIonFluid,1:nIonFluid,1:nI,1:nJ,1:nK) :: fii_IIC, uIonIon2_IIC
    real, dimension(1:nIonFluid,1:nI,1:nJ,1:nK) :: Ti_IC, uIonElec2_IC, fei_IC, fie_IC, &
         nIon_IC, SRho_IC, SRhoUx_IC, SRhoUy_IC, SRhoUz_IC, SP_IC
    real, dimension(1:nNeutral,1:nIonFluid) :: Qexc_II, Qion_II
    real, dimension(1:nNeutral,1:nIonFluid,1:nI,1:nJ,1:nK) :: v_IIC, ve_IIC
    real, dimension(1:nIonFluid,1:nI,1:nJ,1:nK) :: alpha_IC
    real, dimension(1:nIonFluid) :: fiiTot_I, finTot_I, vAdd_I, kinAdd_I, kinSub_I
    real, dimension(1:nIonFluid,1:nNeutral,1:nNeutral,1:nIonFluid,1:nI,1:nJ,1:nK) :: kin_IIIIC

    logical :: DoTest, DoTestMe=.true.
    real :: theta, fenTot, feiTot,logTe
    integer :: i,j,k,iNeutral,jNeutral,iIonFluid,jIonFluid,iTerm,iDim

    !----------------------------------------------------------------------

    ! Do not evaluate any source terms explicitly when running pointimplicit
    if(UsePointImplicit .and. .not. IsPointImplSource) RETURN

    ! Evaluate source terms explicitly even when running pointimplicit
    !if(UsePointImplicit .and. IsPointImplSource) RETURN

    !! Limit region for evaluation for source term evaluation
    !! if(RMin_BLK(iBlock) > 2.) RETURN

    if(iBlock == BlkTest.and.iProc==ProcTest) then
       call set_oktest('user_calc_sources',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    if (iBlock.ne.iNeutralBlockLast) then
       call user_neutral_atmosphere(iBlock)
    end if

    !! Set the source arrays for this block to zero
    SRho_IC        = 0.
    SRhoTerm_IIC   = 0.
    SRhoUx_IC      = 0.
    SRhoUxTerm_IIC = 0.
    SRhoUy_IC      = 0.
    SRhoUyTerm_IIC = 0.
    SRhoUz_IC      = 0.
    SRhoUzTerm_IIC = 0.
    SBx_C          = 0.
    SBy_C          = 0.
    SBz_C          = 0.
    SP_IC          = 0.
    SPTerm_IIC     = 0.
    SPe_C          = 0.
    SPeTerm_IC     = 0.


    ! nElec_C is the electron/ion density in SI units ( n_e=sum(n_i*Zi) )
    do k=1,nK; do j=1,nJ; do i=1,nI
       nIon_IC(1:nIonFluid,i,j,k) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
       nElec_C(i,j,k) = sum(nIon_IC(1:nIonFluid,i,j,k)*ChargeIon_I(1:nIonFluid))
    end do; end do; end do

    !! ion velocity components in SI

    uIon_DIC(1,1:nIonFluid,1:nI,1:nJ,1:nK)=State_VGB(iRhoUxIon_I,1:nI,1:nJ,1:nK,iBlock) / &
         State_VGB(iRhoIon_I,1:nI,1:nJ,1:nK,iBlock)*No2SI_V(UnitU_)
    uIon_DIC(2,1:nIonFluid,1:nI,1:nJ,1:nK)=State_VGB(iRhoUyIon_I,1:nI,1:nJ,1:nK,iBlock) / &
         State_VGB(iRhoIon_I,1:nI,1:nJ,1:nK,iBlock)*No2SI_V(UnitU_)
    uIon_DIC(3,1:nIonFluid,1:nI,1:nJ,1:nK)=State_VGB(iRhoUzIon_I,1:nI,1:nJ,1:nK,iBlock) / &
         State_VGB(iRhoIon_I,1:nI,1:nJ,1:nK,iBlock)*No2SI_V(UnitU_)
    uIonMean_DC(1:3,1:nI,1:nJ,1:nK) = 0.
    do iIonFluid=1,nIonFluid
       uIonMean_DC(1,1:nI,1:nJ,1:nK) = uIonMean_DC(1,1:nI,1:nJ,1:nK)+nIon_IC(iIonFluid,1:nI,1:nJ,1:nK)* &
            uIon_DIC(1,iIonFluid,1:nI,1:nJ,1:nK)/nElec_C(1:nI,1:nJ,1:nK)*ChargeIon_I(iIonFluid)
       uIonMean_DC(2,1:nI,1:nJ,1:nK) = uIonMean_DC(2,1:nI,1:nJ,1:nK)+nIon_IC(iIonFluid,1:nI,1:nJ,1:nK)* &
            uIon_DIC(2,iIonFluid,1:nI,1:nJ,1:nK)/nElec_C(1:nI,1:nJ,1:nK)*ChargeIon_I(iIonFluid)
       uIonMean_DC(3,1:nI,1:nJ,1:nK) = uIonMean_DC(3,1:nI,1:nJ,1:nK)+nIon_IC(iIonFluid,1:nI,1:nJ,1:nK)* &
            uIon_DIC(3,iIonFluid,1:nI,1:nJ,1:nK)/nElec_C(1:nI,1:nJ,1:nK)*ChargeIon_I(iIonFluid) 
    end do

    !! (u_i-u_n)^2 in SI
    do iIonFluid=1,nIonFluid
       do iNeutral=1,nNeutral
          uIonNeu2_IIC(iIonFluid,iNeutral,1:nI,1:nJ,1:nK) = &
               (uIon_DIC(1,iIonFluid,1:nI,1:nJ,1:nK)- &
               UnxNeutral_IG(iNeutral,2-MinI:nI+1-MinI,2-MinJ:nJ+1-MinJ,2-MinK:nK+1-MinK))**2+&
               (uIon_DIC(2,iIonFluid,1:nI,1:nJ,1:nK)- &
               UnyNeutral_IG(iNeutral,2-MinI:nI+1-MinI,2-MinJ:nJ+1-MinJ,2-MinK:nK+1-MinK))**2+&
               (uIon_DIC(3,iIonFluid,1:nI,1:nJ,1:nK)- &
               UnzNeutral_IG(iNeutral,2-MinI:nI+1-MinI,2-MinJ:nJ+1-MinJ,2-MinK:nK+1-MinK))**2
       end do
    end do

    !! (u_i1-u_i2)^2 in SI
    do iIonFluid=1,nIonFluid
       do jIonFluid=1,nIonFluid
          uIonIon2_IIC(iIonFluid,jIonFluid,:,:,:) = &
               (uIon_DIC(1,iIonFluid,1:nI,1:nJ,1:nK)-uIon_DIC(1,jIonFluid,1:nI,1:nJ,1:nK))**2+&
               (uIon_DIC(2,iIonFluid,1:nI,1:nJ,1:nK)-uIon_DIC(2,jIonFluid,1:nI,1:nJ,1:nK))**2+&
               (uIon_DIC(3,iIonFluid,1:nI,1:nJ,1:nK)-uIon_DIC(3,jIonFluid,1:nI,1:nJ,1:nK))**2
       end do
    end do

    if (UseElectronPressure) then
       ! Electron temperature calculated from electron pressure
       ! Ion temperature is calculated from ion pressure
       do k=1,nK; do j=1,nJ; do i=1,nI
          Ti_IC(1:nIonFluid,i,j,k) = State_VGB(iPIon_I,i,j,k,iBlock)*NO2SI_V(UnitP_)/&
               (cBoltzmann*State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*NO2SI_V(UnitN_))
          Te_C(i,j,k) = State_VGB(Pe_,i,j,k,iBlock)*NO2SI_V(UnitP_)/(cBoltzmann* &
               nElec_C(i,j,k))
       end do; end do; end do
    else
       ! Electron temperature calculated from pressure assuming Te_C=Ti_IC*ElectronTemperatureRatio:
       ! p=nkT with n_e=n_i*Z_i (quasi-neutrality), n=n_e+n_i and p=p_e+p_i=p_i*(1+ElectronPressureRatio)
       do k=1,nK; do j=1,nJ; do i=1,nI
          Ti_IC(1:nIonFluid,i,j,k) = State_VGB(iPIon_I,i,j,k,iBlock)*NO2SI_V(UnitP_)/ &
               (cBoltzmann*State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*NO2SI_V(UnitN_))
          Te_C(i,j,k) = State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio/(1.+ElectronPressureRatio)*&
               NO2SI_V(UnitP_)/(cBoltzmann*nElec_C(i,j,k))
       end do; end do; end do
    end if

    do k=1,nK; do j=1,nJ; do i=1,nI
       ! No need to evaluate source terms for cells between conducting core and surface
       if((R_BLK(i,j,k,iBlock) < PlanetRadius).and.(UseResistivePlanet)) CYCLE

       call get_current(i,j,k,iBlock,Current_DC(:,i,j,k))

       ! calculate uElec_DC from Hall velocity -J/(e*n) [m/s]
       uElec_DC(1:3,i,j,k) = uIonMean_DC(1:3,i,j,k)-Current_DC(1:3,i,j,k)/(nElec_C(i,j,k)*Si2No_V(UnitN_)*&
            ElectronCharge)*No2SI_V(UnitU_)

       call calc_electron_collision_rates(Te_C(i,j,k),nElec_C(i,j,k),i,j,k,iBlock,fen_IC(1:nNeutral,i,j,k), &
            fei_IC(1:nIonFluid,i,j,k))
       call user_calc_rates(Ti_IC(1:nIonFluid,i,j,k),Te_C(i,j,k),i,j,k,iBlock,nElec_C(i,j,k),nIon_IC(1:nIonFluid,i,j,k),&
            fin_IIC(1:nIonFluid,1:nNeutral,i,j,k),fii_IIC(1:nIonFluid,1:nIonFluid,i,j,k),fie_IC(1:nIonFluid,i,j,k),&
            alpha_IC(1:nIonFluid,i,j,k),kin_IIIIC(1:nIonFluid,1:nNeutral,1:nNeutral,1:nIonFluid,i,j,k),&
            v_IIC(1:nNeutral,1:nIonFluid,i,j,k),ve_IIC(1:nNeutral,1:nIonFluid,i,j,k),uElec_DC(1:3,i,j,k),&
            uIon_DIC(1:3,1:nIonFluid,i,j,k),Qexc_II(1:nNeutral,1:nIonFluid),Qion_II(1:nNeutral,1:nIonFluid))

       !! Zeroth moment
       !! Sources separated into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"       
       kinAdd_I = 0. ; kinSub_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! addition to individual fluid from charge exchange [1/(m^3*s)]
                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
                   !! subtraction to individual fluid from charge exchange [1/(m^3*s)]
                   kinSub_I(iIonFluid) = kinSub_I(iIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
                end do
             end do
          end do
       end do

       vAdd_I = 0.
       do iNeutral=1,nNeutral
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
       end do

       !! Sources divideded into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
       SRhoTerm_IIC(1,1:nIonFluid,i,j,k) = vAdd_I(1:nIonFluid)*Si2No_V(UnitN_)/Si2No_V(UnitT_)*MassIon_I              !! newly ionized neutrals
       SRhoTerm_IIC(2,1:nIonFluid,i,j,k) = kinAdd_I(1:nIonFluid)*Si2No_V(UnitN_)/Si2No_V(UnitT_)*MassIon_I            !! mass added through ion-neutral charge exchange
       SRhoTerm_IIC(3,1:nIonFluid,i,j,k) = -kinSub_I(1:nIonFluid)*Si2No_V(UnitN_)/Si2No_V(UnitT_)*MassIon_I           !! mass removed through ion-neutral charge exchange
       SRhoTerm_IIC(4,1:nIonFluid,i,j,k) = -alpha_IC(1:nIonFluid,i,j,k)*(nElec_C(i,j,k)* &                            !! loss due to recombination
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitN_)/Si2No_V(UnitT_))*MassIon_I

       !! First moment, x component
       !! d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt combined from zeroth and first moment by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
       fiiTot_I = 0. ; finTot_I = 0. ; vAdd_I = 0.
       do iIonFluid=1,nIonFluid                                                                                       !! momentum transfer by ion-ion collisions
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*&                        !! ion-ion collisions
               (uIon_DIC(1,iIonFluid,i,j,k)-uIon_DIC(1,1:nIonFluid,i,j,k))
       end do                                                                                                         !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*&                         !! ion-neutral collisions
               (UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(1,1:nIonFluid,i,j,k))
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
               (UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(1,1:nIonFluid,i,j,k))
       end do

       kinAdd_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! addition to individual fluid from charge exchange [1/(m^3*s)]
                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
                        (UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(1,jIonFluid,i,j,k))
                end do
             end do
          end do
       end do

       SRhoUxTerm_IIC(1,1:nIonFluid,i,j,k) = (vAdd_I(1:nIonFluid)/Si2No_V(UnitT_)+ &                                  !! newly photoionized neutrals
            kinAdd_I(1:nIonFluid)/Si2No_V(UnitT_))*Si2No_V(UnitN_)*MassIon_I* &                                       !! new ions from charge exchange
            Si2No_V(UnitU_)
       ! Add u_s*drho_s/dt for d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt
       do iIonFluid=1,nIonFluid
          SRhoUxTerm_IIC(1,iIonFluid,i,j,k) = SRhoUxTerm_IIC(1,iIonFluid,i,j,k)+sum(SRhoTerm_IIC(1:4,iIonFluid,i,j,k))&
               *uIon_DIC(1,iIonFluid,i,j,k)*Si2No_V(UnitU_)
       end do
       SRhoUxTerm_IIC(2,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)/ElectronCharge* &             !! current dissipation, ion-electron collisions
            MassIon_I*nIon_IC(1:nIonFluid,i,j,k)/nElec_C(i,j,k)*Current_DC(1,i,j,k)
       SRhoUxTerm_IIC(3,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*MassIon_I*&
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitN_)*(uIon_DIC(1,1:nIonFluid,i,j,k)-uIonMean_DC(1,i,j,k))*Si2No_V(UnitU_)
       SRhoUxTerm_IIC(4,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*fiiTot_I(1:nIonFluid)* &    !! ion-ion collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I
       SRhoUxTerm_IIC(5,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*finTot_I(1:nIonFluid)* &    !! ion neutral collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I

       !! First moment, y component
       !! Sources separated into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
       fiiTot_I = 0. ; finTot_I = 0. ; vAdd_I = 0.
       do iIonFluid=1,nIonFluid                                                                                       !! momentum transfer by ion-ion collisions
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*&                        !! ion-ion collisions
               (uIon_DIC(2,iIonFluid,i,j,k)-uIon_DIC(2,1:nIonFluid,i,j,k))
       end do                                                                                                         !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*&                         !! ion-neutral collisions
               (UnyNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(2,1:nIonFluid,i,j,k))
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
               (UnyNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(2,1:nIonFluid,i,j,k))
       end do


       kinAdd_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! addition to individual fluid from charge exchange [1/(m^3*s)]
                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*&
                        NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
                        (UnyNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(2,jIonFluid,i,j,k))
                end do
             end do
          end do
       end do

       SRhoUyTerm_IIC(1,1:nIonFluid,i,j,k) = (vAdd_I(1:nIonFluid)/Si2No_V(UnitT_)+ &                                  !! newly photoionized neutrals
            kinAdd_I(1:nIonFluid)/Si2No_V(UnitT_))*Si2No_V(UnitN_)*MassIon_I* &                                       !! new ions from charge exchange
            Si2No_V(UnitU_)
       ! Add u_s*drho_s/dt for d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt
       do iIonFluid=1,nIonFluid
          SRhoUyTerm_IIC(1,iIonFluid,i,j,k) = SRhoUyTerm_IIC(1,iIonFluid,i,j,k)+sum(SRhoTerm_IIC(1:4,iIonFluid,i,j,k))&
               *uIon_DIC(2,iIonFluid,i,j,k)*Si2No_V(UnitU_)
       end do
       SRhoUyTerm_IIC(2,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)/ElectronCharge* &             !! current dissipation, ion-electron collisions
            MassIon_I*nIon_IC(1:nIonFluid,i,j,k)/nElec_C(i,j,k)*Current_DC(2,i,j,k)
       SRhoUyTerm_IIC(3,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*MassIon_I*&
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitN_)*(uIon_DIC(2,1:nIonFluid,i,j,k)-uIonMean_DC(2,i,j,k))*Si2No_V(UnitU_)
       SRhoUyTerm_IIC(4,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*fiiTot_I(1:nIonFluid)* &    !! ion-ion collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I
       SRhoUyTerm_IIC(5,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*finTot_I(1:nIonFluid)* &    !! ion neutral collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I

       !! First moment, z component
       !! Sources separated into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
       fiiTot_I = 0. ; finTot_I = 0. ; vAdd_I = 0.
       do iIonFluid=1,nIonFluid                                                                                       !! momentum transfer by ion-ion collisions
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*&                        !! ion-ion collisions
               (uIon_DIC(3,iIonFluid,i,j,k)-uIon_DIC(3,1:nIonFluid,i,j,k))
       end do                                                                                                         !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*&                         !! ion-neutral collisions
               (UnzNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(3,1:nIonFluid,i,j,k))
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
               (UnzNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(3,1:nIonFluid,i,j,k))
       end do

       kinAdd_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! addition to individual fluid from charge exchange [1/(m^3*s)]
                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)* &
                        NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
                        (UnzNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(3,jIonFluid,i,j,k))
                end do
             end do
          end do
       end do

       SRhoUzTerm_IIC(1,1:nIonFluid,i,j,k) = (vAdd_I(1:nIonFluid)/Si2No_V(UnitT_)+ &                                 !! newly photoionized neutrals
            kinAdd_I(1:nIonFluid)/Si2No_V(UnitT_))*Si2No_V(UnitN_)*MassIon_I* &                                      !! new ions from charge exchange
            Si2No_V(UnitU_)
       ! Add u_s*drho_s/dt for d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt
       do iIonFluid=1,nIonFluid
          SRhoUzTerm_IIC(1,iIonFluid,i,j,k) = SRhoUzTerm_IIC(1,iIonFluid,i,j,k)+sum(SRhoTerm_IIC(1:4,iIonFluid,i,j,k))&
               *uIon_DIC(3,iIonFluid,i,j,k)*Si2No_V(UnitU_)
       end do
       SRhoUzTerm_IIC(2,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)/ElectronCharge* &            !! current dissipation, ion-electron collisions
            MassIon_I*nIon_IC(1:nIonFluid,i,j,k)/nElec_C(i,j,k)*Current_DC(3,i,j,k)
       SRhoUzTerm_IIC(3,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*MassIon_I*&
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitN_)*(uIon_DIC(3,1:nIonFluid,i,j,k)-uIonMean_DC(3,i,j,k))*Si2No_V(UnitU_)
       SRhoUzTerm_IIC(4,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*fiiTot_I(1:nIonFluid)* &   !! ion-ion collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I
       SRhoUzTerm_IIC(5,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*finTot_I(1:nIonFluid)* &   !! ion neutral collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I

       ! (u_n-u_e)^2 difference in neutral and electron speeds qubed [m^2/s^2]
       do iNeutral=1,nNeutral
          uNeuElec2_IC(iNeutral,i,j,k) = ((UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uElec_DC(1,i,j,k))**2 &
               +(UnyNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uElec_DC(2,i,j,k))**2 &
               +(UnzNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uElec_DC(3,i,j,k))**2)
       end do

       ! (u_i-u_e)^2 difference in ion and electron speeds qubed [m^2/s^2]
       do iIonFluid=1,nIonFluid
          uIonElec2_IC(iIonFluid,i,j,k) = (uIon_DIC(1,iIonFluid,i,j,k)-uElec_DC(1,i,j,k))**2+&
               (uIon_DIC(2,iIonFluid,i,j,k)-uElec_DC(2,i,j,k))**2+&
               (uIon_DIC(3,iIonFluid,i,j,k)-uElec_DC(3,i,j,k))**2
       end do

       !! Second moment
       !! Sources separated into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"       
       kinSub_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! subtraction to individual fluid from charge exchange [1/(m^3*s)]
                   kinSub_I(iIonFluid) = kinSub_I(iIonFluid) + &!!nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
                end do
             end do
          end do
       end do

       vAdd_I = 0.
       do iNeutral=1,nNeutral
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)*&
               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
       end do

       SPTerm_IIC(1,1:nIonFluid,i,j,k) = -(kinSub_I(1:nIonFluid)/Si2No_V(UnitT_)+ &                                  !! lost ions through charge exchange and recombination
            alpha_IC(1:nIonFluid,i,j,k)*nElec_C(i,j,k)/Si2No_V(UnitT_))*State_VGB(iPIon_I,i,j,k,iBlock)

       fiiTot_I(1:nIonFluid) = 0.                                                                                    !! momentum transfer by ion-ion collisions
       do iIonFluid=1,nIonFluid
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*nIon_IC(1:nIonFluid,i,j,k)*&  
               MassIon_I(1:nIonFluid)/(MassIon_I(1:nIonFluid)+MassIon_I(iIonFluid))*&
               cBoltzmann*(Ti_IC(iIonFluid,i,j,k)-Ti_IC(1:nIonFluid,i,j,k))
       end do

       SPTerm_IIC(2,1:nIonFluid,i,j,k) = 2.*fiiTot_I(1:nIonFluid)/Si2No_V(UnitT_)*Si2No_V(UnitEnergyDens_)        

       finTot_I(1:nIonFluid) = 0.                                                                                    !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*nIon_IC(1:nIonFluid,i,j,k)*&  
               MassIon_I(1:nIonFluid)/(MassIon_I(1:nIonFluid)+NeutralMass_I(iNeutral)/cProtonMass)*&
               cBoltzmann*(TnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-Ti_IC(1:nIonFluid,i,j,k))
       end do

       SPTerm_IIC(3,1:nIonFluid,i,j,k) = 2.*finTot_I(1:nIonFluid)/Si2No_V(UnitT_)*Si2No_V(UnitEnergyDens_)
       SPTerm_IIC(4,1:nIonFluid,i,j,k) = 2.*fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*nIon_IC(1:nIonFluid,i,j,k)*&
            cBoltzmann*(Te_C(i,j,k)-Ti_IC(1:nIonFluid,i,j,k))*Si2No_V(UnitEnergyDens_)
       SPTerm_IIC(5,1:nIonFluid,i,j,k) = 2./3.*fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*cElectronMass*&             !! ion-electron collisional exchange (due to Hall velocity)
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*uIonElec2_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitU_)**2

       fiiTot_I(1:nIonFluid) = 0.                                                                                    !! momentum transfer by ion-ion collisions
       do iIonFluid=1,nIonFluid
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*nIon_IC(1:nIonFluid,i,j,k)*&  
               MassIon_I(1:nIonFluid)*MassIon_I(iIonFluid)/(MassIon_I(1:nIonFluid)+MassIon_I(iIonFluid))*&
               uIonIon2_IIC(1:nIonFluid,iIonFluid,i,j,k)
       end do

       finTot_I(1:nIonFluid) = 0.                                                                                    !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*nIon_IC(1:nIonFluid,i,j,k)*&  
               MassIon_I(1:nIonFluid)*NeutralMass_I(iNeutral)/(MassIon_I(1:nIonFluid)+NeutralMass_I(iNeutral)/cProtonMass)*&
               uIonNeu2_IIC(1:nIonFluid,iNeutral,i,j,k)/cProtonMass
       end do

       SPTerm_IIC(6,1:nIonFluid,i,j,k) = 2./3.*fiiTot_I(1:nIonFluid)/Si2No_V(UnitT_)*Si2No_V(UnitN_)*Si2No_V(UnitU_)**2
       SPTerm_IIC(7,1:nIonFluid,i,j,k) = 2./3.*finTot_I(1:nIonFluid)/Si2No_V(UnitT_)*Si2No_V(UnitN_)*Si2No_V(UnitU_)**2


       do iNeutral=1,nNeutral
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*uIonNeu2_IIC(1:nIonFluid,iNeutral,i,j,k)
       end do
       kinAdd_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! addition to individual fluid from charge exchange [1/(m*s^2)]
                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)*&
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)* &
                        NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*uIonNeu2_IIC(jIonFluid,iNeutral,i,j,k)
                end do
             end do
          end do
       end do

       SPTerm_IIC(8,1:nIonFluid,i,j,k) = 1./3.*(vAdd_I(1:nIonFluid)/Si2No_V(UnitT_)+kinAdd_I(1:nIonFluid)/Si2No_V(UnitT_))*&
            MassIon_I(1:nIonFluid)*Si2No_V(UnitN_)*Si2No_V(UnitU_)**2

       if (UseElectronPressure) then
          SPeTerm_IC(1,i,j,k) = -sum(alpha_IC(1:nIonFluid,i,j,k)*nIon_IC(1:nIonFluid,i,j,k))/ &                           !! lost electrons through recombination
               Si2No_V(UnitT_)*State_VGB(Pe_,i,j,k,iBlock)

          vAdd_I(1:nIonFluid) = 0.
          do iNeutral=1,nNeutral
             vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
                  NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*uNeuElec2_IC(iNeutral,i,j,k)
          end do
          SPeTerm_IC(2,i,j,k) = 1./3.*cElectronMass*sum(vAdd_I)*Si2No_V(UnitRho_)/Si2No_V(UnitT_)*Si2No_V(UnitU_)**2      !! new electrons through photoionized neutrals

          feiTot = 0.
          do iIonFluid=1,nIonFluid
             feiTot = feiTot+fei_IC(iIonFluid,i,j,k)/MassIon_I(iIonFluid)*&
                  (Ti_IC(iIonFluid,i,j,k)-Te_C(i,j,k))
          end do
          SPeTerm_IC(3,i,j,k) = 2.*cElectronMass*Si2No_V(UnitRho_)/Si2No_V(UnitN_)*&                                      !! ion-electron collisional exchange (thermal motion)
               nElec_C(i,j,k)*cBoltzmann*Si2No_V(UnitEnergyDens_)*feiTot/Si2No_V(UnitT_)

          fenTot = 0.
          do iNeutral=1,nNeutral
             fenTot = fenTot+fen_IC(iNeutral,i,j,k)/NeutralMass_I(iNeutral)*&
                  (TnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-Te_C(i,j,k))
          end do
          SPeTerm_IC(4,i,j,k) = 2.*cElectronMass*nElec_C(i,j,k)*cBoltzmann*&                                              !! electron-neutral collisional exchange (thermal motion)
               Si2No_V(UnitEnergyDens_)*fenTot/Si2No_V(UnitT_)

          SPeTerm_IC(5,i,j,k) = 2./3.*sum(fei_IC(1:nIonFluid,i,j,k)*uIonElec2_IC(1:nIonFluid,i,j,k))/ &                   !! ion-electron collisional exchange (due to Hall velocity)
               Si2No_V(UnitT_)*cElectronMass*nElec_C(i,j,k)*Si2No_V(UnitRho_)*Si2No_V(UnitU_)**2

          SPeTerm_IC(6,i,j,k) = 2./3.*sum(fen_IC(1:nNeutral,i,j,k)*uNeuElec2_IC(1:nNeutral,i,j,k))/&                      !! electron-neutral collisional exchange (bulk motion)
               Si2No_V(UnitT_)*cElectronMass*nElec_C(i,j,k)*Si2No_V(UnitRho_)*Si2No_V(UnitU_)**2

          vAdd_I(1:nIonFluid) = 0.
          do iNeutral=1,nNeutral
             vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+((v_IIC(iNeutral,1:nIonFluid,i,j,k)-&
                  ve_IIC(iNeutral,1:nIonFluid,i,j,k))*Qexc_II(iNeutral,1:nIonFluid)- &
                  ve_IIC(iNeutral,1:nIonFluid,i,j,k)*Qion_II(iNeutral,1:nIonFluid))*ChargeIon_I(1:nIonFluid)* &          
                  NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
          end do
          SPeTerm_IC(7,i,j,k) = 2./3.*sum(vAdd_I)*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)                                ! heating of electrons due to ionization excess energy

       end if



       MassLoading(1,i,j,k,iBlock) = SRhoTerm_IIC(1,Op_,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)
       ! MassLoading(2,i,j,k,iBlock) = SRhoTerm_IIC(2,Op_,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)
       MassLoading(2,i,j,k,iBlock) = v_IIC(O2_,Op_,i,j,k)
       ! MassLoading(3,i,j,k,iBlock) = -SRhoTerm_IIC(3,Op_,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)
       MassLoading(3,i,j,k,iBlock) = ve_IIC(O2_,Op_,i,j,k)
       MassLoading(4,i,j,k,iBlock) = -SRhoTerm_IIC(4,Op_,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)
       MassLoading(5,i,j,k,iBlock) = SRhoTerm_IIC(1,O2p_,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)
       ! MassLoading(6,i,j,k,iBlock) = SRhoTerm_IIC(2,O2p_,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)
       MassLoading(6,i,j,k,iBlock) = v_IIC(O2_,O2p_,i,j,k)
       ! MassLoading(7,i,j,k,iBlock) = -SRhoTerm_IIC(3,O2p_,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)
       MassLoading(7,i,j,k,iBlock) = ve_IIC(O2_,O2p_,i,j,k)
       MassLoading(8,i,j,k,iBlock) = -SRhoTerm_IIC(4,O2p_,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)


       !! sum up individual terms
       do iTerm=1,4
          SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(iTerm,1:nIonFluid,i,j,k)
       end do
       ! SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(4,1:nIonFluid,i,j,k)
       do iTerm=1,5
          SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(iTerm,1:nIonFluid,i,j,k)
          SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(iTerm,1:nIonFluid,i,j,k)
          SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(iTerm,1:nIonFluid,i,j,k)
       end do
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(4,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(4,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(4,1:nIonFluid,i,j,k)
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(5,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(5,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(5,1:nIonFluid,i,j,k)
       do iTerm=1,8
          SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(iTerm,1:nIonFluid,i,j,k)
       end do
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(4,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(5,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(6,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(7,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(8,1:nIonFluid,i,j,k)
       if(UseElectronPressure) then
          SPe_C(i,j,k) = sum(SPeTerm_IC(1:7,i,j,k))
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(1,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(2,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(3,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(4,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(5,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(6,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(7,i,j,k)
       end if

       Source_VC(iRhoIon_I   ,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)    + Source_VC(iRhoIon_I   ,i,j,k)
       Source_VC(iRhoUxIon_I ,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)  + Source_VC(iRhoUxIon_I ,i,j,k)
       Source_VC(iRhoUyIon_I ,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)  + Source_VC(iRhoUyIon_I ,i,j,k)
       Source_VC(iRhoUzIon_I ,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)  + Source_VC(iRhoUzIon_I ,i,j,k)
       Source_VC(iPIon_I     ,i,j,k) = SP_IC(1:nIonFluid,i,j,k)      + Source_VC(iPIon_I     ,i,j,k)

       Source_VC(Rho_   ,i,j,k) = sum(SRho_IC(1:nIonFluid,i,j,k))    + Source_VC(Rho_   ,i,j,k)
       Source_VC(rhoUx_ ,i,j,k) = sum(SRhoUx_IC(1:nIonFluid,i,j,k))  + Source_VC(rhoUx_ ,i,j,k)
       Source_VC(rhoUy_ ,i,j,k) = sum(SRhoUy_IC(1:nIonFluid,i,j,k))  + Source_VC(rhoUy_ ,i,j,k)
       Source_VC(rhoUz_ ,i,j,k) = sum(SRhoUz_IC(1:nIonFluid,i,j,k))  + Source_VC(rhoUz_ ,i,j,k)
       Source_VC(Bx_    ,i,j,k) = SBx_C(i,j,k)                       + Source_VC(Bx_    ,i,j,k)
       Source_VC(By_    ,i,j,k) = SBy_C(i,j,k)                       + Source_VC(By_    ,i,j,k)
       Source_VC(Bz_    ,i,j,k) = SBz_C(i,j,k)                       + Source_VC(Bz_    ,i,j,k)
       if(UseElectronPressure) then
          Source_VC(P_     ,i,j,k) = sum(SP_IC(1:nIonFluid,i,j,k))   + Source_VC(P_     ,i,j,k)
          Source_VC(Pe_    ,i,j,k) = SPe_C(i,j,k)                    + Source_VC(Pe_    ,i,j,k)
       else
          Source_VC(P_     ,i,j,k) = sum(SP_IC(1:nIonFluid,i,j,k))*(1.+ElectronPressureRatio) + &
               Source_VC(P_     ,i,j,k)
       end if

    end do;  end do;  end do

    if(DoTestMe) then
       write(*,*)'user_calc_sources:'
       write(*,*)'Inputs: '
       i=iTest ; j=jTest ; k=kTest
       theta=acos((-SW_Ux*Xyz_DGB(x_,i,j,k,iBlock)-SW_Uy*Xyz_DGB(y_,i,j,k,iBlock)&
            -SW_Uz*Xyz_DGB(z_,i,j,k,iBlock))/R_BLK(i,j,k,iBlock)/&
            (SW_Ux**2+SW_Uy**2+SW_Uz**2)**0.5)
123    format (A13,ES25.16,A15,A3,F7.2,A3)
       write(*,123)'x         = ',Xyz_DGB(x_,i,j,k,iBlock)," [rPlanet]"
       write(*,123)'y         = ',Xyz_DGB(y_,i,j,k,iBlock)," [rPlanet]"
       write(*,123)'z         = ',Xyz_DGB(z_,i,j,k,iBlock)," [rPlanet]"
       write(*,123)'r         = ',R_BLK(i,j,k,iBlock)," [rPlanet]"
       write(*,123)'SW_Ux     = ',SW_Ux*No2SI_V(UnitU_)," [m/s]"
       write(*,123)'SW_Uy     = ',SW_Uy*No2SI_V(UnitU_)," [m/s]"
       write(*,123)'SW_Uz     = ',SW_Uz*No2SI_V(UnitU_)," [m/s]"
       write(*,123)'Tmin      = ',Tmin," [K]"
       write(*,*)''
       write(*,*)'Neutrals:'
       do iNeutral=1,nNeutral
          write(*,124)'Neutral species #',iNeutral,': ', NameNeutral_I(iNeutral)," (",&
               NeutralMass_I(iNeutral)/cProtonMass," amu)"
          write(*,123)'n_n       = ',NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)," [m^-3]"
          write(*,123)'m_n       = ',NeutralMass_I(iNeutral)," [kg]"
          write(*,123)'unx       = ',UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)," [m/s]"
          write(*,123)'uny       = ',UnyNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)," [m/s]"
          write(*,123)'unz       = ',UnzNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)," [m/s]"
          write(*,123)'Tn        = ',TnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)," [K]"
       end do
       write(*,*)''
       write(*,*)'Total plasma phase (e- and i+):'
       write(*,123)'Rho       = ',State_VGB(Rho_,i,j,k,iBlock)*No2SI_V(UnitRho_)," [kg/m^3]"
       write(*,123)'uRhox     = ',State_VGB(RhoUx_,i,j,k,iBlock)*No2SI_V(UnitRhoU_)," [kg/(m^2*s)]"
       write(*,123)'uRhoy     = ',State_VGB(RhoUy_,i,j,k,iBlock)*No2SI_V(UnitRhoU_)," [kg/(m^2*s)]"
       write(*,123)'uRhoz     = ',State_VGB(RhoUz_,i,j,k,iBlock)*No2SI_V(UnitRhoU_)," [kg/(m^2*s)]"
       if (UseElectronPressure) then
          write(*,123)'Ptot      = ',(State_VGB(P_,i,j,k,iBlock)+State_VGB(Pe_,i,j,k,iBlock))*&
               No2SI_V(UnitP_)," [kg/(m*s^2)]"
       else
          write(*,123)'Ptot      = ',State_VGB(P_,i,j,k,iBlock)*No2SI_V(UnitP_)," [kg/(m*s^2)]"
       end if
       write(*,123)'Bx        = ',State_VGB(Bx_,i,j,k,iBlock)*No2SI_V(UnitB_)," [T]"
       write(*,123)'By        = ',State_VGB(By_,i,j,k,iBlock)*No2SI_V(UnitB_)," [T]"
       write(*,123)'Bz        = ',State_VGB(Bz_,i,j,k,iBlock)*No2SI_V(UnitB_)," [T]"
       write(*,123)'uMeanx    = ',uIonMean_DC(1,i,j,k)," [m/s]"
       write(*,123)'uMeany    = ',uIonMean_DC(2,i,j,k)," [m/s]"
       write(*,123)'uMeanz    = ',uIonMean_DC(3,i,j,k)," [m/s]"       
       write(*,123)'jx        = ',Current_DC(1,i,j,k)*No2SI_V(UnitJ_)," [A/m^2]"
       write(*,123)'jy        = ',Current_DC(2,i,j,k)*No2SI_V(UnitJ_)," [A/m^2]"
       write(*,123)'jz        = ',Current_DC(3,i,j,k)*No2SI_V(UnitJ_)," [A/m^2]"
       write(*,*)''
       write(*,123)'SRho      = ',sum(SRho_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
            100.*sum(SRho_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/(State_VGB(Rho_,i,j,k,iBlock)),"%)"
       if (State_VGB(RhoUx_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SRhoUx    = ',sum(SRhoUx_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
               " (", 100.*sum(SRhoUx_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/(State_VGB(RhoUx_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SRhoUx    = ',sum(SRhoUx_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
       end if
       if (State_VGB(RhoUy_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SRhoUy    = ',sum(SRhoUy_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
               " (", 100.*sum(SRhoUy_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/(State_VGB(RhoUy_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SRhoUy    = ',sum(SRhoUy_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
       end if
       if (State_VGB(RhoUz_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SRhoUz    = ',sum(SRhoUz_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
               " (", 100.*sum(SRhoUz_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/(State_VGB(RhoUz_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SRhoUz    = ',sum(SRhoUz_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
       end if
       if (State_VGB(Bx_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SBx       = ',SBx_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"," (", &
               100.*SBx_C(i,j,k)*Dt_BLK(iBlock)/(State_VGB(Bx_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SBx       = ',SBx_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"
       end if
       if (State_VGB(By_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SBy       = ',SBy_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"," (", &
               100.*SBy_C(i,j,k)*Dt_BLK(iBlock)/(State_VGB(By_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SBy       = ',SBy_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"
       end if
       if (State_VGB(Bz_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SBz       = ',SBz_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"," (", &
               100.*SBz_C(i,j,k)*Dt_BLK(iBlock)/(State_VGB(Bz_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SBz       = ',SBz_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"
       end if
       if(UseElectronPressure) then
          write(*,123)'SP        = ',sum(SP_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*sum(SP_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/(State_VGB(P_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SP        = ',sum(SP_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitP_)/No2SI_V(UnitT_)*&
               (1+ElectronPressureRatio)," [Pa/s]"," (",100.*sum(SP_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/ &
               (State_VGB(P_,i,j,k,iBlock))*(1+ElectronPressureRatio),"%)"
       end if
       write(*,*)''
       write(*,123)'dt        = ',Dt_BLK(iBlock)*No2SI_V(UnitT_)," [s]"
       write(*,*)''
       write(*,*)'Individual ion fluids:'
       do iIonFluid=1,nIonFluid
          write(*,124)'Ion species     #',iIonFluid,': ',NameFluid_I(iIonFluid+1)," (",&
               MassIon_I(iIonFluid)," amu/",ChargeIon_I(iIonFluid)," e)"
124       format (A17,I2,A3,A7,A3,F5.2,A5,F5.1,A3)
          write(*,123)'Ux        = ',uIon_DIC(1,iIonFluid,i,j,k)," [m/s]"
          write(*,123)'Uy        = ',uIon_DIC(2,iIonFluid,i,j,k)," [m/s]"
          write(*,123)'Uz        = ',uIon_DIC(3,iIonFluid,i,j,k)," [m/s]"
          write(*,123)'ni        = ',nIon_IC(iIonFluid,i,j,k)," [m^-3]"
          write(*,123)'Ti        = ',Ti_IC(iIonFluid,i,j,k)," [K]"
          write(*,123)'Rho       = ',State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRho_)," [kg/m^3]"
          write(*,123)'rhoUx     = ',State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'rhoUy     = ',State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'rhoUz     = ',State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'Pi        = ',State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitP_)," [Pa]"
          write(*,123)'SRho      = ',SRho_IC(iIonFluid,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
               100.*SRho_IC(iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SRhoT1   = ',SRhoTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
               100.*SRhoTerm_IIC(1,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SRhoT2   = ',SRhoTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
               100.*SRhoTerm_IIC(2,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SRhoT3   = ',SRhoTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
               100.*SRhoTerm_IIC(3,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SRhoT4   = ',SRhoTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
               100.*SRhoTerm_IIC(4,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          if (State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock) /= 0.) then
             write(*,123)'SRhoUx    = ',SRhoUx_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUx_IC(iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUxT1 = ',SRhoUxTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUxTerm_IIC(1,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUxT2 = ',SRhoUxTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUxTerm_IIC(2,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUxT3 = ',SRhoUxTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUxTerm_IIC(3,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUxT4 = ',SRhoUxTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUxTerm_IIC(4,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUxT5 = ',SRhoUxTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUxTerm_IIC(5,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          else
             write(*,123)'SRhoUx    = ',SRhoUx_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUxT1 = ',SRhoUxTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUxT2 = ',SRhoUxTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUxT3 = ',SRhoUxTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUxT4 = ',SRhoUxTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUxT5 = ',SRhoUxTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
          end if
          if (State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock) /= 0.) then
             write(*,123)'SRhoUy    = ',SRhoUy_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUy_IC(iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUyT1 = ',SRhoUyTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUyTerm_IIC(1,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUyT2 = ',SRhoUyTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUyTerm_IIC(2,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUyT3 = ',SRhoUyTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUyTerm_IIC(3,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUyT4 = ',SRhoUyTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUyTerm_IIC(4,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUyT5 = ',SRhoUyTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUyTerm_IIC(5,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          else
             write(*,123)'SRhoUy    = ',SRhoUy_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUyT1 = ',SRhoUyTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUyT2 = ',SRhoUyTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUyT3 = ',SRhoUyTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUyT4 = ',SRhoUyTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUyT5 = ',SRhoUyTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
          end if
          if (State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock) /= 0.) then
             write(*,123)'SRhoUz    = ',SRhoUz_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUz_IC(iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUzT1 = ',SRhoUzTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUzTerm_IIC(1,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUzT2 = ',SRhoUzTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUzTerm_IIC(2,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUzT3 = ',SRhoUzTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUzTerm_IIC(3,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUzT4 = ',SRhoUzTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUzTerm_IIC(4,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUzT5 = ',SRhoUzTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUzTerm_IIC(5,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          else
             write(*,123)'SRhoUz    = ',SRhoUz_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUzT1 = ',SRhoUzTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUzT2 = ',SRhoUzTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUzT3 = ',SRhoUzTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUzT4 = ',SRhoUzTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUzT5 = ',SRhoUzTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
          end if
          write(*,123)'SP        = ',SP_IC(iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SP_IC(iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT1     = ',SPTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(1,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT2     = ',SPTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(2,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT3     = ',SPTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(3,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT4     = ',SPTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(4,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT5     = ',SPTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(5,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT6     = ',SPTerm_IIC(6,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(6,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT7     = ',SPTerm_IIC(7,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(7,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT8     = ',SPTerm_IIC(8,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(8,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
       end do
       write(*,*)''
       write(*,*)'Electrons:'
       write(*,123)'n_e       = ',nElec_C(i,j,k)," [m^-3]"
       if (UseElectronPressure) then
          write(*,123)'Pe        = ',State_VGB(Pe_,i,j,k,iBlock)*No2SI_V(UnitP_)," [Pa]"
       else
          write(*,123)'Pe        = ',State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio/&
               (1.+ElectronPressureRatio)*No2SI_V(UnitP_)," [Pa]"
       end if
       write(*,123)'Uex       = ',uElec_DC(1,i,j,k)," [m/s]"
       write(*,123)'Uey       = ',uElec_DC(2,i,j,k)," [m/s]"
       write(*,123)'Uez       = ',uElec_DC(3,i,j,k)," [m/s]"
       write(*,123)'Te        = ',Te_C(i,j,k)," [K]"
       if(UseElectronPressure) then
          if (State_VGB(Pe_,i,j,k,iBlock).gt.0.) then
             write(*,123)'SPe       = ',SPe_C(i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPe_C(i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT1    = ',SPeTerm_IC(1,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(1,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT2    = ',SPeTerm_IC(2,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(2,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT3    = ',SPeTerm_IC(3,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(3,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT4    = ',SPeTerm_IC(4,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(4,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT5    = ',SPeTerm_IC(5,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(5,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT6    = ',SPeTerm_IC(6,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(6,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT7    = ',SPeTerm_IC(7,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(7,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
          else
             write(*,123)'SPe       = ',SPe_C(i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT1    = ',SPeTerm_IC(1,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT2    = ',SPeTerm_IC(2,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT3    = ',SPeTerm_IC(3,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT4    = ',SPeTerm_IC(4,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT5    = ',SPeTerm_IC(5,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT6    = ',SPeTerm_IC(6,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT7    = ',SPeTerm_IC(7,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
          end if
       end if
       write(*,*)''
       write(*,*)'Ion-electron combinations:'
       do iIonFluid=1,nIonFluid
          write(*,*)NameFluid_I(iIonFluid+1), '&  e'
          write(*,123)'fei       = ',fei_IC(iIonFluid,i,j,k)," [1/s]"
          write(*,123)'fie       = ',fie_IC(iIonFluid,i,j,k)," [1/s]"
          write(*,123)'|u_ime|   = ',sqrt(uIonElec2_IC(iIonFluid,i,j,k))," [m/s]"
          write(*,123)'alpha     = ',alpha_IC(iIonFluid,i,j,k)," [m^3/s]"
       end do
       write(*,*)''
       write(*,*)'Ion-ion combinations:'
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             if (iIonFluid.ne.jIonFluid) then
                write(*,*)NameFluid_I(iIonFluid+1), '&  ',NameFluid_I(jIonFluid+1)
                write(*,123)' fii      = ',fii_IIC(iIonFluid,jIonFluid,i,j,k)," [1/s]"
                write(*,123)' |u_imi|  = ',sqrt(uIonIon2_IIC(iIonFluid,jIonFluid,i,j,k))," [m/s]"
             end if
          end do
       end do
       write(*,*)''
       write(*,*)'Ion-neutral combinations:'
       do iIonFluid=1,nIonFluid
          do iNeutral=1,nNeutral
             write(*,*)NameFluid_I(iIonFluid+1), '&  ',NameNeutral_I(iNeutral)
             write(*,123)' v_phio   = ',v_IIC(iNeutral,iIonFluid,i,j,k)-ve_IIC(iNeutral,iIonFluid,i,j,k)," [1/s]"
             write(*,123)' v_eio    = ',ve_IIC(iNeutral,iIonFluid,i,j,k)," [1/s]"
             write(*,123)' v_io     = ',v_IIC(iNeutral,iIonFluid,i,j,k)," [1/s]"
             write(*,123)' fin      = ',fin_IIC(iIonFluid,iNeutral,i,j,k)," [1/s]"
             write(*,123)' |u_imn|  = ',sqrt(uIonNeu2_IIC(iIonFluid,iNeutral,i,j,k))," [m/s]"
             write(*,*)' kin (Ion & Neutral-> Neutral & Ion):'
             do jIonFluid=1,nIonFluid
                do jNeutral=1,nNeutral
                   write(*,*)' ',NameFluid_I(iIonFluid+1),'&  ',NameNeutral_I(iNeutral),'->  ', &
                        NameNeutral_I(jNeutral),'&  ',NameFluid_I(jIonFluid+1),'=',&
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)," [m^3/s]"  
                end do
             end do
          end do
       end do
       write(*,*)''
       write(*,*)'Electron-neutral combinations:'
       do iNeutral=1,nNeutral
          write(*,*)'e & ',NameNeutral_I(iNeutral)
          write(*,123)'fen      = ',fen_IC(iNeutral,i,j,k)," [1/s]"
          write(*,123)'|u_nme|  = ',sqrt(uNeuElec2_IC(iNeutral,i,j,k))," [m/s]"
       end do
       write(*,*)''
    end if

  end subroutine user_calc_sources

  !========================================================================

  subroutine user_update_states(iStage,iBlock)
    use ModAdvance,  ONLY: State_VGB
    use ModPhysics,  ONLY: SW_N, LowDensityRatio, cBoltzmann, ElectronPressureRatio, Si2No_V, &
         No2Si_V, UnitN_, UnitP_, rBody, BodyRho_I, BodyP_I
    use ModEnergy,   ONLY: calc_energy_cell
    use ModGeometry, ONLY: r_BLK, Rmin_BLK, Xyz_DGB
    use ModMain,     ONLY: Body1

    integer,intent(in) :: iStage, iBlock
    integer :: i,j,k,iIonFluid
    real :: r_D(3), dRhoUr_D(3), RhoUr

    real, dimension(1:nI,1:nJ,1:nK) :: nElec_C
    real, dimension(1:nIonFluid,1:nI,1:nJ,1:nK) ::nIon_IC

    !----------------------------------------------------------------------

    call update_states_MHD(iStage,iBlock)

    ! Enforce minimum temperature (pressure), Tmin, if temperatures Ti_IC or Te_C are below

    do k=1,nK; do j=1,nJ; do i=1,nI
       if((R_BLK(i,j,k,iBlock) <= PlanetRadius).and.(UseResistivePlanet)) CYCLE
       do iIonFluid=1,nIonFluid
          ! set minimum mass density (and in these locations Ti = Tmin and vi=vbulkplasma)
          if(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) < SW_n*MassIon_I(iIonFluid)*LowDensityRatio**2) then
             State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) = SW_n*MassIon_I(iIonFluid)*LowDensityRatio**2
             State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) * &
                  State_VGB(RhoUx_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) * &
                  State_VGB(RhoUy_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) * &
                  State_VGB(RhoUz_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)/ &
                  MassIon_I(iIonFluid)*No2SI_V(UnitN_)*cBoltzmann*Tmin*SI2No_V(UnitP_)
          end if
       end do

       State_VGB(Rho_,i,j,k,iBlock) = sum(State_VGB(iRhoIon_I,i,j,k,iBlock))

       nIon_IC(1:nIonFluid,i,j,k) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
       nElec_C(i,j,k) = sum(nIon_IC(1:nIonFluid,i,j,k)*ChargeIon_I(1:nIonFluid))

       do iIonFluid=1,nIonFluid
          ! set minimum pressure
          if(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)*NO2SI_V(UnitP_) < &
               nIon_IC(iIonFluid,i,j,k)*cBoltzmann*Tmin) then
             State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock) = &
                  nIon_IC(iIonFluid,i,j,k)*cBoltzmann*Tmin*SI2No_V(UnitP_)
          end if
       end do

       if(UseElectronPressure) then
          State_VGB(P_,i,j,k,iBlock) = sum(State_VGB(iPIon_I,i,j,k,iBlock))
          if (State_VGB(Pe_,i,j,k,iBlock)*NO2SI_V(UnitP_) < nElec_C(i,j,k)*cBoltzmann*Tmin) then
             State_VGB(Pe_,i,j,k,iBlock) = nElec_C(i,j,k)*cBoltzmann*Tmin*SI2No_V(UnitP_)
          end if
       else
          State_VGB(P_,i,j,k,iBlock) = sum(State_VGB(iPIon_I,i,j,k,iBlock))*(1.+ElectronPressureRatio)
       end if
    end do; end do; end do

    if((Rmin_BLK(iBlock) <= PlanetRadius).and.(UseResistivePlanet)) then
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(R_BLK(i+2,j,k,iBlock) >= PlanetRadius) CYCLE
          if((Body1).and.(R_BLK(i,j,k,iBlock)<rBody)) CYCLE
          State_VGB(iRhoIon_I,i,j,k,iBlock) = PlanetDensity/nIonFluid
          State_VGB(iRhoUxIon_I,i,j,k,iBlock) = 0.0
          State_VGB(iRhoUyIon_I,i,j,k,iBlock) = 0.0
          State_VGB(iRhoUzIon_I,i,j,k,iBlock) = 0.0
          State_VGB(iPIon_I,i,j,k,iBlock) = PlanetPressure/nIonFluid          
          State_VGB(Rho_,i,j,k,iBlock) = PlanetDensity
          State_VGB(P_,i,j,k,iBlock) = PlanetPressure
          State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock) = 0.0
       end do; end do; end do



       do k=MinK,MaxK; do j=MinJ,MaxJ; do i=1,nI
          if(R_BLK(i-1,j,k,iBlock) > PlanetRadius ) CYCLE
          if(R_BLK(i,j,k,iBlock) <= PlanetRadius ) CYCLE

          r_D  = (/ Xyz_DGB(x_,i,j,k,iBlock), Xyz_DGB(y_,i,j,k,iBlock), &
               Xyz_DGB(z_,i,j,k,iBlock) /) / r_BLK(i,j,k,iBlock)

          do iIonFluid=1,nIonFluid
             RhoUr = dot_product(State_VGB(iRhoUxIon_I(iIonFluid):iRhoUzIon_I(iIonFluid),i,j,k,iBlock),r_D)
             if(RhoUr > 0.0) then
                ! ! If flow is out of the planet, remove the radial componet 
                ! ! of the momentum so that the flow is tangential
                ! dRhoUr_D = -r_D*RhoUr

                ! ! Two cell boundary layer inside the planet
                ! !! -1
                ! State_VGB(iRhoIon_I(iIonFluid),i-1,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)
                ! State_VGB(iRhoUxIon_I(iIonFluid):iRhoUzIon_I(iIonFluid),i-1,j,k,iBlock) = &
                !      State_VGB(iRhoUxIon_I(iIonFluid):iRhoUzIon_I(iIonFluid),i,j,k,iBlock) + dRhoUr_D

                ! !! -2
                ! State_VGB(iRhoIon_I(iIonFluid),i-2,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)
                ! State_VGB(iRhoUxIon_I(iIonFluid):iRhoUzIon_I(iIonFluid),i-2,j,k,iBlock) = &
                !      State_VGB(iRhoUxIon_I(iIonFluid):iRhoUzIon_I(iIonFluid),i,j,k,iBlock) + dRhoUr_D

                !! No flow
                ! Two cell boundary layer inside the planet
                !! -1
                ! State_VGB(iRhoIon_I(iIonFluid),i-1,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)
                State_VGB(iRhoIon_I(iIonFluid),i-1,j,k,iBlock) = BodyRho_I(iIonFluid)
                State_VGB(iRhoUxIon_I(iIonFluid):iRhoUzIon_I(iIonFluid),i-1,j,k,iBlock) = 0.0
                State_VGB(iPIon_I(iIonFluid),i-1,j,k,iBlock) = BodyP_I(iIonFluid)

                !! -2
                ! State_VGB(iRhoIon_I(iIonFluid),i-2,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)
                State_VGB(iRhoIon_I(iIonFluid),i-2,j,k,iBlock) = BodyRho_I(iIonFluid)
                State_VGB(iRhoUxIon_I(iIonFluid):iRhoUzIon_I(iIonFluid),i-2,j,k,iBlock) = 0.0
                State_VGB(iPIon_I(iIonFluid),i-2,j,k,iBlock) = BodyP_I(iIonFluid)
             else
                ! If flow is into the planet do nothing so the flow is absorbed
                dRhoUr_D = (/0.0,0.0,0.0/)

                ! Two cell boundary layer inside the planet
                !! -1
                State_VGB(iRhoIon_I(iIonFluid),i-1,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)
                State_VGB(iRhoUxIon_I(iIonFluid):iRhoUzIon_I(iIonFluid),i-1,j,k,iBlock) = &
                     State_VGB(iRhoUxIon_I(iIonFluid):iRhoUzIon_I(iIonFluid),i,j,k,iBlock) + dRhoUr_D
                State_VGB(iPIon_I(iIonFluid),i-1,j,k,iBlock) = State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)

                !! -2
                State_VGB(iRhoIon_I(iIonFluid),i-2,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)
                State_VGB(iRhoUxIon_I(iIonFluid):iRhoUzIon_I(iIonFluid),i-2,j,k,iBlock) = &
                     State_VGB(iRhoUxIon_I(iIonFluid):iRhoUzIon_I(iIonFluid),i,j,k,iBlock) + dRhoUr_D
                State_VGB(iPIon_I(iIonFluid),i-2,j,k,iBlock) = State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)
             end if

          end do
          State_VGB(iPIon_I,i-1,j,k,iBlock) = State_VGB(iPIon_I,i,j,k,iBlock)
          State_VGB(iPIon_I,i-2,j,k,iBlock) = State_VGB(iPIon_I,i,j,k,iBlock)

          State_VGB(Rho_,i-1,j,k,iBlock)   = sum(State_VGB(iRhoIon_I,i-1,j,k,iBlock))
          State_VGB(RhoUx_,i-1,j,k,iBlock) = sum(State_VGB(iRhoUxIon_I,i-1,j,k,iBlock))
          State_VGB(RhoUy_,i-1,j,k,iBlock) = sum(State_VGB(iRhoUyIon_I,i-1,j,k,iBlock))
          State_VGB(RhoUz_,i-1,j,k,iBlock) = sum(State_VGB(iRhoUzIon_I,i-1,j,k,iBlock))
          State_VGB(P_,i-1,j,k,iBlock)     = sum(State_VGB(iPIon_I,i-1,j,k,iBlock))

          State_VGB(Rho_,i-2,j,k,iBlock)   = sum(State_VGB(iRhoIon_I,i-2,j,k,iBlock))
          State_VGB(RhoUx_,i-2,j,k,iBlock) = sum(State_VGB(iRhoUxIon_I,i-2,j,k,iBlock))
          State_VGB(RhoUy_,i-2,j,k,iBlock) = sum(State_VGB(iRhoUyIon_I,i-2,j,k,iBlock))
          State_VGB(RhoUz_,i-2,j,k,iBlock) = sum(State_VGB(iRhoUzIon_I,i-2,j,k,iBlock))
          State_VGB(P_,i-2,j,k,iBlock)     = sum(State_VGB(iPIon_I,i-2,j,k,iBlock))

       end do; end do; end do

    end if

    call calc_energy_cell(iBlock)

  end subroutine user_update_states

  !========================================================================

  subroutine derive_cell_diffusivity(iBlock, i, j, k, TeSI, nIon_I, nElec, EtaSi)
    use ModResistivity,  ONLY: Eta0SI
    use ModAdvance,      ONLY: State_VGB
    use ModVarIndexes,   ONLY: Rho_, Bx_, Bz_
    use ModPhysics,      ONLY: No2Si_V, UnitP_, UnitN_, UnitTemperature_, ElectronPressureRatio, &
         UnitB_, UnitRho_
    use ModConst,        ONLY: cBoltzmann, cElectronMass, cElectronCharge, &
         cMu, cProtonMass, cPi, cEps!, cLightSpeed, cTwoPi
    use ModB0,           ONLY: B0_DGB
    use ModMain,         ONLY: UseB0, iTest, jTest, kTest, BlkTest, ProcTest
    use ModProcMH,       ONLY: iProc

    integer, intent(in)  :: iBlock, i, j, k
    real,    intent(in)  :: TeSI
    real,    intent(in)  :: nIon_I(1:nIonFluid)
    real,    intent(in)  :: nElec
    real,    intent(out) :: EtaSi

    real :: EtaSiColl!, EtaSiSpitzer, lnL
    !    real, save :: SpitzerCoef, EtaPerpSpitzerSi
    logical, save :: FirstCall = .true.
    logical :: DoTest, DoTestMe=.true.
    real :: eeSigma!, B0_D(3)
    real, dimension(nIonFluid) :: fei_I, eiSigma_I
    real, dimension(nNeutral)  :: fen_I, enSigma_I
    integer :: iIonFluid, iNeutral

    !----------------------------------------------------------------------

    if(iBlock==BlkTest.and.i==iTest.and.j==jTest.and.k==kTest.and.iProc==ProcTest) then
       call set_oktest('derive_cell_diffusivity',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    ! Spitzer formulation from Stoecker "Taschenbuch der Physik", Verlag "Harri Deutsch"
    ! lnL = log(1e7*TeSI**1.5/sqrt(nElec))
    ! EtaSiSpitzer = cElectronCharge**2*lnL/(32.*sqrt(2*cPi/cElectronMass*(cBoltzmann*TeSI)**3)*cEps**2)/cMu


    !! Collisional type resisitivity/diffusivity
    call calc_electron_collision_rates(TeSI,nElec,i,j,k,iBlock,fen_I(1:nNeutral),fei_I(1:nIonFluid))
    eiSigma_I(1:nIonFluid) = cElectronCharge**2*nElec/((fei_I(1:nIonFluid)+1E-20)*cElectronMass) 
    enSigma_I(1:nNeutral) = cElectronCharge**2*nElec/((fen_I(1:nNeutral)+1E-20)*cElectronMass)
    !! Eta_G is calculated from both conductivities using Kirchhoff's rule:
    !! 1/sigma_tot = 1/eiSigma_I+1/enSigma_I
    !! The resulting conductivity is close to Spitzer conductivity far from the comet and
    !! decreases due to abundant electron-neutral collisions close to the nucleus
    !! EtaSiColl = 1/(sigma_tot*mu_0) magnetic diffusivity [m^2/s]
    EtaSiColl = (sum(1/eiSigma_I(1:nIonFluid))+sum(1/enSigma_I(1:nNeutral)))/cMu
    !! Total diffusivity [m^2/s]
    EtaSi = Eta0SI + EtaSiColl

    ! TestArray(1,i,j,k,iBlock) = EtaSiColl
    ! TestArray(2,i,j,k,iBlock) = EtaSiSpitzer
    ! TestArray(3,i,j,k,iBlock) = EtaSiSpitzer/EtaSiColl

    if(DoTestMe) then
       write(*,*)'derive_cell_diffusivity:'
       write(*,*)'n_e    = ',nElec," [1/m^3]"
       write(*,*)'Te     = ',TeSI," [K]"
       do iIonFluid=1,nIonFluid
          write(*,*)'e & ',NameFluid_I(iIonFluid+1),':'
          write(*,*)'s_ei  = ',eiSigma_I(iIonFluid)," [1/(Ohm*m)]"
       end do
       do iNeutral=1,nNeutral
          write(*,*)'e & ',NameNeutral_I(iNeutral),':'
          write(*,*)'s_en  = ',enSigma_I(iNeutral)," [1/(Ohm*m)]"
       end do
       write(*,*)'e & e:'
       write(*,*)'s_ee  = ',eeSigma," [1/(Ohm*m)]"
       write(*,*)''
       write(*,*)'Eta0   = ',Eta0Si," [m^2/s]"
       write(*,*)'Eta_en = ',sum(1/enSigma_I(1:nNeutral))/cMu," [m^2/s]"
       write(*,*)'Eta_ei = ',sum(1/eiSigma_I(1:nIonFluid))/cMu," [m^2/s]"
       write(*,*)'Eta_ee = ',1/eeSigma/cMu," [m^2/s]"
       write(*,*)'Eta_eX = ',EtaSiColl," [m^2/s]"
       write(*,*)'EtaTot = ',EtaSi," [m^2/s]"
       write(*,*)''
    end if

  end subroutine derive_cell_diffusivity

  !========================================================================

  subroutine user_set_resistivity(iBlock, Eta_G)

    use ModPhysics,     ONLY: No2Io_V, Io2No_V, No2Si_V, Si2No_V, rBody, &
         UnitN_, UnitTemperature_, UnitX_, UnitT_, UnitP_, ElectronPressureRatio
    use ModProcMH,      ONLY: iProc
    use ModMain,        ONLY: ProcTest, BlkTest, iTest, jTest, kTest, nBlockMax, Body1
    use ModAdvance,     ONLY: State_VGB
    use ModGeometry,    ONLY: Rmin_BLK, R_BLK
    use ModVarIndexes,  ONLY: Rho_, Pe_, P_
    use ModConst,       ONLY: cMu, cBoltzmann, cElectronMass, cElectronCharge
    use ModMultiFluid,  ONLY: MassIon_I
    use ModResistivity, ONLY: Eta0

    integer, intent(in) :: iBlock
    real, intent(out) :: Eta_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK) 

    integer :: i, j, k
    logical :: DoTest, DoTestMe=.true.
    real, dimension(1:nNeutral,MinI:MaxI,MinJ:MaxJ,MinK:MaxK) :: enSigma_IG
    real, dimension(1:nIonFluid,MinI:MaxI,MinJ:MaxJ,MinK:MaxK) :: eiSigma_IG, nIon_IG
    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK) :: Te_G, nElec_G
    real, dimension(1:nNeutral) :: fen_I
    real, dimension(1:nIonFluid) :: fei_I
    integer :: iIonFluid, iNeutral
    real :: EtaSi

    !---------------------------------------------------------------------

    if(iProc==PROCtest .and. iBlock == BlkTest) then
       call set_oktest('user_set_resistivity',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    if (iBlock.ne.iNeutralBlockLast) then
       call user_neutral_atmosphere(iBlock)
    end if

    ! nElec_G is the electron/ion density in SI units (n_e=n_itot)
    nElec_G = 0.
    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
       nIon_IG(1:nIonFluid,i,j,k) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*NO2SI_V(UnitN_)
       nElec_G(i,j,k) = sum(nIon_IG(1:nIonFluid,i,j,k)*ChargeIon_I(1:nIonFluid))
    end do; end do; end do

    if (UseElectronPressure) then
       Te_G = State_VGB(Pe_,:,:,:,iBlock)*NO2SI_V(UnitP_)/(cBoltzmann* &
            nElec_G)
    else
       Te_G(:,:,:) = State_VGB(P_,:,:,:,iBlock)*ElectronPressureRatio/(1.+ElectronPressureRatio)*&
            NO2SI_V(UnitP_)/(cBoltzmann*nElec_G(:,:,:))
    end if

    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
       if((Body1).and.(r_BLK(i,j,k,iBlock) < rBody)) CYCLE
       if((UseResistivePlanet).and.(r_BLK(i,j,k,iBlock) < PlanetRadius)) CYCLE

       call derive_cell_diffusivity(iBlock, i, j, k, Te_G(i,j,k), nIon_IG(1:nIonFluid,i,j,k), nElec_G(i,j,k), EtaSi)
       Eta_G(i,j,k) = EtaSi*SI2No_V(UnitX_)**2/SI2No_V(UnitT_)

    end do; end do; end do

    ! Resistivity inside the body
    if (UseResistivePlanet) then
       do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
          do iLayer=nLayer-1,1,-1
             if(R_BLK(i,j,k,iBlock) < PlanetRadius_I(iLayer+1) ) CYCLE
             if(R_BLK(i,j,k,iBlock) > PlanetRadius_I(iLayer) ) CYCLE
             ! to avoid eta jumps adding Eta_G
             Eta_G(i,j,k) = Eta_G(i,j,k) +Resistivity_I(iLayer+1)+(R_BLK(i,j,k,iBlock)-PlanetRadius_I(iLayer+1))* &
                  ResistivityRate(iLayer)
          end do
       end do; end do; end do
    end if

    if(DoTestMe) then
       write(*,*)'user_set_resistivity:'
       write(*,*)'Te    = ',Te_G(iTest,jTest,kTest)," [K]"
       write(*,*)'n_e   = ',nElec_G(iTest,jTest,kTest)," [m^-3]"
       write(*,*)'Eta   = ',Eta_G(iTest,jTest,kTest)*No2SI_V(UnitX_)**2/No2SI_V(UnitT_)," [m^2/s]"
    end if

  end subroutine user_set_resistivity

  !========================================================================

  subroutine user_material_properties(State_V, i,j,k,iBlock,iDir, &
       EinternalIn, TeIn, NatomicOut, AverageIonChargeOut, &
       EinternalOut, TeOut, PressureOut,   &
       CvOut, GammaOut, HeatCondOut, IonHeatCondOut, TeTiRelaxOut, &
       OpacityPlanckOut_W, OpacityRosselandOut_W, PlanckOut_W, &
       EntropyOut)

    use ModPhysics,      ONLY: No2Si_V, Si2No_V, UnitP_, UnitN_, UnitX_, ElectronPressureRatio, inv_gm1
    use ModVarIndexes,   ONLY: nVar, Rho_, p_, ExtraEInt_
    use ModConst,        ONLY: cElectronCharge, cBoltzmann, cMu, cElectronMass
    use ModAdvance,      ONLY: State_VGB
    use ModMain,         ONLY: ProcTest, iTest, jTest, kTest, BlkTest
    use ModResistivity,  ONLY: Eta0SI
    use ModProcMH,       ONLY: iProc
    use ModGeometry,     ONLY: Xyz_DGB
    use ModNumConst,     ONLY: cTiny

    !------------------------------------------------------------------------
    ! The State_V vector is in normalized units
    real, intent(in) :: State_V(nVar)
    integer, optional, intent(in) :: i, j, k, iBlock, iDir
    real, optional, intent(in)  :: EinternalIn                  ! [J/m^3]
    real, optional, intent(in)  :: TeIn                         ! [K]
    real, optional, intent(out) :: NatomicOut                   ! [1/m^3]
    real, optional, intent(out) :: AverageIonChargeOut          ! dimensionless
    real, optional, intent(out) :: EinternalOut                 ! [J/m^3]
    real, optional, intent(out) :: TeOut                        ! [K]
    real, optional, intent(out) :: PressureOut                  ! [Pa]   
    real, optional, intent(out) :: CvOut                        ! [J/(K*m^3)]  
    real, optional, intent(out) :: GammaOut                     ! dimensionless
    real, optional, intent(out) :: HeatCondOut                  ! [W/(m*K)]   
    real, optional, intent(out) :: IonHeatCondOut               ! [W/(m*K)]
    real, optional, intent(out) :: TeTiRelaxOut                 ! [1/s]  
    real, optional, intent(out) :: OpacityPlanckOut_W(nWave)    ! [1/m] 
    real, optional, intent(out) :: OpacityRosselandOut_W(nWave) ! [1/m] 
    real, optional, intent(out) :: PlanckOut_W(nWave)           ! [J/m^3] 
    real, optional, intent(out) :: EntropyOut

    real, save :: KappaCoeffSI = (cBoltzmann/cElectronCharge)**2/cMu
    real :: nElec, EtaSI, TeSI
    real, dimension(nIonFluid) :: nIon_I
    integer :: iIonFluid
    logical :: DoTest, DoTestMe=.true.

    real :: xmin, xmax, HeatCondFactor, widthmax, widthmin, xMaxyz, xMinyz

    !----------------------------------------------------------------------

    if(iBlock==BlkTest.and.i==iTest.and.j==jTest.and.k==kTest.and.iProc==ProcTest) then
       call set_oktest('user_material_properties',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    nIon_I(1:nIonFluid) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*NO2SI_V(UnitN_)
    nElec = sum(nIon_I(1:nIonFluid)*ChargeIon_I(1:nIonFluid))
    if (UseElectronPressure) then
       TeSI = State_VGB(Pe_,i,j,k,iBlock)*NO2SI_V(UnitP_)/(cBoltzmann*nElec)
    else
       TeSI = State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio/(1.+ElectronPressureRatio)*&
            NO2SI_V(UnitP_)/(cBoltzmann*nElec)
    end if
    if(present(CvOut)) CvOut = cBoltzmann*nElec*inv_gm1
    if(present(TeOut)) TeOut = TeSI
    if(present(AverageIonChargeOut)) AverageIonChargeOut = nElec/sum(nIon_I(1:nIonFluid))
    if(present(NatomicOut)) NatomicOut = sum(nIon_I(1:nIonFluid))
    if(present(TeTiRelaxOut)) TeTiRelaxOut = 0.0
    if(present(HeatCondOut)) then
       if (iBlock.ne.iNeutralBlockLast) then
          call user_neutral_atmosphere(iBlock)
       end if

       call derive_cell_diffusivity(iBlock, i, j, k, TeSI, nIon_I(1:nIonFluid), nElec, EtaSi)
       HeatCondOut = TeSI/EtaSI*KappaCoeffSI

    end if

    if(DoTestMe) then
       write(*,*)'user_material_properties:'
       write(*,*)'n_e    = ',nElec," [1/m^3]"
       write(*,*)'Te     = ',TeSI," [K]"
       if(present(CvOut)) write(*,*)'Cv     = ',CvOut,' [J/(K*m^3)]'
       if(present(HeatCondOut)) then
          write(*,*)'Eta    = ',EtaSI," [m^2/s]"
          write(*,*)'Kappa  = ',HeatCondOut," [W/(m*K)]"
          write(*,*)'Ffree  = ',nElec*sqrt(cBoltzmann*TeSi/cElectronMass)*cBoltzmann*TeSI," [W/m^2]"
       end if
       write(*,*)''
    end if
  end subroutine user_material_properties

  !========================================================================

  subroutine user_init_point_implicit

    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet
    !----------------------------------------------------------------------

    !! Source terms are evaluated explicitly!
    !RETURN

    !! All ion momenta are implicit
    if(UseElectronPressure)then
       allocate(iVarPointImpl_I(5*nIonFluid + 1))
       iVarPointImpl_I(5*nIonFluid + 1) = Pe_
    else
       allocate(iVarPointImpl_I(5*nIonFluid))
    end if

    do iFluid = 1, nIonFluid
       iVarPointImpl_I(5*iFluid-4) = iRhoIon_I(iFluid)
       iVarPointImpl_I(5*iFluid-3) = iRhoUxIon_I(iFluid)
       iVarPointImpl_I(5*iFluid-2) = iRhoUyIon_I(iFluid)
       iVarPointImpl_I(5*iFluid-1) = iRhoUzIon_I(iFluid)
       iVarPointImpl_I(5*iFluid)   = iPIon_I(iFluid)
    end do

    IsPointImplMatrixSet = .false.
    !IsAsymmetric= .false.

  end subroutine user_init_point_implicit


  !========================================================================
  subroutine user_init_session
    use ModPhysics,     ONLY: UnitRho_, UnitP_, UnitX_, Si2No_V, No2Si_V, &
         UnitRhoU_, UnitT_, ElectronPressureRatio
    use ModResistivity, ONLY: Si2NoEta
    use ModGeometry,    ONLY: TypeGeometry
    use CON_planet,     ONLY: RadiusPlanet, MassPlanet
    use ModNumConst,    ONLY: cPi

    logical :: DoTest, DoTestMe=.true.,FirstCall=.true.
    integer :: iBoundary, iIonFluid

    !----------------------------------------------------------------------
    call set_oktest('user_init_session',DoTest,DoTestMe)

!!! Introduce PlanetDensitySi etc., read those and convert 
!!! from that. Initialize these to -1. PlanetRadius should be
!!! in dimensional units.

    if (ElectronPressureRatio.le.0.) call stop_mpi('ERROR: Electron Pressure Ratio > 0 for init!')

    if (UseResistivePlanet) then
       if (TypeGeometry /= 'spherical_lnr') &
            call stop_mpi('ERROR: Correct PARAM.in, need spherical grid.')

       if(PlanetDensitySi < 0.0) &
            PlanetDensitySI  = 3.0*MassPlanet/(4.0*cPi*RadiusPlanet**3)

       if(PlanetRadius < 0.0) &
            PlanetRadius = RadiusPlanet*Si2No_V(UnitX_)

       if(PlanetPressureSi < 0.0) &
            PlanetPressureSi = 1.0e-8*No2Si_V(UnitP_)

       PlanetDensity           = PlanetDensitySi*Si2No_V(UnitRho_)
       PlanetPressure          = PlanetPressureSi*Si2No_V(UnitP_)

       PlanetRadiusSi          = PlanetRadius*No2Si_V(UnitX_)

       PlanetRadiusSi_I = PlanetRadius_I*No2Si_V(UnitX_)
       Resistivity_I = ResistivitySi_I*Si2NoEta
       do iLayer=2,nLayer
          ResistivityRate(iLayer-1) = &
               (Resistivity_I(iLayer) - Resistivity_I(iLayer-1))/&
               (PlanetRadius_I(iLayer) - PlanetRadius_I(iLayer-1))
       end do
    end if

    if(DoTestMe.and.FirstCall) then

       FirstCall = .false.
       if(UseResistivePlanet) then
          write(*,*) 'Resistiv Planet Model'       
          write(*,*) 'Planet density     = ',PlanetDensitySi
          write(*,*) 'Planet pressure    = ',PlanetPressureSi
          write(*,*) 'Planet radius      = ',PlanetRadiusSi
          if(nLayer > 0 ) then
             write(*,*) ''
             write(*,*) '   |-------- Planet Resistivity Profile -----|'
             write(*,*) '       Radius [m]          Magnetic diffusivity [m^2/s]'
             do iLayer =1,nLayer 
                write(*,"(A7,E10.3,A15,E10.3)") " ",PlanetRadiusSi_I(iLayer)," ",&
                     ResistivitySi_I(iLayer)
             end do
          else
             write(*,*) 'Conducting Planet (eta =0)'
          end if
       end if
       write(*,*) ''
       write(*,*)'Units'
       write(*,*)'No2SI_V(UnitRho_)   =',No2SI_V(UnitRho_)
       write(*,*)'No2SI_V(UnitRhoU_)  =',No2SI_V(UnitRhoU_)
       write(*,*)'No2SI_V(UnitP_)     =',No2SI_V(UnitP_)
       write(*,*)'No2SI_V(UnitT_)     =',No2SI_V(UnitT_)

    end if

  end subroutine user_init_session


  !========================================================================
  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional,&
       PlotVar_G, PlotVarBody, UsePlotVarBody,&
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModAdvance,    ONLY: State_VGB, RhoUx_, RhoUy_, RhoUz_
    use ModPhysics,    ONLY: No2Si_V, Si2No_V, UnitP_, UnitN_, UnitU_, UnitT_, &
         ElectronCharge, ElectronPressureRatio, UnitX_
    use ModVarIndexes, ONLY: Rho_, P_, Pe_
    use ModConst,      ONLY: cBoltzmann
    use ModCurrent,    ONLY: get_current
    use ModMultiFluid, ONLY: MassIon_I
    use ModMain,       ONLY: Dt_BLK
    use BATL_lib,      ONLY: CellVolume_GB

    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(-1:nI+2, -1:nJ+2, -1:nK+2)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    integer :: i, j, k, iIonFluid
    real :: nElec
    real, dimension(3)           :: Current_I, uIonMean_I
    real, dimension(nIonFluid)   :: nIon_I
    real, dimension(3,nIonFluid) :: uIon_I


    !--------------------------------------------------------------------------

    IsFound = .true.

    if (iBlock.ne.iNeutralBlockLast) then
       call user_neutral_atmosphere(iBlock)
    end if

    select case(NameVar)
    case('nn1')
       NameIdlUnit = '1/cm^3'
       NameTecUnit = '[1/cm^3]'
       !do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = NnNeutral_IG(1,i-MinI+1,j-MinJ+1,k-MinK+1)/1E6
       end do; end do; end do

    case('unx1')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[cm/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = UnxNeutral_IG(1,i-MinI+1,j-MinJ+1,k-MinK+1)/1E3 !! x direction
       end do; end do; end do

    case('uny1')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[cm/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = UnyNeutral_IG(1,i-MinI+1,j-MinJ+1,k-MinK+1)/1E3 !! y direction
       end do; end do; end do

    case('unz1')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[km/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = UnzNeutral_IG(1,i-MinI+1,j-MinJ+1,k-MinK+1)/1E3 !! z direction
       end do; end do; end do

    case('tn1')
       NameIdlUnit = 'K'
       NameTecUnit = '[K]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = TnNeutral_IG(1,i-MinI+1,j-MinJ+1,k-MinK+1)
       end do; end do; end do

    case('te')
       NameIdlUnit = 'K'
       NameTecUnit = '[K]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          nIon_I(1:nIonFluid) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
          nElec = sum(nIon_I(1:nIonFluid)*ChargeIon_I(1:nIonFluid))
          if(UseElectronPressure)then
             PlotVar_G(i,j,k) = State_VGB(Pe_,i,j,k,iBlock)*No2SI_V(UnitP_)/&
                  (cBoltzmann*nElec)
          else
             PlotVar_G(i,j,k) = State_VGB(P_,i,j,k,iBlock)*No2SI_V(UnitP_)*ElectronPressureRatio/&
                  (1.+ElectronPressureRatio)/(cBoltzmann*nElec)
          end if
       end do; end do; end do

    case('ti1')
       NameIdlUnit = 'K'
       NameTecUnit = '[K]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = State_VGB(iPIon_I(1),i,j,k,iBlock)*NO2SI_V(UnitP_)/ &
               (cBoltzmann*State_VGB(iRhoIon_I(1),i,j,k,iBlock)/MassIon_I(1)*NO2SI_V(UnitN_))
       end do; end do; end do

    case('uex')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[km/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          call get_current(i,j,k,iBlock,Current_I)
          nIon_I(1:nIonFluid) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
          uIon_I(1,1:nIonFluid) = State_VGB(iRhoUxIon_I,i,j,k,iBlock) / &
               State_VGB(iRhoIon_I,i,j,k,iBlock)*No2SI_V(UnitU_)
          nElec = sum(nIon_I(1:nIonFluid)*ChargeIon_I(1:nIonFluid))
          uIonMean_I(1) = 0.
          do iIonFluid=1,nIonFluid
             uIonMean_I(1) = uIonMean_I(1)+nIon_I(iIonFluid)* &
                  uIon_I(1,iIonFluid)/nElec*ChargeIon_I(iIonFluid)
          end do
          PlotVar_G(i,j,k) = (uIonMean_I(1)-Current_I(1)/(nElec*Si2No_V(UnitN_)*&
               ElectronCharge)*No2SI_V(UnitU_))/1E3
       end do; end do; end do

    case('uey')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[km/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          call get_current(i,j,k,iBlock,Current_I)
          nIon_I(1:nIonFluid) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
          uIon_I(2,1:nIonFluid) = State_VGB(iRhoUyIon_I,i,j,k,iBlock) / &
               State_VGB(iRhoIon_I,i,j,k,iBlock)*No2SI_V(UnitU_)
          nElec = sum(nIon_I(1:nIonFluid)*ChargeIon_I(1:nIonFluid))
          uIonMean_I(2) = 0.
          do iIonFluid=1,nIonFluid
             uIonMean_I(2) = uIonMean_I(2)+nIon_I(iIonFluid)* &
                  uIon_I(2,iIonFluid)/nElec*ChargeIon_I(iIonFluid)
          end do
          PlotVar_G(i,j,k) = (uIonMean_I(2)-Current_I(2)/(nElec*Si2No_V(UnitN_)*&
               ElectronCharge)*No2SI_V(UnitU_))/1E3
       end do; end do; end do

    case('uez')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[km/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          call get_current(i,j,k,iBlock,Current_I)
          nIon_I(1:nIonFluid) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
          uIon_I(3,1:nIonFluid) = State_VGB(iRhoUzIon_I,i,j,k,iBlock) / &
               State_VGB(iRhoIon_I,i,j,k,iBlock)*No2SI_V(UnitU_)
          nElec = sum(nIon_I(1:nIonFluid)*ChargeIon_I(1:nIonFluid))
          uIonMean_I(3) = 0.
          do iIonFluid=1,nIonFluid
             uIonMean_I(3) = uIonMean_I(3)+nIon_I(iIonFluid)* &
                  uIon_I(3,iIonFluid)/nElec*ChargeIon_I(iIonFluid)
          end do
          PlotVar_G(i,j,k) = (uIonMean_I(3)-Current_I(3)/(nElec*Si2No_V(UnitN_)*&
               ElectronCharge)*No2SI_V(UnitU_))/1E3
       end do; end do; end do

    case('dt')
       NameIdlUnit = 's'
       NameTecUnit = '[s]'
       PlotVar_G(:,:,:) = Dt_BLK(iBlock)*No2SI_V(UnitT_)

    case('mlop')
       NameIdlUnit = 'kg*m^-3s^-1'   
       NameTecUnit = '[kg*m^-3s^-1]'
       do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
          !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = MassLoading(1,i,j,k,iBlock)
       end do; end do; end do

    case('ionop')
       ! total ionization rate of O2 + hv or e -> Op + O + (2)e 
       NameIdlUnit = '1/s'   
       NameTecUnit = '[1/s]'
       do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
          !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = MassLoading(2,i,j,k,iBlock)
       end do; end do; end do

    case('eionop')
       ! electron impact ionization rate O2 + e -> Op + O + 2e 
       NameIdlUnit = ' '   
       NameTecUnit = '[ ]'
       do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
          !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = MassLoading(3,i,j,k,iBlock)
       end do; end do; end do

    case('recop')
       NameIdlUnit = 'kg*m^-3s^-1'   
       NameTecUnit = '[kg*m^-3s^-1]'
       do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
          !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = MassLoading(4,i,j,k,iBlock)
       end do; end do; end do

    case('mlo2p')
       NameIdlUnit = 'kg*m^-3s^-1'   
       NameTecUnit = '[kg*m^-3s^-1]'
       do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
          !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = MassLoading(5,i,j,k,iBlock)
       end do; end do; end do
    case('iono2p')
       ! total ionization rate of O2 + hv or e -> O2p + (2)e 
       NameIdlUnit = '1/s'   
       NameTecUnit = '[1/s]'
       do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
          !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = MassLoading(6,i,j,k,iBlock)
       end do; end do; end do
    case('eiono2p')
       ! electron impact ionization rate O2 + e -> Op + 2e 
       NameIdlUnit = '1/s'   
       NameTecUnit = '[1/s]'
       do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
          !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = MassLoading(7,i,j,k,iBlock)
       end do; end do; end do
    case('reco2p')
       NameIdlUnit = 'kg*m^-3s^-1'   
       NameTecUnit = '[kg*m^-3s^-1]'
       do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
          !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = MassLoading(8,i,j,k,iBlock)
       end do; end do; end do
    case('cellvol')
       NameIdlUnit = 'm^3'   
       NameTecUnit = '[m^3]'
       do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
          !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = CellVolume_GB(i,j,k,iBlock)*No2SI_V(UnitX_)**3
       end do; end do; end do
       ! case('fluxlim')
       !    NameIdlUnit = ' '   
       !    NameTecUnit = '[ ]'
       !    do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
       !       PlotVar_G(i,j,k) = FluxLimited_GB(i,j,k,iBlock)
       !    end do; end do; end do
    case default
       IsFound = .false.
    end select

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

  end subroutine user_set_plot_var



  !========================================================================

  subroutine user_preset_conditions(i,j,k,iBlock)
    ! This is applied as initial conditions and in the upstream boundary for the semi-implicit heat conduction
    use ModAdvance,  ONLY: P_, Pe_, State_VGB
    use ModPhysics,  ONLY: SW_N, SW_Ux, SW_Uy, SW_Uz, SW_T_Dim, LowDensityRatio, ElectronPressureRatio, Io2No_V, &
         UnitTemperature_

    integer,intent(in) :: i, j, k, iBlock

    State_VGB(OpRho_,i,j,k,iBlock)     = SW_n*MassIon_I(Op_)
    State_VGB(OpRhoUx_,i,j,k,iBlock)   = SW_n*MassIon_I(Op_)*SW_Ux
    State_VGB(OpRhoUy_,i,j,k,iBlock)   = SW_n*MassIon_I(Op_)*SW_Uy
    State_VGB(OpRhoUz_,i,j,k,iBlock)   = SW_n*MassIon_I(Op_)*SW_Uz
    State_VGB(OpP_,i,j,k,iBlock)       = SW_n*SW_T_dim*Io2No_V(UnitTemperature_)

    State_VGB(O2pRho_,i,j,k,iBlock)    = SW_n*LowDensityRatio*MassIon_I(O2p_)
    State_VGB(O2pRhoUx_,i,j,k,iBlock)  = SW_n*LowDensityRatio*MassIon_I(O2p_)*SW_Ux
    State_VGB(O2pRhoUy_,i,j,k,iBlock)  = SW_n*LowDensityRatio*MassIon_I(O2p_)*SW_Uy
    State_VGB(O2pRhoUz_,i,j,k,iBlock)  = SW_n*LowDensityRatio*MassIon_I(O2p_)*SW_Uz
    State_VGB(O2pP_,i,j,k,iBlock)      = SW_n*LowDensityRatio*SW_T_dim*Io2No_V(UnitTemperature_)

    State_VGB(Rho_,i,j,k,iBlock)       = sum(State_VGB(iRhoIon_I,i,j,k,iBlock))
    State_VGB(RhoUx_,i,j,k,iBlock)     = sum(State_VGB(iRhoUxIon_I,i,j,k,iBlock))
    State_VGB(RhoUy_,i,j,k,iBlock)     = sum(State_VGB(iRhoUyIon_I,i,j,k,iBlock))
    State_VGB(RhoUz_,i,j,k,iBlock)     = sum(State_VGB(iRhoUzIon_I,i,j,k,iBlock))

    if(UseElectronPressure) then
       State_VGB(P_,i,j,k,iBlock)      = sum(State_VGB(iPIon_I,i,j,k,iBlock))
       State_VGB(Pe_,i,j,k,iBlock)     = State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio
    else
       State_VGB(P_,i,j,k,iBlock)      = sum(State_VGB(iPIon_I,i,j,k,iBlock))* &
            (1.+ElectronPressureRatio)
    end if

  end subroutine user_preset_conditions

  !========================================================================

  subroutine user_set_ICs(iBlock)
    use ModIO,       ONLY: restart
    use ModProcMH,   ONLY: iProc
    use ModMain,     ONLY: iTest, jTest, kTest, ProcTest, BlkTest
    use ModAdvance,  ONLY: P_, Pe_, State_VGB
    use ModPhysics,  ONLY: ElectronPressureRatio, No2Si_V, UnitRho_, UnitRhoU_, UnitP_, rBody, SW_Bz
    use ModConst,    ONLY: cBoltzmann
    use ModGeometry, ONLY: R_BLK

    integer, intent(in) :: iBlock

    logical :: DoTest, DoTestMe=.true.
    integer :: i, j, k, iIonFluid
    ! !-------------------------------------------------------------------------
    if(iProc==PROCtest .and. iBlock==BLKtest)then
       call set_oktest('user_set_ICs', DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    !if (iBlock.ne.iNeutralBlockLast) then
    !   call user_neutral_atmosphere(iBlock)
    !end if

    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI

       call user_preset_conditions(i,j,k,iBlock)

       !! inside the planet/moon
       if(.not.UseResistivePlanet) CYCLE
       if(R_BLK(i+2,j,k,iBlock) > PlanetRadius) CYCLE
       State_VGB(iRhoIon_I,i,j,k,iBlock) = PlanetDensity/nIonFluid
       State_VGB(iRhoUxIon_I,i,j,k,iBlock) = 0.0
       State_VGB(iRhoUyIon_I,i,j,k,iBlock) = 0.0
       State_VGB(iRhoUzIon_I,i,j,k,iBlock) = 0.0
       State_VGB(iPIon_I,i,j,k,iBlock) = PlanetPressure/nIonFluid
       State_VGB(Rho_,i,j,k,iBlock) = PlanetDensity
       State_VGB(P_,i,j,k,iBlock) = PlanetPressure
       State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock) = 0.0
       ! State_VGB(Bx_:By_,i,j,k,iBlock) = 0.0
       ! State_VGB(Bz_,i,j,k,iBlock) = SW_Bz                            !! Bz component is assumed to be fully diffused through

    end do; end do ; end do


    if(DoTestMe) then
       i=iTest ; j=jTest ; k=kTest
123    format (A13,ES25.16,A15)
       do iIonFluid=1,nIonFluid
          write(*,*)'Ion species #',iIonFluid,': ',NameFluid_I(iIonFluid+1)
          write(*,123)'Rho       = ',State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRho_)," [kg/m^3]"
          write(*,123)'rhoUx     = ',State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'rhoUy     = ',State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'rhoUz     = ',State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'Pi        = ',State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitP_)," [Pa]"
       end do
       write(*,*)''
       write(*,*)'Total:'
       write(*,123)'Rho       = ',State_VGB(Rho_,i,j,k,iBlock)*No2SI_V(UnitRho_)," [kg/m^3]"
       write(*,123)'rhoUx     = ',State_VGB(RhoUx_,i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
            " [kg/(m^2*s)]"
       write(*,123)'rhoUy     = ',State_VGB(RhoUy_,i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
            " [kg/(m^2*s)]"
       write(*,123)'rhoUz     = ',State_VGB(RhoUz_,i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
            " [kg/(m^2*s)]"
       if (UseElectronPressure) then
          write(*,123)'PiTot     = ',State_VGB(P_,i,j,k,iBlock)*No2SI_V(UnitP_)," [Pa]"
          write(*,123)'Pe        = ',State_VGB(Pe_,i,j,k,iBlock)*No2SI_V(UnitP_)," [Pa]"
       else
          write(*,123)'PiTot     = ',State_VGB(P_,i,j,k,iBlock)/(1.+ElectronPressureRatio)* &
               No2SI_V(UnitP_)," [Pa]"
          write(*,123)'Pe        = ',State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio/&
               (1.+ElectronPressureRatio)*No2SI_V(UnitP_)," [Pa]"
       end if

    end if

  end subroutine user_set_ICs

  !========================================================================

  subroutine user_set_cell_boundary(iBlock,iSide, TypeBc, IsFound)

    use ModAdvance,  ONLY: State_VGB
    use ModImplicit, ONLY: StateSemi_VGB, iTeImpl
    use ModSize,     ONLY: nI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModPhysics,  ONLY: Si2No_V, UnitTemperature_, SW_n, SW_Ux, SW_Uy, SW_Uz, &
         LowDensityRatio, ElectronPressureRatio, SW_T_dim, Io2No_V
    use ModGeometry, ONLY: Xyz_DGB, ExtraBc_, IsBoundaryCell_GI

    integer,          intent(in)  :: iBlock, iSide
    character(len=*), intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer:: i, j, k
    real :: TeSi

    character(len=*), parameter :: NameSub = 'user_set_cell_boundary'
    !-------------------------------------------------------------------
    IsFound = .true.

    if(TypeBc == 'usersemi')then

       ! Overwrite outermost ghost cells in spherical mesh
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = nI+1, MaxI
          if(.not.IsBoundaryCell_GI(i,j,k,ExtraBc_)) CYCLE

          !call user_preset_conditions(i,j,k,iBlock)
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, TeOut=TeSi)
          StateSemi_VGB(iTeImpl,i,j,k,iBlock) = TeSi*Si2No_V(UnitTemperature_)

       end do; end do; end do

       RETURN
    elseif(TypeBc == 'usersemilinear')then
       RETURN
    end if


  end subroutine user_set_cell_boundary

  !============================================================================

  subroutine user_get_log_var(VarValue, TypeVar, Radius)

    use ModMain,       ONLY: Dt_BLK, BLKtest, ProcTest, Unused_B
    use ModPhysics,    ONLY: No2SI_V, UnitT_, UnitX_
    use ModProcMH,     ONLY: iProc
    use ModGeometry,   ONLY: R_BLK
    use BATL_lib,      ONLY: CellVolume_GB

    real, intent(out)             :: VarValue
    character (len=*), intent(in) :: TypeVar
    real, intent(in), optional    :: Radius

    integer :: i, j, k, iBlock
    real :: rMax = 4.0 ! Maximum distance to evaluate mass loading

    character (len=*), parameter  :: NameSub = 'user_get_log_var'
    !-------------------------------------------------------------------------
    select case(TypeVar)
    case('dtpnt')
       ! local time step
       VarValue = 0.
       if (iProc == ProcTest) VarValue = Dt_BLK(BLKtest)*No2SI_V(UnitT_)

    case('mlop')
       ! Mass loading of Op ions [kg/s] 

       VarValue = 0.0
       do iBlock=1,nBLK
          if(Unused_B(iBlock))CYCLE
          do k=1,nK; do j=1,nJ; do i=1,nI

             ! ommit cells outside 4 rE (no massloading that far)
             if (R_BLK(i,j,k,iBlock) > rMax) CYCLE
             ! ommit cells inside body
             if (R_BLK(i,j,k,iBlock) < 1.0) CYCLE

             ! Sum up local mass loading
             VarValue = VarValue+ &
                  MassLoading(1,i,j,k,iBlock)*CellVolume_GB(i,j,k,iBlock)

             ! Test (should return volume 4*pi/3*(4^3-1^3)=263.89)
             ! VarValue = VarValue+ &
             !      CellVolume_GB(i,j,k,iBlock)
             
          end do; end do; end do
       end do
       VarValue = VarValue*No2SI_V(UnitX_)**3
       !VarValue = VarValue !Test

    case('recop')
       ! Recombination of Op ions [kg/s] 

       VarValue = 0.0
       do iBlock=1,nBLK
          if(Unused_B(iBlock))CYCLE
          do k=1,nK; do j=1,nJ; do i=1,nI

             ! ommit cells outside 4 rE (no massloading that far)
             if (R_BLK(i,j,k,iBlock) > rMax) CYCLE
             ! ommit cells inside body
             if (R_BLK(i,j,k,iBlock) < 1.0) CYCLE

             ! Sum up local mass loading
             VarValue = VarValue+ &
                  MassLoading(4,i,j,k,iBlock)*CellVolume_GB(i,j,k,iBlock)
             
          end do; end do; end do
       end do
       VarValue = VarValue*No2SI_V(UnitX_)**3

    case('reco2p')
       ! Recombination of O2p ions [kg/s]

       VarValue = 0.0
       do iBlock=1,nBLK
          if(Unused_B(iBlock))CYCLE
          do k=1,nK; do j=1,nJ; do i=1,nI

             ! ommit cells outside 4 rE (no massloading that far)
             if (R_BLK(i,j,k,iBlock) > rMax) CYCLE
             ! ommit cells inside body
             if (R_BLK(i,j,k,iBlock) < 1.0) CYCLE

             ! Sum up local mass loading
             VarValue = VarValue+ &
                  MassLoading(8,i,j,k,iBlock)*CellVolume_GB(i,j,k,iBlock)
             
          end do; end do; end do
       end do
       VarValue = VarValue*No2SI_V(UnitX_)**3

    case('mlo2p')
       ! Mass loading of O2p ions [kg/s] 

       VarValue = 0.0
       do iBlock=1,nBLK
          if(Unused_B(iBlock))CYCLE
          do k=1,nK; do j=1,nJ; do i=1,nI

             ! ommit cells outside 4 rE (no massloading that far)
             if (R_BLK(i,j,k,iBlock) > rMax) CYCLE
             ! ommit cells inside body
             if (R_BLK(i,j,k,iBlock) < 1.0) CYCLE

             ! Sum up local mass loading
             VarValue = VarValue+ &
                  MassLoading(5,i,j,k,iBlock)*CellVolume_GB(i,j,k,iBlock)
             
          end do; end do; end do
       end do
       VarValue = VarValue*No2SI_V(UnitX_)**3

    case default
       VarValue = -7777.0
    end select

  end subroutine user_get_log_var

  !============================================================================

  subroutine user_set_face_boundary(VarsGhostFace_V)

    use ModSize,         ONLY: x_
    use ModVarIndexes,   ONLY: nVar, Bx_, Bz_
    use ModFaceBoundary, ONLY: TimeBc, iFace, jFace, kFace, FaceCoords_D, &
         iBoundary, VarsTrueFace_V
    use ModSolarwind,    ONLY: get_solar_wind_point
    !    use ModB0,           ONLY: B0_DX
    use ModMain,         ONLY: body1_
    use ModPhysics,      ONLY: LowDensityRatio, SW_Ux, SW_Uy, SW_Uz, SW_n, SW_T_dim, &
         ElectronPressureRatio, UnitTemperature_, Io2No_V, SW_Bx, SW_By, SW_Bz, &
         NO2SI_V, UnitP_, UnitRho_, UnitRhoU_, UnitU_, UnitB_, UnitN_, Io2SI_V, &
         BodyRho_I, BodyP_I, BodyNDim_I, BodyTDim_I

    logical :: FirstCall = .true., DoTest, DoTestMe=.true.
    integer :: iIonFluid
    real    :: UdotR(nIonFluid), URefl_D(1:3,nIonFluid), rFace, cos_theta
    !real    :: BdotR, BRefl_D(1:3)

    real, intent(out):: VarsGhostFace_V(nVar)

    !------------------------------------------------------------------------

    if(DoTestMe.and.FirstCall) then
       call set_oktest('user_set_face_boundary',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    ! write(*,*)'iBoundary = ',iBoundary

    !! Outer boundaries
    if(iBoundary >= 0) then

       call get_solar_wind_point(TimeBc, FaceCoords_D(x_),VarsGhostFace_V)

       ! Magnetospheric oxygen ions
       VarsGhostFace_V(OpRho_)     = SW_n*MassIon_I(Op_)*Io2No_V(UnitRho_)
       VarsGhostFace_V(OpRhoUx_)   = SW_Ux
       VarsGhostFace_V(OpRhoUy_)   = SW_Uy
       VarsGhostFace_V(OpRhoUz_)   = SW_Uz
       VarsGhostFace_V(OpP_)       = SW_n*Io2No_V(UnitN_)*SW_T_dim*Io2No_V(UnitTemperature_) ! cBoltzmann is in Io2No_V(UnitTemperature_)

       ! O2p pick-up ions
       VarsGhostFace_V(O2pRho_)   = SW_n*MassIon_I(O2p_)*LowDensityRatio*Io2No_V(UnitRho_)
       VarsGhostFace_V(O2pRhoUx_) = SW_Ux
       VarsGhostFace_V(O2pRhoUy_) = SW_Uy
       VarsGhostFace_V(O2pRhoUz_) = SW_Uz
       VarsGhostFace_V(O2pP_)     = SW_n*Io2No_V(UnitN_)*LowDensityRatio*SW_T_dim*Io2No_V(UnitTemperature_) ! cBoltzmann is in Io2No_V(UnitTemperature_)

       VarsGhostFace_V(Rho_)       = sum(VarsGhostFace_V(iRhoIon_I))
       VarsGhostFace_V(RhoUx_)     = SW_Ux
       VarsGhostFace_V(RhoUy_)     = SW_Uy
       VarsGhostFace_V(RhoUz_)     = SW_Uz

       ! VarsGhostFace_V(Bx_:Bz_)    = VarsGhostFace_V(Bx_:Bz_) - B0_DX(:,iFace, jFace, kFace)
       ! VarsGhostFace_V(Bx_)        = SW_Bx - B0_DX(1,iFace, jFace, kFace)
       ! VarsGhostFace_V(By_)        = SW_By - B0_DX(2,iFace, jFace, kFace)
       ! VarsGhostFace_V(Bz_)        = SW_Bz - B0_DX(3,iFace, jFace, kFace)

       if(UseElectronPressure) then
          VarsGhostFace_V(P_)      = sum(VarsGhostFace_V(iPIon_I))
          VarsGhostFace_V(Pe_)     = VarsGhostFace_V(P_)*ElectronPressureRatio
       else
          VarsGhostFace_V(P_)      = sum(VarsGhostFace_V(iPIon_I))*(1.+ElectronPressureRatio)
       end if

    else if ((iBoundary <= body1_).and.(.not. UseResistivePlanet)) then
    !! Body boundaries

       !! Projection length of U_ and B_ on the local surface radius vector
       !! in units of the surface radius vector [rBody]

       do iIonFluid=1,nIonFluid
          UdotR(iIonFluid) = dot_product(VarsTrueFace_V(iRhoUx_I(iIonFluid):iRhoUz_I(iIonFluid)), &
               FaceCoords_D)/dot_product(FaceCoords_D,FaceCoords_D)
          !! Projection vectors

          URefl_D(1:3,iIonFluid) = UdotR(iIonFluid)*FaceCoords_D(1:3)
       end do
       !BdotR = dot_product(VarsTrueFace_V(Bx_:Bz_),FaceCoords_D)/ & 
       !     dot_product(FaceCoords_D,FaceCoords_D) 

       !! Projection vectors
       !BRefl_D = BdotR*FaceCoords_D

       !! Floating boundary conditions allowing inflow 
       VarsGhostFace_V = VarsTrueFace_V

       !! Bz component propagated through moon, Bx and By didn't  
       !  VarsGhostFace_V(Bx_:By_) = 0.0 
       !  VarsGhostFace_V(Bz_)     = SW_Bz 

       !! set outward flux body value (Europa's surface not considered as plasma source)
       !! leave inward flux untouched       
       do iIonFluid=1,nIonFluid
          if (UdotR(iIonFluid) > 0.0) then
             VarsGhostFace_V(iUx_I(iIonFluid):iUz_I(iIonFluid)) = 0.0
             VarsGhostFace_V(iRhoUxIon_I(iIonFluid):iRhoUzIon_I(iIonFluid)) = 0.0
             VarsGhostFace_V(iRhoIon_I(iIonFluid)) = BodyRho_I(iIonFluid)
             VarsGhostFace_V(iPIon_I(iIonFluid)) = BodyP_I(iIonFluid)
          endif
       end do

       ! write(*,*)'UseFixedBoundary = ', UseFixedBoundary ! public vs nonpublic ???

       if (UseFixedBoundary) then
          ! rFace = sqrt(dot_product(FaceCoords_D,FaceCoords_D))
          ! cos_theta=sum(Dir_I(1:3)*FaceCoords_D(1:3))/(rFace*Dir_length)
          ! if (cos_theta >= 0.) then

          ! Use density and temperature from #UPSTREAMBOUNDARY in PARAM.in file
          if (sum(Dir_I(1:3)*FaceCoords_D(1:3)) >= 0.) then
             do iIonFluid=1,nIonFluid
                VarsGhostFace_V(iRhoIon_I(iIonFluid)) = BodyNDimUp_I(iIonFluid)*Io2No_V(UnitN_)* &
                     MassIon_I(iIonFluid)
                VarsGhostFace_V(iPIon_I(iIonFluid)) = BodyNDimUp_I(iIonFluid)*Io2No_V(UnitN_)* &
                     BodyTDimUp_I(iIonFluid)*Io2No_V(UnitTemperature_)
                ! write(*,*)'Ion species #',iIonFluid,': ',NameFluid_I(iIonFluid+1)
                ! write(*,*)'VarsGhostFaceRho   = ',VarsGhostFace_V(iRhoIon_I(iIonFluid))*No2SI_V(UnitRho_)," [kg/m^3]"
                ! write(*,*)'VarsGhostFaceP     = ',VarsGhostFace_V(iPIon_I(iIonFluid))*No2SI_V(UnitP_)," [Pa]"
                ! write(*,*)'mi                 = ',MassIon_I(iIonFluid)
                ! write(*,*)'BodyNDimUp         = ',BodyNDimUp_I(iIonFluid)
                ! write(*,*)'BodyTDimUp         = ',BodyTDimUp_I(iIonFluid)
             end do
          else
             do iIonFluid=1,nIonFluid
                VarsGhostFace_V(iRhoIon_I(iIonFluid)) = BodyNDimDo_I(iIonFluid)*Io2No_V(UnitN_)* &
                     MassIon_I(iIonFluid)
                VarsGhostFace_V(iPIon_I(iIonFluid)) = BodyNDimDo_I(iIonFluid)*Io2No_V(UnitN_)* &
                     BodyTDimDo_I(iIonFluid)*Io2No_V(UnitTemperature_)
                ! write(*,*)'Ion species #',iIonFluid,': ',NameFluid_I(iIonFluid+1)
                ! write(*,*)'VarsGhostFaceRho   = ',VarsGhostFace_V(iRhoIon_I(iIonFluid))*No2SI_V(UnitRho_)," [kg/m^3]"
                ! write(*,*)'VarsGhostFaceP     = ',VarsGhostFace_V(iPIon_I(iIonFluid))*No2SI_V(UnitP_)," [Pa]"
                ! write(*,*)'mi                 = ',MassIon_I(iIonFluid)
                ! write(*,*)'BodyNDimDo         = ',BodyNDimDo_I(iIonFluid)
                ! write(*,*)'BodyTDimDo         = ',BodyTDimDo_I(iIonFluid)
             end do
          end if
       end if
       VarsGhostFace_V(Rho_)   = sum(VarsGhostFace_V(iRhoIon_I))
       VarsGhostFace_V(RhoUx_) = sum(VarsGhostFace_V(iRhoIon_I)*VarsGhostFace_V(iRhoUxIon_I))/ &
            sum(VarsGhostFace_V(iRhoIon_I))
       VarsGhostFace_V(RhoUy_) = sum(VarsGhostFace_V(iRhoIon_I)*VarsGhostFace_V(iRhoUyIon_I))/ &
            sum(VarsGhostFace_V(iRhoIon_I))
       VarsGhostFace_V(RhoUz_) = sum(VarsGhostFace_V(iRhoIon_I)*VarsGhostFace_V(iRhoUzIon_I))/ &
            sum(VarsGhostFace_V(iRhoIon_I))
       VarsGhostFace_V(P_)     = sum(VarsGhostFace_V(iPIon_I))

    end if

    if(DoTestMe) then

       !FirstCall = .false.
       write(*,*)'Boundary No =',iBoundary
       write(*,*)'VarsGhostFaces'
       do iIonFluid=1,nIonFluid
          write(*,*)'Ion species #',iIonFluid,': ',NameFluid_I(iIonFluid+1)
          write(*,*)'FaceCoordsX        = ',FaceCoords_D(1)," [RE]"
          write(*,*)'FaceCoordsY        = ',FaceCoords_D(2)," [RE]"
          write(*,*)'FaceCoordsZ        = ',FaceCoords_D(3)," [RE]"
          write(*,*)'VarsGhostFaceRho   = ',VarsGhostFace_V(iRhoIon_I(iIonFluid))*No2SI_V(UnitRho_)," [kg/m^3]"
          write(*,*)'VarsGhostFaceUx    = ',VarsGhostFace_V(iRhoUxIon_I(iIonFluid))*No2SI_V(UnitU_)," [m/s]"
          write(*,*)'VarsGhostFaceUy    = ',VarsGhostFace_V(iRhoUyIon_I(iIonFluid))*No2SI_V(UnitU_)," [m/s]"
          write(*,*)'VarsGhostFaceUz    = ',VarsGhostFace_V(iRhoUzIon_I(iIonFluid))*No2SI_V(UnitU_)," [m/s]"
          write(*,*)'VarsGhostFaceP     = ',VarsGhostFace_V(iPIon_I(iIonFluid))*No2SI_V(UnitP_)," [Pa]"
       end do
       write(*,*)''
       write(*,*)'Total ion fluid:'
       write(*,*)'VarsGhostFaceRho   = ',sum(VarsGhostFace_V(iRhoIon_I)*No2SI_V(UnitRho_))," [kg/m^3]"
       write(*,*)'VarsGhostFaceUx    = ',sum(VarsGhostFace_V(iRhoIon_I)*VarsGhostFace_V(iRhoUxIon_I))/ &
            sum(VarsGhostFace_V(iRhoIon_I))*No2SI_V(UnitU_)," [m/s]"
       write(*,*)'VarsGhostFaceUx    = ',sum(VarsGhostFace_V(iRhoIon_I)*VarsGhostFace_V(iRhoUyIon_I))/ &
            sum(VarsGhostFace_V(iRhoIon_I))*No2SI_V(UnitU_)," [m/s]"
       write(*,*)'VarsGhostFaceUx    = ',sum(VarsGhostFace_V(iRhoIon_I)*VarsGhostFace_V(iRhoUzIon_I))/ &
            sum(VarsGhostFace_V(iRhoIon_I))*No2SI_V(UnitU_)," [m/s]"
       write(*,*)'VarsGhostFaceP   = ',VarsGhostFace_V(P_)*No2SI_V(UnitP_)," [Pa]"
       if (UseElectronPressure) then
          write(*,*) ''
          write(*,*)'VarsGhostFacePe  = ',VarsGhostFace_V(Pe_)*No2SI_V(UnitP_)," [Pa]"
       end if
       write(*,*)''
       write(*,*)'Magnetic field:'
       write(*,*)'VarsGhostFaceBx    = ',VarsGhostFace_V(Bx_)*No2SI_V(UnitB_)," [T]"
       write(*,*)'VarsGhostFaceBy    = ',VarsGhostFace_V(By_)*No2SI_V(UnitB_)," [T]"
       write(*,*)'VarsGhostFaceBz    = ',VarsGhostFace_V(Bz_)*No2SI_V(UnitB_)," [T]"
       write(*,*)''

    end if

  end subroutine user_set_face_boundary

  !============================================================================

  subroutine user_set_boundary_cells(iBlock)

    use ModGeometry, ONLY: ExtraBc_, IsBoundaryCell_GI, Xyz_DGB, x1, x2
    use ModPhysics,  ONLY: SW_Ux    

    implicit none

    integer, intent(in):: iBlock

    character (len=*), parameter :: Name='user_set_boundary_cells'

    !--------------------------------------------------------------------------
    ! For inflow in positive x direction
    if (SW_Ux > 0. ) then
       IsBoundaryCell_GI(:,:,:,ExtraBc_) = Xyz_DGB(x_,:,:,:,iBlock) < x1
    else
       ! For inflow in negative x direction
       IsBoundaryCell_GI(:,:,:,ExtraBc_) = Xyz_DGB(x_,:,:,:,iBlock) > x2
    end if

  end subroutine user_set_boundary_cells

end module ModUser
