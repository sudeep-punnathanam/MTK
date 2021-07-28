module rosenbluth
  use, intrinsic :: ieee_exceptions
  use, intrinsic :: ieee_features, only: ieee_underflow_flag
  use consts, only: PR, PI, TWOPI, strlen, lstrlen, MAX_NO_OF_SPECIES, Dashed_Line
  use variables, only: SetReducedUnits, CurrentCoordinates, CurrentSimulationCell, TrialCoordinates, &
    MainCellList, NumberOfSimulations, DisplayFrequency, NumberOfProductionCycles, beta, Temperature, &
    K_B, MAX_NO_OF_SIMULATIONS, REDUCED_UNITS, NumberOfSystems
  use random, only: RandomNumber, UniformRandomRotationMatrix, RandomUnitVector, SetRandomNumberSeed
  use storage, only: StorageInteractions, UpdateStorageInteractions, ResetStorageInteractions, &
    Operator(+), InteractionAccumulator
  use config, only: AllocateMemoryForCoordinates, ReAllocateMemoryForCoordinates
  use atoms_and_molecules, only: ReadAtomInfo, ReadMoleculeInfo, Molecule, GetCenterOfMass, NumberOfSpecies
  use inter, only: MoleculeSystemShortRangePairwiseInteraction, AtomSystemShortRangePairwiseInteraction, &
    AllocateMemoryToSeparationVectors
  use utils, only: ErrorMessage, RotateVectors, cross_product, split, ReadStringFromFile, OpenFile, CloseFile
  use simcell, only: RandomPointInCell, ApplyBoundaryCondition, SetSimulationCellBoxLengthsAndBoxAngles
  use forcefield, only: ReadPairwiseInteractions, ReadIntraMolecularPotentials, &
    PotentialModel_info, Molecule_Intra, GetBondLength, GetBendAngle, &
    AngleInteraction, TorsionInteraction, TwoBodyNonBondedInteraction, GetRosenbluthWeightForBond, &
    GetRosenbluthWeightForAngle, IntraPairInteraction, IntraCoulInteraction
  implicit none
  private
  save

  integer, parameter :: NumberOfTrialInsertions=10, NumberOfTorsionTrials=100, NumberOfB_AngleTrials=1000
  integer, parameter :: MAX_NO_OF_BRANCH_POINTS=4, MAX_INTRA_PAIRS=20

  type :: CBMC_Step_Info
    character(len=strlen) :: SamplingType='rgs'
    integer :: NumberOfAtoms=1
    integer :: GroupNumber=1
    integer, dimension(MAX_NO_OF_BRANCH_POINTS) :: NumberOfTorsions=0
    integer, dimension(MAX_NO_OF_BRANCH_POINTS) :: NumberOfIntraPairs=0
    integer, dimension(MAX_NO_OF_BRANCH_POINTS) :: NumberOfIntraCouls=0

    integer, dimension(MAX_NO_OF_BRANCH_POINTS) :: Atom=0
    integer, dimension(MAX_NO_OF_BRANCH_POINTS) :: BondNumber=0
    integer, dimension(MAX_NO_OF_BRANCH_POINTS) :: A_AngleNumber=0
    integer, dimension(MAX_NO_OF_BRANCH_POINTS,MAX_NO_OF_BRANCH_POINTS) :: B_AngleNumber=0
    integer, dimension(MAX_NO_OF_BRANCH_POINTS,MAX_NO_OF_BRANCH_POINTS) :: TorsionNumber=0
    integer, dimension(MAX_NO_OF_BRANCH_POINTS,MAX_INTRA_PAIRS) :: IntraPairNumber=0
    integer, dimension(MAX_NO_OF_BRANCH_POINTS,MAX_INTRA_PAIRS) :: IntraCoulNumber=0

    integer :: Atom_p1=0, Atom_p2=0
  end type CBMC_Step_Info

  type :: CBMC_Growth_Info
    integer :: NumberOfSteps=1
    real(PR) :: IdealRosenbluthWeight=1._PR
    type(CBMC_Step_Info), dimension(:), allocatable :: Step
  end type CBMC_Growth_Info

  type :: CBMC_Move_Set
    integer :: NumberOfGrowthSequences=1
    type(CBMC_Growth_Info), dimension(:), allocatable :: GrowthSequence
  end type CBMC_Move_Set

  type(CBMC_Move_Set), dimension(MAX_NO_OF_SPECIES) :: CBMC_Move

  !=====================================================================================================================
  public :: CBMC_Move
  public :: CalculateIdealGasRosenbluthWeight, ReadIdealGasRosenbluthWeight
  public :: ReadCBMC_Moves, CBMC_grow, CBMC_regrow

contains
  !=================================================================================================
  subroutine CalculateIdealGasRosenbluthWeight(PROJECT_DIR,SimulationInput,LogUnitNo,error)
    character(len=lstrlen), intent(in)               :: PROJECT_DIR
    character(len=lstrlen), dimension(:), intent(in) :: SimulationInput
    integer, intent(in)                              :: LogUnitNo
    logical, intent(Out)                             :: error

    integer :: MoveType, CycleLength, ParticleMoveWeight
    logical :: overlap, EquilibrationStage, RestartMode
    character(len=lstrlen) :: string
    type(StorageInteractions) :: Interactions
    real(PR) :: energy, rosen, rosen_sum, rosen_avg, rosensq_sum, rosensq_avg, rosen_error
    integer :: iseed

    integer :: unitno, simno, cycleno, lineno, sys, spc, mol
    integer :: StartLine, SequenceNumber
    integer :: nmoles, natoms
    real(PR) :: rvalue
    real(pr), dimension(3) :: boxlength, boxangle
    logical :: existing_molecule
    

    error=.false.
    !---------------------------------------------------------------------------------------------------
    !** Reduced Units
    call SetReducedUnits(SimulationInput)

    !---------------------------------------------------------------------------------------------------
    !** Set Random Number Seed
    call ReadStringFromFile(string,'random number seed',SimulationInput,'Simulation Input',lineno,error=error)
    if(error)return
    read(string,*)iseed
    call SetRandomNumberSeed(iseed)

    !---------------------------------------------------------------------------------------------------
    !** Atoms
    call ReadAtomInfo(PROJECT_DIR,LogUnitNo,error)
    if(error)return

    !---------------------------------------------------------------------------------------------------
    !** Molecules
    call ReadMoleculeInfo(PROJECT_DIR,SimulationInput,error)
    if(error)return
    !---------------------------------------------------------------------------------------------------
    !** CBMC Moves
    call ReadCBMC_Moves(PROJECT_DIR,LogUnitNo,error)
    if(error)return

    !---------------------------------------------------------------------------------------------------
    !** Forcefield
    call ReadPairwiseInteractions(PROJECT_DIR,LogUnitNo,error)
    if(error)return
    call ReadIntraMolecularPotentials(PROJECT_DIR,LogUnitNo,error)
    if(error)return
    call AllocateMemoryToSeparationVectors(error)
    if(error)return

    !----------------------------------------------------------------------------------------------------
    !** Number Of Simulations
    call ReadStringFromFile(string,'number of simulations',SimulationInput,'Simulation Input',lineno,error=error)
    if(error)return
    read(string,*)NumberOfSimulations
    if(NumberOfSimulations > MAX_NO_OF_SIMULATIONS)then
      write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
        'Number of simulations exceeds MAX_NO_OF_SIMULATIONS:',MAX_NO_OF_SIMULATIONS
      error=.true.
      return
    end if

    !----------------------------------------------------------------------------------------------------
    !** Number Of Production Cycles
    call ReadStringFromFile(string,'number of production cycles',SimulationInput,'Simulation Input',lineno,error=error)
    if(error)return
    read(string,*)NumberOfProductionCycles

    !----------------------------------------------------------------------------------------------------
    !** Display Frequency
    call ReadStringFromFile(string,'display frequency',SimulationInput,'Simulation Input',lineno,error=error)
    if(error)return
    read(string,*)DisplayFrequency

    !----------------------------------------------------------------------------------------------------
    !**       Simulation Cells
    !----------------------------------------------------------------------------------------------------
    sys=1
    boxlength=(/1000._PR, 1000._PR, 1000._PR/)
    boxangle=(/90._PR, 90._PR, 90._PR/)*PI/180._PR
    call SetSimulationCellBoxLengthsAndBoxAngles(CurrentSimulationCell(sys),boxlength,boxangle)

    !----------------------------------------------------------------------------------------------------
    !**       State Variables
    !----------------------------------------------------------------------------------------------------
    !** Read Temperature
    StartLine=1
    do
      call ReadStringFromFile(string,'temperature',SimulationInput(StartLine:),'Simulation Input',lineno)
      if(lineno == 0)exit
      read(string,*)simno,sys,rvalue
      if(simno > NumberOfSimulations)then
        write(ErrorMessage,'(2a,i5,4x,a,i2,a,i2)')__FILE__,':',__LINE__, &
          'Simulation number, ',simno,', exceeds Number of Simulations :',NumberOfSimulations
        error=.true.
        return
      else if(sys > 1 )then
        write(ErrorMessage,'(2a,i5,4x,a)')__FILE__,':',__LINE__, &
          'Only one System allowed.'
        error=.true.
        return
      end if
      Temperature(simno,sys)=rvalue
      if(REDUCED_UNITS)then
        beta(simno,sys)=1._PR/rvalue
      else
        beta(simno,sys)=1._PR/(K_B*rvalue)
      end if
      StartLine=StartLine+lineno
    end do
      
    sys=1
    spc=1
    mol=1
    existing_molecule=.false.

    nmoles=1
    natoms=Molecule(spc)%NumberOfAtoms
    !** Allocate memory to coordinates
    call AllocateMemoryForCoordinates(CurrentCoordinates(spc,sys),nmoles,natoms,error)
    if(error)return
    call AllocateMemoryForCoordinates(TrialCoordinates(spc,sys),nmoles,natoms,error)
    if(error)return

    call OpenFile('IdealGasRosenbluthWeight.out','write','rewind',unitno,error)
    if(error)return
    write(unitno,'(4a20)')'Temperature','Sequence Number','Rosenbluth Weight','Error'
    write(unitno,'(a)')trim(Dashed_Line)
    do simno=1,NumberOfSimulations
      do SequenceNumber=1,CBMC_Move(spc)%NumberOfGrowthSequences
        !** Production Stage **
        rosen_sum=0._PR
        do cycleno=1,NumberOfProductionCycles
          rosen=CBMC_grow(SequenceNumber,sys,spc,mol,beta(simno,sys),existing_molecule,overlap)
          !** Sample Averages
          if(.not. overlap)then
            rosen_sum=rosen_sum+rosen
            rosensq_sum=rosensq_sum+rosen**2
          end if

          if(mod(cycleno,DisplayFrequency)==0)then
            write(LogUnitNo,*)cycleno,rosen
            call flush(LogUnitNo)
          end if
        end do
        rosen_avg=rosen_sum/real(NumberOfProductionCycles,PR)
        rosensq_avg=rosensq_sum/real(NumberOfProductionCycles,PR)
        rosen_error=2._PR*sqrt((rosensq_avg-rosen_avg**2)/real(NumberOfProductionCycles-1,PR))
        write(unitno,'(f20.4,i20,2es20.4)')Temperature(simno,sys), SequenceNumber, rosen_avg,rosen_error
      end do
    end do
    call CloseFile(unitno,error)
    if(error)return

  end subroutine CalculateIdealGasRosenbluthWeight

  !=====================================================================================================================
  subroutine ReadCBMC_Moves(PROJECT_DIR,LogUnitNo,error)
    character(len=lstrlen), intent(in) :: PROJECT_DIR
    integer, intent(in)                :: LogUnitNo
    logical, intent(out)               :: error

    character(len=lstrlen) :: filename
    type(CBMC_Move_Set)    :: Move
    type(CBMC_Growth_Info) :: GrowthSequence
    type(CBMC_Step_Info)   :: Step
    integer :: spc, atomno, unitno, indx
    integer :: SequenceNumber, StepNumber, nsteps
    integer :: ierror
    character(len=lstrlen) :: line
    integer :: nfields_1,nfields_2
    character(len=strlen), dimension(20) :: field_1,field_2
    logical :: lopen

    error=.false.
    do spc=1,NumberOfSpecies
      filename=trim(PROJECT_DIR)//'/'//adjustl(trim(Molecule(spc)%Name))//'.mol'
      call OpenFile(filename,'read','rewind',unitno,error)
      if(error)return

      !** Determine Number of Growth Sequences and allocate space
      SequenceNumber=0
      do
        read(unit=unitno,fmt='(a)',iostat=ierror)line
        if(ierror == -1)exit
        line=adjustl(line)
        indx=index(line,'CBMC Sequence')
        if(indx == 1)then
          SequenceNumber=SequenceNumber+1
        end if
      end do
      if(SequenceNumber == 0)then
        !** For Rigid Molecules
        CBMC_Move(spc)%NumberOfGrowthSequences=1
        allocate(CBMC_Move(spc)%GrowthSequence(1),STAT=ierror)
        if(ierror /=0)then
          error=.true.
          write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
            'Failed to allocate memory. STAT = ',ierror
          return
        end if
        CBMC_Move(spc)%GrowthSequence(1)%NumberOfSteps=1
        allocate(CBMC_Move(spc)%GrowthSequence(1)%Step(1),STAT=ierror)
        if(ierror /=0)then
          error=.true.
          write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
            'Failed to allocate memory. STAT = ',ierror
          return
        end if
        CBMC_Move(spc)%GrowthSequence(1)%Step(1)%SamplingType='rgs'
        CBMC_Move(spc)%GrowthSequence(1)%Step(1)%GroupNumber=1
      else
        !** For Flexible Molecules
        Move%NumberOfGrowthSequences=SequenceNumber
        allocate(Move%GrowthSequence(Move%NumberOfGrowthSequences),STAT=ierror)
        if(ierror /=0)then
          error=.true.
          write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
            'Failed to allocate memory. STAT = ',ierror
          return
        end if

        !** Read Growth Sequences
        rewind(unitno)
        SequenceNumber=0
        do
          read(unit=unitno,fmt='(a)',iostat=ierror)line
          if(ierror == -1)exit
          line=adjustl(line)
          indx=index(line,'CBMC Sequence')
          if(indx == 1)then
            SequenceNumber=SequenceNumber+1
            read(unitno,*)nsteps
            GrowthSequence%NumberOfSteps=nsteps
            allocate(GrowthSequence%Step(nsteps),STAT=ierror)
            if(ierror /=0)then
              error=.true.
              write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
                'Failed to allocate memory. STAT = ',ierror
              return
            end if

            do StepNumber=1,GrowthSequence%NumberOfSteps
              read(unitno,'(a)')line
              nfields_1=split(line,field_1,' ')
              read(field_1(1),*)Step%SamplingType
              select case (Step%SamplingType)
              case ('rgs')
                !** Rigid Group Sampling
                read(field_1(2),*)Step%GroupNumber
              case ('fas')
                !** First Atom Sampling
                Step%NumberOfAtoms=1
                read(field_1(2),*)Step%Atom(1)
              case ('ss')
                !** Sphere Sampling
                Step%NumberOfAtoms=1
                read(field_1(2),*)Step%Atom(1)
                read(field_1(3),*)Step%Atom_p1
                read(field_1(4),*)Step%BondNumber(1)
              case ('bps')
                read(field_1(2),*)Step%NumberOfAtoms
                read(field_1(3),*)Step%Atom_p1
                read(field_1(4),*)Step%Atom_p2
                do atomno=1,Step%NumberOfAtoms
                  read(unitno,'(a)')line
                  nfields_2=split(line,field_2)
                  read(field_2(1),*)Step%Atom(atomno)
                  read(field_2(2),*)Step%BondNumber(atomno)
                  read(field_2(3),*)Step%A_AngleNumber(atomno)

                  read(field_2(4),*)Step%B_AngleNumber(atomno,1:atomno-1)
                  read(field_2(5),*)Step%NumberOfTorsions(atomno)
                  read(field_2(6),*)Step%TorsionNumber(atomno,1:Step%NumberOfTorsions(atomno))
                  read(field_2(7),*)Step%NumberOfIntraPairs(atomno)
                  read(field_2(8),*)Step%IntraPairNumber(atomno,1:Step%NumberOfIntraPairs(atomno))
                  read(field_2(9),*)Step%NumberOfIntraCouls(atomno)
                  read(field_2(10),*)Step%IntraCoulNumber(atomno,1:Step%NumberOfIntraCouls(atomno))
                end do
              case default
                error=.true.
                write(ErrorMessage,'(2a,i5,4x,2a)')__FILE__,':',__LINE__, &
                  'Unknown CBMC Step Sampling Type= ',adjustl(trim(Step%SamplingType))
                return
              end select
              GrowthSequence%Step(StepNumber)=Step
              call ResetCBMC_Step_Info(Step)
            end do
            Move%GrowthSequence(SequenceNumber)=GrowthSequence
            deallocate(GrowthSequence%Step,STAT=ierror)
            if(ierror /=0)then
              error=.true.
              write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
                'Failed to deallocate memory. STAT = ',ierror
              return
            end if
          else
          end if
        end do
        CBMC_Move(spc)=Move
        deallocate(Move%GrowthSequence,STAT=ierror)
        if(ierror /=0)then
          error=.true.
          write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
            'Failed to deallocate memory. STAT = ',ierror
          return
        end if
      end if

      call CloseFile(unitno,error)
      if(error)return
    end do

  end subroutine ReadCBMC_Moves

  !=====================================================================================================================
  subroutine ReadIdealGasRosenbluthWeight(SimulationInput,error)
    character(len=lstrlen), dimension(:), intent(in) :: SimulationInput
    logical, intent(out) :: error

    character(len=lstrlen) :: string='',filename=''
    integer :: StartLine, lineno
    integer :: sys, spc
    integer :: simno,seqno
    character(len=strlen) :: spcname=''
    real(PR) :: rvalue

    error=.false.

    !** Read Ideal Gas Rosenbluth Weight
    StartLine=1
    do
      call ReadStringFromFile(string,'IdealRFactor',SimulationInput(StartLine:),'Simulation Input',lineno)
      if(lineno == 0)exit
      read(string,*)spcname,simno,sys,seqno,rvalue
      if(simno > NumberOfSimulations)then
        write(ErrorMessage,'(2a,i5,4x,a,i2,a,i2)')__FILE__,':',__LINE__, &
          'Simulation number, ',simno,', exceeds Number of Simulations :',NumberOfSimulations
        error=.true.
        return
      else if(sys > NumberOfSystems)then
        write(ErrorMessage,'(2a,i5,4x,a,i2,a,i2)')__FILE__,':',__LINE__, &
          'System number, ',sys,', exceeds Number of Systems :',NumberOfSystems
        error=.true.
        return
      else if(seqno > CBMC_Move(spc)%NumberOfGrowthSequences)then
        write(ErrorMessage,'(2a,i5,4x,a,i2,3a,i2)')__FILE__,':',__LINE__, &
          'Sequence number, ',seqno,', for species ',trim(Molecule(spc)%Name), &
          'exceeds Number of Growth Sequences: ',CBMC_Move(spc)%NumberOfGrowthSequences
        error=.true.
        return
      end if
      CBMC_Move(spc)%GrowthSequence(seqno)%IdealRosenbluthWeight=rvalue
      StartLine=StartLine+lineno
    end do

  end subroutine ReadIdealGasRosenbluthWeight

  !=====================================================================================================================
  function CBMC_grow(SequenceNumber,sys,spc,mol,beta,existing_molecule,overlap) result(RosenbluthWeight)
    integer, intent(in)                            :: SequenceNumber
    integer, intent(in)                            :: sys,spc,mol
    real(PR), intent(in)                           :: beta
    logical, intent(in)                            :: existing_molecule
    logical, intent(out)                           :: overlap
    real(PR)                                       :: RosenbluthWeight

    integer :: StepNumber
    type(CBMC_Step_Info) :: Step

    overlap=.false.
    !** Copy the current coordinates to trial coordinates.
    TrialCoordinates(sys,spc)%Positions(:,mol,:)=CurrentCoordinates(sys,spc)%Positions(:,mol,:)
    TrialCoordinates(sys,spc)%CenterOfMass(:,mol)=CurrentCoordinates(sys,spc)%CenterOfMass(:,mol)

    RosenbluthWeight=1._PR
    do StepNumber=1,CBMC_Move(spc)%GrowthSequence(SequenceNumber)%NumberOfSteps
      Step=CBMC_Move(spc)%GrowthSequence(SequenceNumber)%Step(StepNumber)
      select case(trim(Step%SamplingType))
      case ('rgs')
        RosenbluthWeight=RosenbluthWeight*RigidGroupSampling(Step,sys,spc,mol,beta,existing_molecule,overlap)
      case ('fas')
        RosenbluthWeight=RosenbluthWeight*FirstAtomSampling(Step,sys,spc,mol,beta,existing_molecule,overlap)
      case ('ss')
        RosenbluthWeight=RosenbluthWeight*SphereSampling(Step,sys,spc,mol,beta,existing_molecule,overlap)
      case ('bps')
        RosenbluthWeight=RosenbluthWeight*BranchPointSampling(Step,sys,spc,mol,beta,existing_molecule,overlap)
      end select
      if(overlap)return
    end do
    TrialCoordinates(sys,spc)%CenterOfMass(:,mol)=GetCenterOfMass(spc,TrialCoordinates(sys,spc)%Positions(:,mol,:))
  end function CBMC_grow

  !=====================================================================================================================
  function CBMC_regrow(SequenceNumber,sys,spc,mol,beta,overlap,error) result(RosenbluthWeightRatio)
    integer, intent(in)    :: SequenceNumber
    integer, intent(in)    :: sys,spc,mol
    real(PR), intent(in)   :: beta
    logical, intent(out)   :: overlap
    logical, intent(out)   :: error
    real(PR)               :: RosenbluthWeightRatio

    logical :: existing_molecule
    integer :: StepNumber, StartStepNumber
    real(PR) :: RosenbluthWeight_old,RosenbluthWeight_new
    type(CBMC_Step_Info) :: Step

    error=.false.
    !** Copy the current coordinates to trial coordinates.
    TrialCoordinates(sys,spc)%Positions(:,mol,:)=CurrentCoordinates(sys,spc)%Positions(:,mol,:)
    TrialCoordinates(sys,spc)%CenterOfMass(:,mol)=CurrentCoordinates(sys,spc)%CenterOfMass(:,mol)

    !** Randomly choose a StartStepNumber
    StartStepNumber=int(RandomNumber()*CBMC_Move(spc)%GrowthSequence(SequenceNumber)%NumberOfSteps)+1

    existing_molecule=.true.
    RosenbluthWeight_old=1._PR
    do StepNumber=StartStepNumber,CBMC_Move(spc)%GrowthSequence(SequenceNumber)%NumberOfSteps
      Step=CBMC_Move(spc)%GrowthSequence(SequenceNumber)%Step(StepNumber)
      select case(Step%SamplingType)
      case ('rgs')
        RosenbluthWeight_old=RosenbluthWeight_old*RigidGroupSampling(Step,sys,spc,mol,beta,existing_molecule,overlap)
      case ('fas')
        RosenbluthWeight_old=RosenbluthWeight_old*FirstAtomSampling(Step,sys,spc,mol,beta,existing_molecule,overlap)
      case ('ss')
        RosenbluthWeight_old=RosenbluthWeight_old*SphereSampling(Step,sys,spc,mol,beta,existing_molecule,overlap)
      case ('bps')
        RosenbluthWeight_old=RosenbluthWeight_old*BranchPointSampling(Step,sys,spc,mol,beta,existing_molecule,overlap)
      end select
      if(overlap)then
        write(ErrorMessage,'(2a,i4,4x,a)')__FILE__,':',__LINE__, &
          'Something wrong! Existing configuration contains overlapping atoms'
        error=.true.
        return
      end if
    end do

    existing_molecule=.false.
    RosenbluthWeight_new=1._PR
    do StepNumber=StartStepNumber,CBMC_Move(spc)%GrowthSequence(SequenceNumber)%NumberOfSteps
      Step=CBMC_Move(spc)%GrowthSequence(SequenceNumber)%Step(StepNumber)
      select case(Step%SamplingType)
      case ('rgs')
        RosenbluthWeight_new=RosenbluthWeight_new*RigidGroupSampling(Step,sys,spc,mol,beta,existing_molecule,overlap)
      case ('fas')
        RosenbluthWeight_new=RosenbluthWeight_new*FirstAtomSampling(Step,sys,spc,mol,beta,existing_molecule,overlap)
      case ('ss')
        RosenbluthWeight_new=RosenbluthWeight_new*SphereSampling(Step,sys,spc,mol,beta,existing_molecule,overlap)
      case ('bps')
        RosenbluthWeight_new=RosenbluthWeight_new*BranchPointSampling(Step,sys,spc,mol,beta,existing_molecule,overlap)
      end select
      if(overlap)return
    end do

    TrialCoordinates(sys,spc)%CenterOfMass(:,mol)=GetCenterOfMass(spc,TrialCoordinates(sys,spc)%Positions(:,mol,:))
    RosenbluthWeightRatio=RosenbluthWeight_new/RosenbluthWeight_old
  end function CBMC_regrow

  !=====================================================================================================================
  function RigidGroupSampling(Step,sys,spc,mol,beta,existing_molecule,overlap) result(RosenbluthWeight)
    type(CBMC_Step_Info), intent(in)         :: Step
    integer, intent(in)                      :: sys, spc, mol
    real(PR), intent(in)                     :: beta
    logical, intent(in)                      :: existing_molecule
    logical, intent(out)                     :: overlap
    real(PR)                                 :: RosenbluthWeight

    logical :: Success
    integer :: grp,natoms
    integer :: TrialNumber
    real(PR), dimension(3) :: vec
    real(PR), dimension(3,3) :: RotationMatrix
    real(PR), dimension(:,:,:), allocatable :: TrialPosition
    real(PR), dimension(NumberOfTrialInsertions) :: rosen
    type(StorageInteractions) :: TrialInteractions
    
    grp=Step%GroupNumber
    natoms=Molecule(spc)%Group(grp)%NumberOfAtoms
    allocate(TrialPosition(3,natoms,NumberOfTrialInsertions))

    Success=.false.
    !** First Trial: Current coordinates for deletion and new coordinates for insertion
    TrialNumber=1
    if(existing_molecule)then
      vec=TrialCoordinates(spc,sys)%CenterOfMass(:,mol)
      TrialPosition(:,:,TrialNumber)=TrialCoordinates(spc,sys)%Positions(:,mol,:)
    else
      !** Get a random point within the cell for the Center of Mass
      call RandomPointInCell(CurrentSimulationCell(sys),vec)
      TrialPosition(:,:,TrialNumber)=Molecule(spc)%Group(grp)%ReferencePosition+spread(vec,2,natoms)
      if(natoms > 1)then
        call UniformRandomRotationMatrix(RotationMatrix)
        call RotateVectors(TrialPosition(:,:,TrialNumber),RotationMatrix,(/0._PR, 0._PR, 0._PR/))
        call ApplyBoundaryCondition(CurrentSimulationCell(sys),natoms,TrialPosition(:,:,TrialNumber), vec)
      end if
    end if
    !** Calculate Energy from Short Range Interactions
    call MoleculeSystemShortRangePairwiseInteraction(TrialPosition(1:3,1:natoms,TrialNumber),vec, &
        CurrentCoordinates(:,sys),CurrentSimulationCell(sys),MainCellList(sys),spc,mol,TrialInteractions,overlap)
    if(overlap)then
      rosen(TrialNumber)=0._PR
      if(existing_molecule)return        !** Something wrong as existing atoms are overlapping
    else
      Success=.true.
      call UpdateStorageInteractions(TrialInteractions)
      rosen(TrialNumber)=exp(-beta*TrialInteractions%Energy%Total)
    end if
    if(natoms == 1)then
      TrialCoordinates(spc,sys)%Positions(:,mol,:)=TrialPosition(:,:,TrialNumber)
      RosenbluthWeight=rosen(TrialNumber)
      return
    end if

    !** Rest of the trials
    do TrialNumber=2,NumberOfTrialInsertions
      TrialPosition(:,:,TrialNumber)=Molecule(spc)%Group(grp)%ReferencePosition+spread(vec,2,natoms)
      call UniformRandomRotationMatrix(RotationMatrix)
      call RotateVectors(TrialPosition(:,:,TrialNumber),RotationMatrix,vec)
      call ApplyBoundaryCondition(CurrentSimulationCell(sys),natoms,TrialPosition(:,:,TrialNumber), vec)

      !** Calculate Energy from Short Range Interactions
      call MoleculeSystemShortRangePairwiseInteraction(TrialPosition(1:3,1:natoms,TrialNumber),vec, &
        CurrentCoordinates(:,sys),CurrentSimulationCell(sys),MainCellList(sys),spc,mol,TrialInteractions,overlap)
      if(overlap)then
        rosen(TrialNumber)=0._PR
      else
        Success=.true.
        call UpdateStorageInteractions(TrialInteractions)
        rosen(TrialNumber)=exp(-beta*TrialInteractions%Energy%Total)
      end if 
    end do
    overlap=(.not. Success)
    if(overlap)return

    !** Choose a Trial Position
    if(existing_molecule)then
      TrialNumber=1
    else
      TrialNumber=ChosenTrialNumber(NumberOfTrialInsertions,rosen)
      TrialCoordinates(spc,sys)%Positions(1:3,mol,1:natoms)=TrialPosition(1:3,1:natoms,TrialNumber)
    end if
    RosenbluthWeight=sum(rosen)/real(NumberOfTrialInsertions,PR)

  end function RigidGroupSampling

  !=====================================================================================================================
  function FirstAtomSampling(Step,sys,spc,mol,beta,existing_molecule,overlap) result(RosenbluthWeight)
    type(CBMC_Step_Info), intent(in)         :: Step
    integer, intent(in)                      :: sys, spc, mol
    real(PR), intent(in)                     :: beta
    logical, intent(in)                      :: existing_molecule
    logical, intent(out)                     :: overlap
    real(PR)                                 :: RosenbluthWeight

    integer :: atm
    integer :: TrialNumber
    real(PR), dimension(3) :: vec
    real(PR), dimension(NumberOfTrialInsertions) :: rosen
    type(StorageInteractions) :: TrialInteractions
    
    atm=Step%Atom(1)

    if(existing_molecule)then
      vec=TrialCoordinates(spc,sys)%Positions(:,mol,atm)
    else
      !** Get a random point within the cell for the Center of Mass
      call RandomPointInCell(CurrentSimulationCell(sys),vec)
      TrialCoordinates(spc,sys)%Positions(:,mol,atm)=vec
    end if
    !** Calculate Energy from Short Range Interactions
    call AtomSystemShortRangePairwiseInteraction(vec,vec,CurrentCoordinates(:,sys),CurrentSimulationCell(sys), &
      MainCellList(sys),spc,mol,atm,TrialInteractions,overlap)
    if(overlap)return
    call UpdateStorageInteractions(TrialInteractions)
    RosenbluthWeight=exp(-beta*TrialInteractions%Energy%Total)

  end function FirstAtomSampling

  !=====================================================================================================================
  function SphereSampling(Step,sys,spc,mol,beta,existing_molecule,overlap) result(RosenbluthWeight)
    type(CBMC_Step_Info), intent(in)         :: Step
    integer, intent(in)                      :: sys, spc, mol
    real(PR), intent(in)                     :: beta
    logical, intent(in)                      :: existing_molecule
    logical, intent(out)                     :: overlap
    real(PR)                                 :: RosenbluthWeight

    logical :: Success
    integer :: atm
    integer :: TrialNumber
    real(PR) :: BondLength
    real(PR), dimension(3) :: vec, Position_p1
    real(PR), dimension(3,NumberOfTrialInsertions) :: TrialPosition
    real(PR), dimension(NumberOfTrialInsertions) :: rosen
    type(StorageInteractions) :: TrialInteractions
    
    atm=Step%Atom(1)
    Position_p1=TrialCoordinates(spc,sys)%Positions(:,mol,Step%Atom_p1)
    if(existing_molecule)then
      vec=TrialCoordinates(spc,sys)%Positions(:,mol,atm)-Position_p1
      BondLength=norm2(vec)
    else
      BondLength=GetBondLength(Molecule_intra(spc)%Bond(Step%BondNumber(1)),beta)
    end if
    RosenbluthWeight=GetRosenbluthWeightForBond(Molecule_intra(spc)%Bond(Step%BondNumber(1)),beta,Bondlength)

    Success=.false.
    !** First Trial: Current coordinates for deletion and new coordinates for insertion
    TrialNumber=1
    if(existing_molecule)then
      vec=TrialCoordinates(spc,sys)%Positions(:,mol,atm)
    else
      !** Get a random point on surface of sphere
      call RandomUnitVector(vec)
      vec=Position_p1+BondLength*vec
    end if
    TrialPosition(:,TrialNumber)=vec
    !** Calculate Energy from Short Range Interactions
    call AtomSystemShortRangePairwiseInteraction(vec,vec,CurrentCoordinates(:,sys),CurrentSimulationCell(sys), &
      MainCellList(sys),spc,mol,atm,TrialInteractions,overlap)
    if(overlap)then
      rosen(TrialNumber)=0._PR
      if(existing_molecule)return        !** Something wrong as existing atoms are overlapping
    else
      Success=.true.
      call UpdateStorageInteractions(TrialInteractions)
      rosen(TrialNumber)=exp(-beta*TrialInteractions%Energy%Total)
    end if

    !** Rest of the trials
    do TrialNumber=2,NumberOfTrialInsertions
      !** Get a random point on surface of sphere
      call RandomUnitVector(vec)
      vec=Position_p1+BondLength*vec
      TrialPosition(:,TrialNumber)=vec

      !** Calculate Energy from Short Range Interactions
      call AtomSystemShortRangePairwiseInteraction(vec,vec,CurrentCoordinates(:,sys),CurrentSimulationCell(sys), &
        MainCellList(sys),spc,mol,atm,TrialInteractions,overlap)
      if(overlap)then
        rosen(TrialNumber)=0._PR
      else
        Success=.true.
        call UpdateStorageInteractions(TrialInteractions)
        rosen(TrialNumber)=exp(-beta*TrialInteractions%Energy%Total)
      end if
    end do
    overlap=(.not. Success)
    if(overlap)return

    !** Choose a Trial Position
    if(existing_molecule)then
      TrialNumber=1
    else
      TrialNumber=ChosenTrialNumber(NumberOfTrialInsertions,rosen)
      TrialCoordinates(spc,sys)%Positions(1:3,mol,atm)=TrialPosition(:,TrialNumber)
    end if
    RosenbluthWeight=RosenbluthWeight*sum(rosen)/real(NumberOfTrialInsertions,PR)

  end function SphereSampling

  !=====================================================================================================================
  function BranchPointSampling(Step,sys,spc,mol,beta,existing_molecule,overlap) result(RosenbluthWeight)
    type(CBMC_Step_Info), intent(in)         :: Step
    integer, intent(in)                      :: sys, spc, mol
    real(PR), intent(in)                     :: beta
    logical, intent(in)                      :: existing_molecule
    logical, intent(out)                     :: overlap
    real(PR)                                 :: RosenbluthWeight

    logical :: Success
    integer :: i, j, TrialNumber, TorsionTrialNumber
    integer :: atm, a1, a2, a4
    real(PR) :: rijsq
    real(PR) :: energy, B_Theta, phi, TorsionRosenbluthWeight
    type(PotentialModel_info) :: Torsion, IntraPair
    real(PR), dimension(3) :: vec, vecZ, vecX, vecY, vecX1, vecY1
    real(PR), dimension(3) :: Position_p1, Position_p2, Position_p3, r3
    real(PR), dimension(MAX_NO_OF_BRANCH_POINTS) :: Bondlength, Theta, phi_B
    real(PR), dimension(3,MAX_NO_OF_BRANCH_POINTS) :: r2
    real(PR), dimension(NumberOfB_AngleTrials) :: rosenB, phi_trial
    real(PR), dimension(NumberOfTorsionTrials) :: rosenT
    real(PR), dimension(NumberOfTrialInsertions) :: rosen, phi_T
    type(StorageInteractions) :: TrialInteractions
    type(InteractionAccumulator) :: Accumulator
    logical :: underflow

    !** Coordinates of previous two atoms
    Position_p1=TrialCoordinates(spc,sys)%Positions(:,mol,Step%Atom_p1)
    Position_p2=TrialCoordinates(spc,sys)%Positions(:,mol,Step%Atom_p2)
    !** Get a set of orthogonal vectors, i.e. vecX, vecY, vecZ
    vecZ=Position_p1-Position_p2
    vecZ=vecZ/norm2(vecZ)
    if(abs(vecZ(1)) > 1e-5_pr .or. abs(vecZ(2)) > 1e-5_pr)then
      vecX(1)=-vecZ(2)
      vecX(2)=vecZ(1)
      vecX(3)=0._PR !vecZ(3)
    else
      vecX(1)=-vecZ(3)
      vecX(2)=0._PR !vecZ(2)
      vecX(3)=vecZ(1)
    end if
    vecX=vecX/norm2(vecX)
    vecY=cross_product(vecZ,vecX)

    !**Obtain Bondlength and theta for 1st atom
    !----------------------------------------------------------------------------------------
    if(existing_molecule)then
      !** Calculate bond length and theta from current position
      vec=TrialCoordinates(spc,sys)%Positions(:,mol,Step%Atom(1))-Position_p1 
      Bondlength(1)=norm2(vec)
      Theta(1)=acos(dot_product(vec,-vecZ))
      !** Get a set of orthogonal vectors, i.e. vecX1 and vecY1 such that phi_B(1) = 0
      vecY1=cross_product(vecZ,vec)
      vecX1=cross_product(vecY,vecZ)
    else
      !** Generate bond length and theta from intramolecular potential
      Bondlength(1)=GetBondLength(Molecule_intra(spc)%Bond(Step%BondNumber(1)),beta)
      Theta(1)=GetBendAngle(Molecule_intra(spc)%Angle(Step%A_AngleNumber(1)),beta)
    end if
    RosenbluthWeight=GetRosenbluthWeightForBond(Molecule_intra(spc)%Bond(Step%BondNumber(1)),beta,Bondlength(1))
    RosenbluthWeight=RosenbluthWeight*GetRosenbluthWeightForAngle(Molecule_intra(spc)%Angle(Step%A_AngleNumber(1)), &
      beta,Theta(1))

    !** Sample B_Angles and calculate Rosenbluth Weight
    !----------------------------------------------------------------------------------------
    phi_B(1)=0._PR
    r2(:,1)=(/sin(Theta(1)), 0._PR, -cos(Theta(1)) /)
    !** Branch Point Sampling
    do i=2,Step%NumberOfAtoms
      TrialNumber=1
      if(existing_molecule)then
        !** Get bond length, theta and phi from current position
        vec=TrialCoordinates(spc,sys)%Positions(:,mol,Step%Atom(i))-Position_p1 
        Bondlength(i)=norm2(vec)
        Theta(i)=acos(dot_product(vec,-vecZ))
        phi_trial(TrialNumber)=atan2(dot_product(vec,vecY1),dot_product(vec,vecX1))
      else
        !** Generate bond length and theta from intramolecular potential
        Bondlength(i)=GetBondLength(Molecule_intra(spc)%Bond(Step%BondNumber(i)),beta)
        Theta(i)=GetBendAngle(Molecule_intra(spc)%Angle(Step%A_AngleNumber(i)),beta)
        !** Choose phi uniformly between -pi an pi
        phi_trial(TrialNumber)=TWOPI*RandomNumber()-PI
      end if
      RosenbluthWeight=RosenbluthWeight*GetRosenbluthWeightForBond(Molecule_intra(spc)%Bond(Step%BondNumber(i)), &
        beta,Bondlength(i))
      RosenbluthWeight=RosenbluthWeight*GetRosenbluthWeightForAngle(Molecule_intra(spc)%Angle(Step%A_AngleNumber(i)), &
        beta,Theta(i))
      r3=(/sin(Theta(i))*cos(phi_trial(TrialNumber)), sin(Theta(i))*sin(phi_trial(TrialNumber)), -cos(Theta(i)) /)
      energy=0._PR
      do j=1,i-1
        B_Theta=acos(dot_product(r2(:,j),r3))
        energy=energy+AngleInteraction(Molecule_intra(spc)%Angle(Step%B_AngleNumber(i,j)),B_Theta)
      end do
      rosenB(TrialNumber)=exp(-beta*energy)
      call ieee_get_flag(ieee_underflow,underflow)
      if(underflow)then
        call ieee_set_flag(ieee_underflow,.false.)
        rosenB(TrialNumber)=0._PR
      end if

      do TrialNumber=2,NumberOfB_AngleTrials
        !** Choose phi uniformly between -pi an pi
        phi_trial(TrialNumber)=TWOPI*RandomNumber()-PI
        r3=(/sin(Theta(i))*cos(phi_trial(TrialNumber)), sin(Theta(i))*sin(phi_trial(TrialNumber)), -cos(Theta(i)) /)
        energy=0._PR
        do j=1,i-1
          B_Theta=acos(dot_product(r2(:,j),r3))
          energy=energy+AngleInteraction(Molecule_intra(spc)%Angle(Step%B_AngleNumber(i,j)),B_Theta)
        end do
        rosenB(TrialNumber)=exp(-beta*energy)
        call ieee_get_flag(ieee_underflow,underflow)
        if(underflow)then
          call ieee_set_flag(ieee_underflow,.false.)
          rosenB(TrialNumber)=0._PR
        end if
      end do
      if(existing_molecule)then
        TrialNumber=1
      else
        TrialNumber=ChosenTrialNumber(NumberOfB_AngleTrials,rosenB)
      end if
      phi_B(i)=phi_trial(TrialNumber)
      r2(:,i)=(/sin(Theta(i))*cos(phi_B(i)), sin(Theta(i))*sin(phi_B(i)), -cos(Theta(i)) /)
      RosenbluthWeight=RosenbluthWeight*sum(rosenB)/real(NumberOfB_AngleTrials)
    end do

    !** Sample Torsion Angles and compute Torsion Rosenbluth Weight for 1st Trial
    !----------------------------------------------------------------------------------------
    TrialNumber=1
    if(all(Step%NumberOfTorsions == 0))then
      phi_T(TrialNumber)=RandomNumber()*TWOPI
      TorsionRosenbluthWeight=1._PR
    else
      !** First Torsion Trial
      TorsionTrialNumber=1
      !** Choose an angle phi between -pi and pi
      if(.not. existing_molecule)phi_trial(TorsionTrialNumber)=RandomNumber()*TWOPI-PI
      !** Calculate Torsion Energy
      energy=0._PR
      do i=1,Step%NumberOfAtoms
        if(existing_molecule)then
          r3=TrialCoordinates(spc,sys)%Positions(:,mol,Step%Atom(i))-Position_p1 
        else
          r3=sin(Theta(i))*cos(phi_B(i)+phi_trial(TorsionTrialNumber))*vecX+ &
             sin(Theta(i))*sin(phi_B(i)+phi_trial(TorsionTrialNumber))*vecY- &
             cos(Theta(i))*vecZ
        end if
        do j=1,Step%NumberOfTorsions(i)
          Torsion=Molecule_intra(spc)%Torsion(Step%TorsionNumber(i,j))
          a1=Molecule_Intra(spc)%TorsionAtom(Step%TorsionNumber(i,j),1)
          a4=Molecule_Intra(spc)%TorsionAtom(Step%TorsionNumber(i,j),4)
          if(Step%Atom(i) == a1)then
            Position_p3=TrialCoordinates(spc,sys)%Positions(:,mol,a4)
          else
            Position_p3=TrialCoordinates(spc,sys)%Positions(:,mol,a1)
          end if
          vec=Position_p2-Position_p3
          vecY1=cross_product(vec,vecZ)
          vecY1=vecY1/norm2(vecY1)
          vecX1=cross_product(vecY1,vecZ)
          phi=atan2(dot_product(r3,vecY1),dot_product(r3,vecX1))
          energy=energy+TorsionInteraction(Torsion,phi)
        end do
      end do
      rosenT(TorsionTrialNumber)=exp(-beta*energy)
      
      !** Rest of Torsion Trials for First Trial
      do TorsionTrialNumber=2,NumberOfTorsionTrials
        !** Choose an angle phi between -pi and pi
        phi_trial(TorsionTrialNumber)=RandomNumber()*TWOPI-PI
        !** Calculate Torsion Energy
        energy=0._PR
        do i=1,Step%NumberOfAtoms
          r3=sin(Theta(i))*cos(phi_B(i)+phi_trial(TorsionTrialNumber))*vecX+ &
             sin(Theta(i))*sin(phi_B(i)+phi_trial(TorsionTrialNumber))*vecY- &
             cos(Theta(i))*vecZ
          do j=1,Step%NumberOfTorsions(i)
            Torsion=Molecule_intra(spc)%Torsion(Step%TorsionNumber(i,j))
            a1=Molecule_Intra(spc)%TorsionAtom(Step%TorsionNumber(i,j),1)
            a4=Molecule_Intra(spc)%TorsionAtom(Step%TorsionNumber(i,j),4)
            if(Step%Atom(i) == a1)then
              Position_p3=TrialCoordinates(spc,sys)%Positions(:,mol,a4)
            else
              Position_p3=TrialCoordinates(spc,sys)%Positions(:,mol,a1)
            end if
            vec=Position_p2-Position_p3
            vecY1=cross_product(vec,vecZ)
            vecY1=vecY1/norm2(vecY1)
            vecX1=cross_product(vecY1,vecZ)
            phi=atan2(dot_product(r3,vecY1),dot_product(r3,vecX1))
            energy=energy+TorsionInteraction(Torsion,phi)
          end do
        end do
        rosenT(TorsionTrialNumber)=exp(-beta*energy)
      end do
      !** Choose Torsion Angle
      if(.not. existing_molecule)then
        TorsionTrialNumber=ChosenTrialNumber(NumberOfTorsionTrials,rosenT)
        phi_T(TrialNumber)=phi_trial(TorsionTrialNumber)
      end if
      TorsionRosenbluthWeight=sum(rosenT)/real(NumberOfTorsionTrials,PR)
    end if
     
    Success=.false.
    !** Compute position and calculate LJ energy for First trial
    energy=0._PR
    loop1:do i=1,Step%NumberOfAtoms
      atm=Step%Atom(i)
      if(existing_molecule)then
        !** Use current position
        vec=TrialCoordinates(spc,sys)%Positions(:,mol,atm)
      else
        !** Compute Trial Position
        vec=sin(Theta(i))*cos(phi_B(i)+phi_T(TrialNumber))*vecX+ &
            sin(Theta(i))*sin(phi_B(i)+phi_T(TrialNumber))*vecY- &
            cos(Theta(i))*vecZ
        vec=Position_p1+vec*BondLength(i)
      end if
      !** Calculate Intra Pair Interactions
      do j=1,Step%NumberOfIntraPairs(i)
        IntraPair=Molecule_intra(spc)%IntraPair(Step%IntraPairNumber(i,j))
        a1=Molecule_intra(spc)%IntraPairAtom(Step%IntraPairNumber(i,j),1)
        a2=Molecule_intra(spc)%IntraPairAtom(Step%IntraPairNumber(i,j),2)
        if(atm == a2)then
          rijsq=sum((TrialCoordinates(spc,sys)%Positions(:,mol,a1)-vec)**2)
        else
          rijsq=sum((TrialCoordinates(spc,sys)%Positions(:,mol,a2)-vec)**2)
        end if
        energy=energy+IntraPairInteraction(Molecule_intra(spc)%IntraPair(i),rijsq,overlap)
        if(overlap)exit loop1
      end do

      !** Calculate Intermolecular Interactions
      call AtomSystemShortRangePairwiseInteraction(vec,vec,CurrentCoordinates(:,sys),CurrentSimulationCell(sys), &
        MainCellList(sys),spc,mol,atm,TrialInteractions,overlap)
      if(overlap)exit loop1
      call UpdateStorageInteractions(TrialInteractions)
      energy=energy+TrialInteractions%Energy%Total
    end do loop1
    !** RosenbluthWeight
    if(overlap)then
      rosen(TrialNumber)=0._PR
      if(existing_molecule)return                    !** something wrong existing atoms overlapping
    else
      Success=.true.
      rosen(TrialNumber)=exp(-beta*energy)*TorsionRosenbluthWeight
      call ieee_get_flag(ieee_underflow,underflow)
      if(underflow)then
        call ieee_set_flag(ieee_underflow,.false.)
        rosen(TrialNumber)=0._PR
      end if
    end if

    !** Sample Torsion Angles and compute Torsion Rosenbluth Weight for rest of the trials
    !----------------------------------------------------------------------------------------
    do TrialNumber=2,NumberOfTrialInsertions
      if(all(Step%NumberOfTorsions == 0))then
        phi_T(TrialNumber)=RandomNumber()*TWOPI
        TorsionRosenbluthWeight=1._PR
      else
        do TorsionTrialNumber=1,NumberOfTorsionTrials
          !** Choose an angle phi between -pi and pi
          phi_trial(TorsionTrialNumber)=RandomNumber()*TWOPI-PI
          !** Calculate Torsion Energy
          energy=0._PR
          do i=1,Step%NumberOfAtoms
            atm=Step%Atom(i)
            r3=sin(Theta(i))*cos(phi_B(i)+phi_trial(TorsionTrialNumber))*vecX+ &
               sin(Theta(i))*sin(phi_B(i)+phi_trial(TorsionTrialNumber))*vecY- &
               cos(Theta(i))*vecZ
            do j=1,Step%NumberOfTorsions(i)
              Torsion=Molecule_intra(spc)%Torsion(Step%TorsionNumber(i,j))
              a1=Molecule_Intra(spc)%TorsionAtom(Step%TorsionNumber(i,j),1)
              a4=Molecule_Intra(spc)%TorsionAtom(Step%TorsionNumber(i,j),4)
              if(atm == a1)then
                Position_p3=TrialCoordinates(spc,sys)%Positions(:,mol,a4)
              else
                Position_p3=TrialCoordinates(spc,sys)%Positions(:,mol,a1)
              end if
              vec=Position_p2-Position_p3
              vecY1=cross_product(vec,vecZ)
              vecY1=vecY1/norm2(vecY1)
              vecX1=cross_product(vecY1,vecZ)
              phi=atan2(dot_product(r3,vecY1),dot_product(r3,vecX1))
              energy=energy+TorsionInteraction(Torsion,phi)
            end do
          end do
          rosenT(TorsionTrialNumber)=exp(-beta*energy)
        end do
        if(.not. existing_molecule)then
          !** Choose Torsion Angle
          TorsionTrialNumber=ChosenTrialNumber(NumberOfTorsionTrials,rosenT)
          phi_T(TrialNumber)=phi_trial(TorsionTrialNumber)
        end if
        TorsionRosenbluthWeight=sum(rosenT)/real(NumberOfTorsionTrials,PR)
      end if
       
      !** Compute position and calculate LJ energy
      energy=0._PR
      loop2:do i=1,Step%NumberOfAtoms
        atm=Step%Atom(i)
        vec=sin(Theta(i))*cos(phi_B(i)+phi_T(TrialNumber))*vecX+ &
            sin(Theta(i))*sin(phi_B(i)+phi_T(TrialNumber))*vecY- &
            cos(Theta(i))*vecZ
        vec=Position_p1+vec*BondLength(i)
        !** Calculate Intra Pair Interactions
        do j=1,Step%NumberOfIntraPairs(i)
          IntraPair=Molecule_intra(spc)%IntraPair(Step%IntraPairNumber(i,j))
          a1=Molecule_intra(spc)%IntraPairAtom(Step%IntraPairNumber(i,j),1)
          a2=Molecule_intra(spc)%IntraPairAtom(Step%IntraPairNumber(i,j),2)
          if(atm == a2)then
            rijsq=sum((TrialCoordinates(spc,sys)%Positions(:,mol,a1)-vec)**2)
          else
            rijsq=sum((TrialCoordinates(spc,sys)%Positions(:,mol,a2)-vec)**2)
          end if
          energy=energy+IntraPairInteraction(Molecule_intra(spc)%IntraPair(i),rijsq,overlap)
          if(overlap)exit loop2
        end do

        !** Calculate Intermolecular Interactions
        call AtomSystemShortRangePairwiseInteraction(vec,vec,CurrentCoordinates(:,sys),CurrentSimulationCell(sys), &
          MainCellList(sys),spc,mol,atm,TrialInteractions,overlap)
        if(overlap)exit loop2
        call UpdateStorageInteractions(TrialInteractions)
        energy=energy+TrialInteractions%Energy%Total
      end do loop2
      !** RosenbluthWeight
      if(overlap)then
        rosen(TrialNumber)=0._PR
      else
        Success=.true.
        rosen(TrialNumber)=exp(-beta*energy)*TorsionRosenbluthWeight
        call ieee_get_flag(ieee_underflow,underflow)
        if(underflow)then
          call ieee_set_flag(ieee_underflow,.false.)
          rosen(TrialNumber)=0._PR
        end if
      end if
    end do
    overlap=(.not. Success)
    if(overlap)return

    !** Choose Trial Number
    if(existing_molecule)then
      TrialNumber=1
    else
      TrialNumber=ChosenTrialNumber(NumberOfTrialInsertions,rosen)
      do i=1,Step%NumberOfAtoms
        atm=Step%Atom(i)
        vec=sin(Theta(i))*cos(phi_B(i)+phi_T(TrialNumber))*vecX+ &
            sin(Theta(i))*sin(phi_B(i)+phi_T(TrialNumber))*vecY- &
            cos(Theta(i))*vecZ
        vec=Position_p1+vec*Bondlength(i)
        TrialCoordinates(spc,sys)%Positions(1:3,mol,atm)=vec
      end do
    end if
    RosenbluthWeight=RosenbluthWeight*sum(rosen)/real(NumberOfTrialInsertions,PR)

  end function BranchPointSampling

  !=================================================================================================
  function ChosenTrialNumber(NumberOfTrials,rosen) result(TrialNumber)
    integer, intent(in)                :: NumberOfTrials
    real(PR), dimension(:), intent(in) :: rosen
    integer                            :: TrialNumber

    real(PR) :: xrnd,wt

    xrnd=RandomNumber()*sum(rosen)
    
    wt=0._PR
    do TrialNumber=1,NumberOfTrials
      wt=wt+rosen(TrialNumber)
      if(xrnd < wt)return
    end do

  end function ChosenTrialNumber

  !=================================================================================================
  subroutine ResetCBMC_Step_Info(Step)
    type(CBMC_Step_Info), intent(out) :: Step
  end subroutine ResetCBMC_Step_Info

end module rosenbluth
