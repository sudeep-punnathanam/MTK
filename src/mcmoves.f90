module mcmoves
!!$  use, intrinsic :: ieee_exceptions
!!$  use, intrinsic :: ieee_features, only: ieee_underflow_flag
!!$  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use consts, only: PR, PI, strlen, lstrlen, MAX_NO_OF_SPECIES, MAX_NO_OF_SYSTEMS
  use variables, only: CutOffDistance, NumberOfSystems, TotalNumberOfMolecules, TotalMass, beta, Pressure, Fugacity, &
    CurrentCoordinates, CurrentSimulationCell, CurrentInteractions, TrialCoordinates, TrialSimulationCell, CurrentKfactor, CurrentStructureFactor, &
    MainCellList
  use utils, only: ReadStringFromFile, ReadString, lowercase, ErrorMessage, Matrix_x_Vector
  use atoms_and_molecules, only: GetSpeciesNumber, Molecule, NumberOfSpecies, GetCenterOfMass
  use random, only: RandomNumber, RandomRotationMatrix, UniformRandomRotationMatrix
  use storage, only: StorageInteractions, ResetStorageInteractions, UpdateStorageInteractions, operator(+), operator(-)
  use inter, only: MoleculeSystemShortRangePairwiseInteraction, MoleculeSystemLongRangePairwiseInteraction, &
    TotalShortRangePairwiseInteraction, TotalLongRangePairwiseInteraction, IntraMolecularInteraction
  use simcell, only: ApplyBoundaryCondition, SetSimulationCellBoxLengths
  use config, only: ReAllocateMemoryForCoordinates
  use rosenbluth, only: CBMC_Move, CBMC_grow, CBMC_regrow
  use ewaldsum, only: KfactorInfo, StructureFactorInfo, CalculateMoleculeStructureFactor, EwaldFourierInteraction, &
    CalculateKfactors, CalculateTotalStructureFactor, NumberOfKVectors, COULOMB_INTERACTION, &
    TotalIntraMolecularFourierInteraction, IntramolecularFourier
  use cell_list, only: UpdateCellList, InsertMoleculeToCellList, DeleteMoleculeFromCellList, CreateCellList, ResizeCellListArrays
  implicit none
  public
  save

  !=============================================================================================================================
  !                           Monte Carlo Moves
  !=============================================================================================================================
  !** Move Definition
  type :: MoveParameters
    integer                                         :: Weight=0,CumulativeWeight=0
    real(PR), dimension(5)                          :: arg=0._PR
    integer, dimension(5)                           :: iarg=0
    integer                                         :: Success=0,Attempts=0
  end type MoveParameters

  !=============================================================================================================================
  !** Particle Moves
  type(MoveParameters), dimension(MAX_NO_OF_SPECIES,MAX_NO_OF_SYSTEMS)    :: Translation
  type(MoveParameters), dimension(MAX_NO_OF_SPECIES,MAX_NO_OF_SYSTEMS)    :: Rotation
  type(MoveParameters), dimension(MAX_NO_OF_SPECIES,MAX_NO_OF_SYSTEMS)    :: Insertion
  type(MoveParameters), dimension(MAX_NO_OF_SPECIES,MAX_NO_OF_SYSTEMS)    :: Deletion
  type(MoveParameters), dimension(MAX_NO_OF_SPECIES,MAX_NO_OF_SYSTEMS)    :: CutAndRegrow

  !** Volume Move
  type(MoveParameters), dimension(MAX_NO_OF_SYSTEMS)                      :: VolumeChange
  integer, dimension(MAX_NO_OF_SPECIES,MAX_NO_OF_SYSTEMS)                 :: TotalWeights=0

  type(StructureFactorInfo), private :: TrialStructureFactor, OldStructureFactor, NewStructureFactor
  type(KfactorInfo) :: TrialKfactor

contains
  !=============================================================================================================================
  subroutine ReadMoveInformation(SimulationInput,error)
    character(len=lstrlen), dimension(:), intent(in) :: SimulationInput
    logical, intent(out) :: error

    integer :: StartLine, StartLine_2, lineno, pos
    integer :: spc,sys,spc1,spc2
    integer :: weight
    character(len=lstrlen) :: string
    character(len=strlen) :: spcname,spc1name,spc2name

    error=.false.

    StartLine=1
    sys=0
    scell:do
      call ReadStringFromFile(string,'simulation cell',SimulationInput(StartLine:),'Simulation Input',lineno)
      if(lineno /=0)then
        sys=sys+1
      else
        exit scell
      end if
      StartLine=StartLine+lineno

      !** Read molecule information
      StartLine_2=StartLine
      spc=0
      spc_loop:do
        call ReadStringFromFile(string,'species',SimulationInput(StartLine_2:),'Simulation Input',lineno,'simulation cell')
        if(lineno == 0)exit spc_loop
        pos=ReadString(string,'name',string)
        spcname=adjustl(string)
        spc=GetSpeciesNumber(spcname,error)
        if(error)return

        StartLine_2=StartLine_2+lineno
        call ReadStringFromFile(string,'shape',SimulationInput(StartLine_2:),'Simulation Input',lineno,'species')
        if(lineno /= 0)Molecule(spc)%Shape=lowercase(trim(adjustl(string)))
 
        !** Read Translation Move Parameters
        call ReadStringFromFile(string,'translation weight',SimulationInput(StartLine_2:),'Simulation Input',lineno,'species')
        if(lineno /= 0)read(string,*)Translation(spc,sys)%Weight
        call ReadStringFromFile(string,'maximum displacement',SimulationInput(StartLine_2:),'Simulation Input',lineno,'species')
        if(lineno /= 0)read(string,*)Translation(spc,sys)%arg(1)
 
        !** Read Rotation Move Parameters
        call ReadStringFromFile(string,'rotation weight',SimulationInput(StartLine_2:),'Simulation Input',lineno,'species')
        if(lineno /= 0)read(string,*)Rotation(spc,sys)%Weight
        Rotation(spc,sys)%arg(1)=PI
        call ReadStringFromFile(string,'rotation atom',SimulationInput(StartLine_2:),'Simulation Input',lineno,'species')
        if(lineno /= 0)read(string,*)Rotation(spc,sys)%iarg(1)
 
        !** Read Insertion and Deletion Move Parameters
        call ReadStringFromFile(string,'insertion weight',SimulationInput(StartLine_2:),'Simulation Input',lineno,'species')
        if(lineno /= 0)read(string,*)Insertion(spc,sys)%Weight
        if(lineno /= 0)Deletion(spc,sys)%Weight=Insertion(spc,sys)%Weight
 
        !** Read Cut and Regrow Move Parameters
        call ReadStringFromFile(string,'cut and regrow weight',SimulationInput(StartLine_2:),'Simulation Input',lineno,'species')
        if(lineno /= 0)read(string,*)CutAndRegrow(spc,sys)%Weight
      end do spc_loop

      call ReadStringFromFile(string,'volume change weight',SimulationInput(StartLine:),'Simulation Input',lineno,'simulation cell')
      if(lineno /= 0)then
        read(string,*)VolumeChange(sys)%Weight
        VolumeChange(sys)%arg(1)=CurrentSimulationCell(sys)%Volume/100._PR
      end if


      !** Normalize Move Probabilities
      do spc=1,MAX_NO_OF_SPECIES
        Translation(spc,sys)%CumulativeWeight= Translation(spc,sys)%Weight
        Rotation(spc,sys)%CumulativeWeight= Translation(spc,sys)%CumulativeWeight + Rotation(spc,sys)%Weight
        Insertion(spc,sys)%CumulativeWeight= Rotation(spc,sys)%CumulativeWeight + Insertion(spc,sys)%Weight
        Deletion(spc,sys)%CumulativeWeight= Insertion(spc,sys)%CumulativeWeight + Deletion(spc,sys)%Weight
        CutAndRegrow(spc,sys)%CumulativeWeight= Deletion(spc,sys)%CumulativeWeight + CutAndRegrow(spc,sys)%Weight

        TotalWeights(spc,sys)=CutAndRegrow(spc,sys)%CumulativeWeight
      end do
    end do scell

  end subroutine ReadMoveInformation

  !=============================================================================================================================
  !** Particle Move
  !=============================================================================================================================
  subroutine ParticleMove(simno,error)
    integer, intent(in) :: simno
    logical, intent(Out) :: error

    integer :: spc,sys
    integer :: i

integer :: mol
logical :: overlap
real(PR) :: rosen

    !** Select a species for particle move
    call SelectSystemAndSpecies(spc,sys)

    error=.false.
    i=int(RandomNumber()*real(TotalWeights(spc,sys),PR))
    !**  Translation
    if(i < Translation(spc,sys)%CumulativeWeight)then
!!$write(unit=6,fmt='(5x,a15)',advance='no')'Translation'
      call TranslationMove(spc,sys,simno,error)
!!$write(unit=6,fmt='(5x,e15.4)',advance='yes')CurrentInteractions(sys)%Energy%vdW
      return
    end if

    !**  Rotation
    if(i < Rotation(spc,sys)%CumulativeWeight)then
!!$write(unit=6,fmt='(5x,a15)',advance='no')'Rotation'
      call RotationMove(spc,sys,simno,error)
!!$write(unit=6,fmt='(5x,e15.4)',advance='yes')CurrentInteractions(sys)%Energy%vdW
      return
    end if

    !**  Insertion
    if(i < Insertion(spc,sys)%CumulativeWeight)then
!!$write(unit=6,fmt='(5x,a15)',advance='no')'Insertion'
      call InsertionMove(spc,sys,simno,error)
!!$write(unit=6,fmt='(5x,e15.4)',advance='yes')CurrentInteractions(sys)%Energy%vdW
      return
    end if

    !**  Deletion
    if(i < Deletion(spc,sys)%CumulativeWeight)then
!!$write(unit=6,fmt='(5x,a15)',advance='no')'Deletion'
      call DeletionMove(spc,sys,simno,error)
!!$write(unit=6,fmt='(5x,e15.4)',advance='yes')CurrentInteractions(sys)%Energy%vdW
      return
    end if

    !**  Cut and Regrow
    if(i < CutAndRegrow(spc,sys)%CumulativeWeight)then
!!$write(unit=6,fmt='(5x,a15)',advance='no')'CutAndRegrow'
      call CutAndRegrowMove(spc,sys,simno,error)
!!$write(unit=6,fmt='(5x,e15.4)',advance='yes')CurrentInteractions(sys)%Energy%vdW
      return
    end if

  end subroutine ParticleMove

  !=============================================================================================================================
  !**  Translation Move
  !=============================================================================================================================
  subroutine TranslationMove(spc,sys,simno,error)
    integer, intent(in)  :: spc,sys,simno
    logical, intent(Out) :: error

    integer  :: mol,atm
    integer  :: nmoles,natoms
    logical  :: overlap
    real(PR) :: drmax
    real(PR) :: arg
    real(PR) :: d1, d2, d3
    type(StorageInteractions) :: OldInteractions, NewInteractions, DeltaInteractions, EwaldInteractions
    logical :: underflow

    error=.false.
    Translation(spc,sys)%Attempts= Translation(spc,sys)%Attempts+1

    nmoles=CurrentCoordinates(spc,sys)%NumberOfMolecules
    if(nmoles == 0)return
    mol=min(int(RandomNumber()*real(nmoles,PR))+1,nmoles)
    natoms=Molecule(spc)%NumberOfAtoms
    drmax=Translation(spc,sys)%arg(1)

    !** Create Trial Configuration
    d1=(2._PR*RandomNumber()-1._PR)*drmax
    d2=(2._PR*RandomNumber()-1._PR)*drmax
    d3=(2._PR*RandomNumber()-1._PR)*drmax
    do atm=1,natoms
      TrialCoordinates(spc,sys)%Positions(1,1,atm)=CurrentCoordinates(spc,sys)%Positions(1,mol,atm)+d1
      TrialCoordinates(spc,sys)%Positions(2,1,atm)=CurrentCoordinates(spc,sys)%Positions(2,mol,atm)+d2
      TrialCoordinates(spc,sys)%Positions(3,1,atm)=CurrentCoordinates(spc,sys)%Positions(3,mol,atm)+d3
    end do
    TrialCoordinates(spc,sys)%CenterOfMass(1,1)=CurrentCoordinates(spc,sys)%CenterOfMass(1,mol)+d1
    TrialCoordinates(spc,sys)%CenterOfMass(2,1)=CurrentCoordinates(spc,sys)%CenterOfMass(2,mol)+d2
    TrialCoordinates(spc,sys)%CenterOfMass(3,1)=CurrentCoordinates(spc,sys)%CenterOfMass(3,mol)+d3
    call ApplyBoundaryCondition(CurrentSimulationCell(sys),natoms,TrialCoordinates(spc,sys)%Positions(1:3,1,1:natoms), &
      TrialCoordinates(spc,sys)%CenterOfMass(1:3,1))

    !** Calculate Change in Energy
    call MoleculeSystemShortRangePairwiseInteraction(CurrentCoordinates(spc,sys)%Positions(1:3,mol,1:natoms),CurrentCoordinates(spc,sys)%CenterOfMass(1:3,mol), &
      CurrentCoordinates(:,sys),CurrentSimulationCell(sys),MainCellList(sys),spc,mol,OldInteractions,overlap)
    if(overlap)then
      write(ErrorMessage,'(2a,i4,4x,a)')__FILE__,':',__LINE__, &
        'Something wrong! Existing configuration contains overlapping atoms'
      error=.true.
      return
    end if
    call MoleculeSystemShortRangePairwiseInteraction(TrialCoordinates(spc,sys)%Positions(1:3,1,1:natoms),TrialCoordinates(spc,sys)%CenterOfMass(1:3,1), &
      CurrentCoordinates(:,sys),CurrentSimulationCell(sys),MainCellList(sys),spc,mol,NewInteractions,overlap)
    if(overlap)return
    DeltaInteractions=NewInteractions-OldInteractions

    if(Molecule(spc)%HasPartialCharges)then
      call CalculateMoleculeStructureFactor(sys,spc,CurrentCoordinates(spc,sys)%Positions(1:3,mol,1:natoms),CurrentCoordinates(spc,sys)%CenterOfMass(1:3,mol), &
        CurrentSimulationCell(sys),OldStructureFactor)
      call CalculateMoleculeStructureFactor(sys,spc,TrialCoordinates(spc,sys)%Positions(1:3,1,1:natoms),TrialCoordinates(spc,sys)%CenterOfMass(1:3,1), &
        CurrentSimulationCell(sys),NewStructureFactor)
      call UpdateStructureFactor(sys,TrialStructureFactor,CurrentStructureFactor(sys),NewStructureFactor,OldStructureFactor)
      call EwaldFourierInteraction(sys,CurrentSimulationCell(sys)%Volume,CurrentKfactor(sys),TrialStructureFactor,EwaldInteractions)
      EwaldInteractions%Energy%CoulFourier=EwaldInteractions%Energy%CoulFourier-CurrentInteractions(sys)%Energy%CoulFourier
      EwaldInteractions%Virial%CoulFourier=EwaldInteractions%Virial%CoulFourier-CurrentInteractions(sys)%Virial%CoulFourier
      DeltaInteractions=DeltaInteractions+EwaldInteractions
    end if

    call UpdateStorageInteractions(DeltaInteractions)

    !** Accept or Reject the move
    arg=exp(-beta(simno,sys)*DeltaInteractions%Energy%Total)
!!$    call ieee_get_flag(ieee_underflow,underflow)
!!$    if(underflow)then
!!$      call ieee_set_flag(ieee_underflow,.false.)
!!$      arg=0._PR
!!$    end if
    if(RandomNumber() < arg)then ! .or. DeltaInteractions%Energy%Total < 0._PR)then
      Translation(spc,sys)%Success=Translation(spc,sys)%Success+1
      if(Molecule(spc)%HasPartialCharges)CurrentStructureFactor(sys)=TrialStructureFactor
      CurrentInteractions(sys)=CurrentInteractions(sys)+DeltaInteractions
      call UpdateStorageInteractions(CurrentInteractions(sys))

      if(MainCellList(sys)%active)call UpdateCellList(spc,mol,CurrentCoordinates(spc,sys)%Positions(:,mol,:), &
       TrialCoordinates(spc,sys)%Positions(:,1,:),CurrentSimulationCell(sys),MainCellList(sys))

      do atm=1,natoms
        CurrentCoordinates(spc,sys)%Positions(1,mol,atm)=TrialCoordinates(spc,sys)%Positions(1,1,atm)
        CurrentCoordinates(spc,sys)%Positions(2,mol,atm)=TrialCoordinates(spc,sys)%Positions(2,1,atm)
        CurrentCoordinates(spc,sys)%Positions(3,mol,atm)=TrialCoordinates(spc,sys)%Positions(3,1,atm)
      end do
      CurrentCoordinates(spc,sys)%CenterOfMass(1,mol)=TrialCoordinates(spc,sys)%CenterOfMass(1,1)
      CurrentCoordinates(spc,sys)%CenterOfMass(2,mol)=TrialCoordinates(spc,sys)%CenterOfMass(2,1)
      CurrentCoordinates(spc,sys)%CenterOfMass(3,mol)=TrialCoordinates(spc,sys)%CenterOfMass(3,1)
    end if
  end subroutine TranslationMove

  !=============================================================================================================================
  !** Rotation Move
  !=============================================================================================================================
  subroutine RotationMove(spc,sys,simno,error)
    integer, intent(in)  :: spc,sys,simno
    logical, intent(Out) :: error

    integer  :: mol,atm
    integer  :: nmoles,natoms
    logical  :: overlap
    real(PR) :: theta_max
    real(PR) :: arg
    real(PR), dimension(3) :: RotatePoint, vec
    real(PR), dimension(3,3) :: RotationMatrix
    type(StorageInteractions) :: OldInteractions, NewInteractions, DeltaInteractions, EwaldInteractions
    logical :: underflow

    error=.false.
    Rotation(spc,sys)%Attempts=Rotation(spc,sys)%Attempts+1

    nmoles=CurrentCoordinates(spc,sys)%NumberOfMolecules
    if(nmoles == 0)return
    mol=min(int(RandomNumber()*real(nmoles,PR))+1,nmoles)
    natoms=Molecule(spc)%NumberOfAtoms

    !** Create Trial Configuration
    !** Obtain a Rotation Matrix for a random rotation
    theta_max=Rotation(spc,sys)%arg(1)
    call RandomRotationMatrix(RotationMatrix,theta_max)
    !** Rotate the molecule about its Rotation atom
    if(Rotation(spc,sys)%iarg(1) == 0)then
      !** Rotate about center of mass
      RotatePoint=CurrentCoordinates(spc,sys)%CenterOfMass(:,mol)
    else
      !** Rotate about atom iarg(1)
      RotatePoint=CurrentCoordinates(spc,sys)%Positions(:,mol,Rotation(spc,sys)%iarg(1))
    end if
    !** Rotate atoms of the molecule
    do atm=1,natoms
      vec(1:3)=CurrentCoordinates(spc,sys)%Positions(1:3,mol,atm)-RotatePoint(1:3)
      vec=Matrix_x_Vector(RotationMatrix,vec)
      TrialCoordinates(spc,sys)%Positions(1:3,1,atm)=RotatePoint(1:3)+vec(1:3)
    end do
    TrialCoordinates(spc,sys)%CenterOfMass(1:3,1)=GetCenterOfMass(spc,TrialCoordinates(spc,sys)%Positions(1:3,1,1:natoms))

    call ApplyBoundaryCondition(CurrentSimulationCell(sys),natoms,TrialCoordinates(spc,sys)%Positions(1:3,1,1:natoms), &
      TrialCoordinates(spc,sys)%CenterOfMass(1:3,1))
 
    !** Calculate Change in Energy
    call MoleculeSystemShortRangePairwiseInteraction(CurrentCoordinates(spc,sys)%Positions(1:3,mol,1:natoms),CurrentCoordinates(spc,sys)%CenterOfMass(1:3,mol), &
      CurrentCoordinates(:,sys),CurrentSimulationCell(sys),MainCellList(sys),spc,mol,OldInteractions,overlap)
    if(overlap)then
      write(ErrorMessage,'(2a,i4,4x,a)')__FILE__,':',__LINE__, &
        'Something wrong! Existing configuration contains overlapping atoms'
      error=.true.
      return
    end if
    call MoleculeSystemShortRangePairwiseInteraction(TrialCoordinates(spc,sys)%Positions(1:3,1,1:natoms),TrialCoordinates(spc,sys)%CenterOfMass(1:3,1), &
      CurrentCoordinates(:,sys),CurrentSimulationCell(sys),MainCellList(sys),spc,mol,NewInteractions,overlap)
    if(overlap)return
    DeltaInteractions=NewInteractions-OldInteractions

    if(Molecule(spc)%HasPartialCharges)then
      call CalculateMoleculeStructureFactor(sys,spc,CurrentCoordinates(spc,sys)%Positions(1:3,mol,1:natoms),CurrentCoordinates(spc,sys)%CenterOfMass(1:3,mol), &
        CurrentSimulationCell(sys),OldStructureFactor)
      call CalculateMoleculeStructureFactor(sys,spc,TrialCoordinates(spc,sys)%Positions(1:3,1,1:natoms),TrialCoordinates(spc,sys)%CenterOfMass(1:3,1), &
        CurrentSimulationCell(sys),NewStructureFactor)
      call UpdateStructureFactor(sys,TrialStructureFactor,CurrentStructureFactor(sys),NewStructureFactor,OldStructureFactor)
      call EwaldFourierInteraction(sys,CurrentSimulationCell(sys)%Volume,CurrentKfactor(sys),TrialStructureFactor,EwaldInteractions)
      EwaldInteractions%Energy%CoulFourier=EwaldInteractions%Energy%CoulFourier-CurrentInteractions(sys)%Energy%CoulFourier
      EwaldInteractions%Virial%CoulFourier=EwaldInteractions%Virial%CoulFourier-CurrentInteractions(sys)%Virial%CoulFourier
      DeltaInteractions=DeltaInteractions+EwaldInteractions
    end if

    call UpdateStorageInteractions(DeltaInteractions)

    !** Accept or Reject the move
    arg=exp(-beta(simno,sys)*DeltaInteractions%Energy%Total)
!!$    call ieee_get_flag(ieee_underflow,underflow)
!!$    if(underflow)then
!!$      call ieee_set_flag(ieee_underflow,.false.)
!!$      arg=0._PR
!!$    end if
    if(RandomNumber() < arg)then ! .or. DeltaInteractions%Energy%Total < 0._PR)then
      Rotation(spc,sys)%Success=Rotation(spc,sys)%Success+1
      if(Molecule(spc)%HasPartialCharges)CurrentStructureFactor(sys)=TrialStructureFactor
      CurrentInteractions(sys)=CurrentInteractions(sys)+DeltaInteractions
      call UpdateStorageInteractions(CurrentInteractions(sys))

      if(MainCellList(sys)%active)call UpdateCellList(spc,mol,CurrentCoordinates(spc,sys)%Positions(:,mol,:), &
       TrialCoordinates(spc,sys)%Positions(:,1,:),CurrentSimulationCell(sys),MainCellList(sys))

      do atm=1,natoms
        CurrentCoordinates(spc,sys)%Positions(1,mol,atm)=TrialCoordinates(spc,sys)%Positions(1,1,atm)
        CurrentCoordinates(spc,sys)%Positions(2,mol,atm)=TrialCoordinates(spc,sys)%Positions(2,1,atm)
        CurrentCoordinates(spc,sys)%Positions(3,mol,atm)=TrialCoordinates(spc,sys)%Positions(3,1,atm)
      end do
      CurrentCoordinates(spc,sys)%CenterOfMass(1,mol)=TrialCoordinates(spc,sys)%CenterOfMass(1,1)
      CurrentCoordinates(spc,sys)%CenterOfMass(2,mol)=TrialCoordinates(spc,sys)%CenterOfMass(2,1)
      CurrentCoordinates(spc,sys)%CenterOfMass(3,mol)=TrialCoordinates(spc,sys)%CenterOfMass(3,1)
    end if
  end subroutine RotationMove

  !=============================================================================================================================
  !**  Insertion Move
  !=============================================================================================================================
  subroutine InsertionMove(spc,sys,simno,error)
    integer, intent(in)  :: spc,sys,simno
    logical, intent(Out) :: error

    integer :: SequenceNumber
    integer :: natoms, nmoles, mol
    real(PR) :: RosenbluthWeight, arg, energy, deltaU
    logical :: overlap, existing_molecule
    type(StorageInteractions) :: Interactions, LongRangeInteractions, EwaldInteractions

    error=.false.
    Insertion(spc,sys)%Attempts=Insertion(spc,sys)%Attempts+1

    nmoles=CurrentCoordinates(spc,sys)%NumberOfMolecules
    natoms=Molecule(spc)%NumberOfAtoms
    mol=nmoles+1

    !** Check and Increase Storage
    call ReAllocateMemoryForCoordinates(CurrentCoordinates(spc,sys),nmoles+1,natoms,error)
    if(error)return
    call ReAllocateMemoryForCoordinates(TrialCoordinates(spc,sys),nmoles+1,natoms,error)
    if(error)return

    !** Calculate Change in Energy and Rosenbluth Weight
    existing_molecule=.false.
    SequenceNumber=int(RandomNumber()*CBMC_Move(spc)%NumberOfGrowthSequences)+1
    RosenbluthWeight=CBMC_grow(SequenceNumber,sys,spc,mol,beta(simno,sys),existing_molecule,overlap)
    if(overlap)return
    CurrentCoordinates(spc,sys)%Positions(:,mol,:)=TrialCoordinates(spc,sys)%Positions(:,mol,:)
    CurrentCoordinates(spc,sys)%CenterOfMass(:,mol)=GetCenterOfMass(spc,TrialCoordinates(spc,sys)%Positions(:,mol,:))

    !** Calculate Short Range van der Waals and Ewald Real Space Interactions
    call MoleculeSystemShortRangePairwiseInteraction(CurrentCoordinates(spc,sys)%Positions(1:3,mol,1:natoms), &
      CurrentCoordinates(spc,sys)%CenterOfMass(1:3,mol), CurrentCoordinates(:,sys),CurrentSimulationCell(sys), &
      MainCellList(sys),spc,mol,Interactions,overlap)

    !** Calculate Intermolecular Interactions
    Interactions=Interactions+IntraMolecularInteraction(CurrentCoordinates(spc,sys)%Positions(:,mol,:),spc,mol,overlap)
    call UpdateStorageInteractions(Interactions)

    !** Calculate Long Range van der Waals and Ewald Fourier Interactions
    call MoleculeSystemLongRangePairwiseInteraction(spc,mol,CurrentCoordinates(:,sys)%NumberOfMolecules,CurrentSimulationCell(sys)%Volume,LongRangeInteractions)
    call UpdateStorageInteractions(LongRangeInteractions)
    Interactions=Interactions+LongRangeInteractions
    deltaU=LongRangeInteractions%Energy%Total

    if(Molecule(spc)%HasPartialCharges)then
      call CalculateMoleculeStructureFactor(sys,spc,CurrentCoordinates(spc,sys)%Positions(1:3,mol,1:natoms),CurrentCoordinates(spc,sys)%CenterOfMass(1:3,mol), &
        CurrentSimulationCell(sys),NewStructureFactor)
      call AddToTheStructureFactor(sys,TrialStructureFactor,CurrentStructureFactor(sys),NewStructureFactor)
      call EwaldFourierInteraction(sys,CurrentSimulationCell(sys)%Volume,CurrentKfactor(sys),TrialStructureFactor,EwaldInteractions)
      call IntramolecularFourier(spc,CurrentCoordinates(spc,sys)%Positions(1:3,mol,1:natoms),energy)
      EwaldInteractions%Energy%IntraFourier=energy
      EwaldInteractions%Energy%CoulFourier=EwaldInteractions%Energy%CoulFourier-CurrentInteractions(sys)%Energy%CoulFourier
      EwaldInteractions%Virial%CoulFourier=EwaldInteractions%Virial%CoulFourier-CurrentInteractions(sys)%Virial%CoulFourier
      Interactions=Interactions+EwaldInteractions
      deltaU=deltaU+EwaldInteractions%Energy%CoulFourier-EwaldInteractions%Energy%IntraFourier
    end if

    RosenbluthWeight=RosenbluthWeight*exp(-beta(simno,sys)*deltaU)
    call UpdateStorageInteractions(Interactions)

    !** Accept or Reject the move
    arg=beta(simno,sys)*Fugacity(spc,simno,sys)*CurrentSimulationCell(sys)%Volume/real(nmoles+1,PR)
    arg=arg*RosenbluthWeight/CBMC_Move(spc)%GrowthSequence(SequenceNumber)%IdealRosenbluthWeight
    if(RandomNumber() < arg)then
      Insertion(spc,sys)%Success=Insertion(spc,sys)%Success+1
      if(Molecule(spc)%HasPartialCharges)CurrentStructureFactor(sys)=TrialStructureFactor
      CurrentInteractions(sys)=CurrentInteractions(sys)+Interactions
      call UpdateStorageInteractions(CurrentInteractions(sys))

      if(MainCellList(sys)%active)then
        call InsertMoleculeToCellList(spc,mol,CurrentCoordinates(spc,sys)%Positions(:,mol,:),CurrentSimulationCell(sys),MainCellList(sys),error)
        if(error)return
      end if

      CurrentCoordinates(spc,sys)%NumberOfMolecules=CurrentCoordinates(spc,sys)%NumberOfMolecules+1
      TotalNumberOfMolecules(sys)=TotalNumberOfMolecules(sys)+1
      TotalMass(sys)=TotalMass(sys)+Molecule(spc)%MolecularWeight
    end if

  end subroutine InsertionMove

  !=============================================================================================================================
  !**  Deletion Move
  !=============================================================================================================================
  subroutine DeletionMove(spc,sys,simno,error)
    integer, intent(in)  :: spc,sys,simno
    logical, intent(Out) :: error

    integer :: SequenceNumber
    integer :: natoms, nmoles, mol
    logical :: overlap, existing_molecule
    real(PR) :: RosenbluthWeight, arg, energy, deltaU
    type(StorageInteractions) :: Interactions, LongRangeInteractions, EwaldInteractions

    error=.false.
    Deletion(spc,sys)%Attempts=Deletion(spc,sys)%Attempts+1

    nmoles=CurrentCoordinates(spc,sys)%NumberOfMolecules
    natoms=Molecule(spc)%NumberOfAtoms
    if(nmoles == 0)return
    !** Select a molecule at random
    mol=min(int(RandomNumber()*real(nmoles,PR))+1,nmoles)

    !** Calculate Change in Energy and Rosenbluth Weight
    existing_molecule=.true.
    SequenceNumber=int(RandomNumber()*CBMC_Move(spc)%NumberOfGrowthSequences)+1
    RosenbluthWeight=CBMC_grow(SequenceNumber,sys,spc,mol,beta(simno,sys),existing_molecule,overlap)
    if(overlap)then
      write(ErrorMessage,'(2a,i4,4x,a)')__FILE__,':',__LINE__, &
        'Something wrong! Existing configuration contains overlapping atoms'
      error=.true.
      return
    end if

    !** Calculate Short Range van der Waals and Ewald Real Space Interactions
    call MoleculeSystemShortRangePairwiseInteraction(CurrentCoordinates(spc,sys)%Positions(1:3,mol,1:natoms), &
      CurrentCoordinates(spc,sys)%CenterOfMass(1:3,mol), CurrentCoordinates(:,sys),CurrentSimulationCell(sys), &
      MainCellList(sys),spc,mol,Interactions,overlap)
    call UpdateStorageInteractions(Interactions)

    !** Calculate Intermolecular Interactions
    Interactions=Interactions+IntraMolecularInteraction(CurrentCoordinates(spc,sys)%Positions(:,mol,:),spc,mol,overlap)
    call UpdateStorageInteractions(Interactions)

    !** Calculate Long Range van der Waals
    call MoleculeSystemLongRangePairwiseInteraction(spc,nmoles,CurrentCoordinates(:,sys)%NumberOfMolecules,CurrentSimulationCell(sys)%Volume,LongRangeInteractions)
    if(error)return
    call UpdateStorageInteractions(LongRangeInteractions)
    Interactions=Interactions+LongRangeInteractions
    deltaU=LongRangeInteractions%Energy%Total

    !** Calculate Ewald Fourier Interactions
    if(Molecule(spc)%HasPartialCharges)then
      call CalculateMoleculeStructureFactor(sys,spc,CurrentCoordinates(spc,sys)%Positions(1:3,mol,1:natoms),CurrentCoordinates(spc,sys)%CenterOfMass(1:3,mol), &
        CurrentSimulationCell(sys),OldStructureFactor)
      call RemoveFromTheStructureFactor(sys,TrialStructureFactor,CurrentStructureFactor(sys),OldStructureFactor)
      call EwaldFourierInteraction(sys,CurrentSimulationCell(sys)%Volume,CurrentKfactor(sys),TrialStructureFactor,EwaldInteractions)
      call IntramolecularFourier(spc,CurrentCoordinates(spc,sys)%Positions(1:3,mol,1:natoms),energy)
      EwaldInteractions%Energy%IntraFourier=energy
      EwaldInteractions%Energy%CoulFourier=CurrentInteractions(sys)%Energy%CoulFourier-EwaldInteractions%Energy%CoulFourier
      EwaldInteractions%Virial%CoulFourier=CurrentInteractions(sys)%Virial%CoulFourier-EwaldInteractions%Virial%CoulFourier
      Interactions=Interactions+EwaldInteractions
      deltaU=deltaU+EwaldInteractions%Energy%CoulFourier-EwaldInteractions%Energy%IntraFourier
    end if

    RosenbluthWeight=RosenbluthWeight*exp(-beta(simno,sys)*deltaU)
    call UpdateStorageInteractions(Interactions)

    !** Accept or Reject the move
    arg=real(nmoles,PR)/(beta(simno,sys)*Fugacity(spc,simno,sys)*CurrentSimulationCell(sys)%Volume)
    arg=arg*CBMC_Move(spc)%GrowthSequence(SequenceNumber)%IdealRosenbluthWeight/RosenbluthWeight
    if(RandomNumber() < arg)then
      Deletion(spc,sys)%Success=Deletion(spc,sys)%Success+1
      if(Molecule(spc)%HasPartialCharges)CurrentStructureFactor(sys)=TrialStructureFactor
      CurrentInteractions(sys)=CurrentInteractions(sys)-Interactions
      call UpdateStorageInteractions(CurrentInteractions(sys))

      if(MainCellList(sys)%active)call DeleteMoleculeFromCellList(spc,mol,nmoles,CurrentCoordinates(spc,sys)%Positions(:,mol,:), &
        CurrentCoordinates(spc,sys)%Positions(:,nmoles,:),CurrentSimulationCell(sys),MainCellList(sys))

      CurrentCoordinates(spc,sys)%NumberOfMolecules=CurrentCoordinates(spc,sys)%NumberOfMolecules-1
      TotalNumberOfMolecules(sys)=TotalNumberOfMolecules(sys)-1
      TotalMass(sys)=TotalMass(sys)-Molecule(spc)%MolecularWeight
      CurrentCoordinates(spc,sys)%Positions(1:3,mol,1:natoms)=CurrentCoordinates(spc,sys)%Positions(1:3,nmoles,1:natoms)
      CurrentCoordinates(spc,sys)%CenterOfMass(1:3,mol)=CurrentCoordinates(spc,sys)%CenterOfMass(1:3,nmoles)
    end if

  end subroutine DeletionMove

  !=============================================================================================================================
  !**  Cut and Regrow Move
  !=============================================================================================================================
  subroutine CutAndRegrowMove(spc,sys,simno,error)
    integer, intent(in)  :: spc,sys,simno
    logical, intent(Out) :: error

    integer :: SequenceNumber
    integer :: natoms, nmoles, mol
    real(PR) :: RosenbluthWeightRatio, arg, old_energy, new_energy
    logical :: overlap
    type(StorageInteractions) :: OldInteractions, NewInteractions, DeltaInteractions, EwaldInteractions
    real(PR) :: deltaU

    error=.false.
    CutAndRegrow(spc,sys)%Attempts=CutAndRegrow(spc,sys)%Attempts+1

    nmoles=CurrentCoordinates(spc,sys)%NumberOfMolecules
    if(nmoles == 0)return
    mol=min(int(RandomNumber()*real(nmoles,PR))+1,nmoles)
    natoms=Molecule(spc)%NumberOfAtoms

    !** Calculate Change in Energy and Rosenbluth Weight
    SequenceNumber=int(RandomNumber()*CBMC_Move(spc)%NumberOfGrowthSequences)+1
    RosenbluthWeightRatio=CBMC_regrow(SequenceNumber,sys,spc,mol,beta(simno,sys),overlap,error)
    if(error)return
    if(overlap)return
    TrialCoordinates(spc,sys)%CenterOfMass(:,mol)=GetCenterOfMass(spc,TrialCoordinates(spc,sys)%Positions(:,mol,:))

    !** Old Interactions
    !** Calculate Short Range van der Waals and Ewald Real Space Interactions
    call MoleculeSystemShortRangePairwiseInteraction(CurrentCoordinates(spc,sys)%Positions(1:3,mol,1:natoms),CurrentCoordinates(spc,sys)%CenterOfMass(1:3,mol), &
      CurrentCoordinates(:,sys),CurrentSimulationCell(sys),MainCellList(sys),spc,mol,OldInteractions,overlap)
    !** Calculate Intermolecular Interactions
    OldInteractions=OldInteractions+IntraMolecularInteraction(CurrentCoordinates(spc,sys)%Positions(:,mol,:),spc,mol,overlap)
    call UpdateStorageInteractions(OldInteractions)

    !** New Interactions
    !** Calculate Short Range van der Waals and Ewald Real Space Interactions
    call MoleculeSystemShortRangePairwiseInteraction(TrialCoordinates(spc,sys)%Positions(1:3,mol,1:natoms),TrialCoordinates(spc,sys)%CenterOfMass(1:3,mol), &
      CurrentCoordinates(:,sys),CurrentSimulationCell(sys),MainCellList(sys),spc,mol,NewInteractions,overlap)
    !** Calculate Intermolecular Interactions
    NewInteractions=NewInteractions+IntraMolecularInteraction(TrialCoordinates(spc,sys)%Positions(:,mol,:),spc,mol,overlap)
    call UpdateStorageInteractions(NewInteractions)

    !** Interaction Difference
    DeltaInteractions=NewInteractions-OldInteractions
    deltaU=0._PR

    !** Calculate Ewald Fourier Interactions
    if(Molecule(spc)%HasPartialCharges)then
      call CalculateMoleculeStructureFactor(sys,spc,CurrentCoordinates(spc,sys)%Positions(1:3,mol,1:natoms),CurrentCoordinates(spc,sys)%CenterOfMass(1:3,mol), &
        CurrentSimulationCell(sys),OldStructureFactor)
      call CalculateMoleculeStructureFactor(sys,spc,TrialCoordinates(spc,sys)%Positions(1:3,1,1:natoms),TrialCoordinates(spc,sys)%CenterOfMass(1:3,1), &
        CurrentSimulationCell(sys),NewStructureFactor)
      call UpdateStructureFactor(sys,TrialStructureFactor,CurrentStructureFactor(sys),NewStructureFactor,OldStructureFactor)
      call EwaldFourierInteraction(sys,CurrentSimulationCell(sys)%Volume,CurrentKfactor(sys),TrialStructureFactor,EwaldInteractions)
      call IntramolecularFourier(spc,CurrentCoordinates(spc,sys)%Positions(1:3,mol,1:natoms),old_energy)
      call IntramolecularFourier(spc,TrialCoordinates(spc,sys)%Positions(1:3,mol,1:natoms),new_energy)
      EwaldInteractions%Energy%IntraFourier=new_energy-old_energy
      EwaldInteractions%Energy%CoulFourier=EwaldInteractions%Energy%CoulFourier-CurrentInteractions(sys)%Energy%CoulFourier
      EwaldInteractions%Virial%CoulFourier=EwaldInteractions%Virial%CoulFourier-CurrentInteractions(sys)%Virial%CoulFourier
      DeltaInteractions=DeltaInteractions+EwaldInteractions
      deltaU=EwaldInteractions%Energy%CoulFourier-EwaldInteractions%Energy%IntraFourier
    end if
    call UpdateStorageInteractions(DeltaInteractions)

    !** Accept or Reject the move
    arg=RosenbluthWeightRatio*exp(-beta(simno,sys)*deltaU)
    if(RandomNumber() < arg)then
      CutAndRegrow(spc,sys)%Success=CutAndRegrow(spc,sys)%Success+1
      if(Molecule(spc)%HasPartialCharges)CurrentStructureFactor(sys)=TrialStructureFactor
      CurrentInteractions(sys)=CurrentInteractions(sys)+DeltaInteractions
      call UpdateStorageInteractions(CurrentInteractions(sys))

      if(MainCellList(sys)%active)then
        call InsertMoleculeToCellList(spc,mol,CurrentCoordinates(spc,sys)%Positions(:,mol,:),CurrentSimulationCell(sys),MainCellList(sys),error)
        if(error)return
      end if

      CurrentCoordinates(spc,sys)%Positions(:,mol,:)=TrialCoordinates(spc,sys)%Positions(:,mol,:)
      CurrentCoordinates(spc,sys)%CenterOfMass(:,mol)=TrialCoordinates(spc,sys)%CenterOfMass(:,mol)
    end if

  end subroutine CutAndRegrowMove

  !=============================================================================================================================
  !** Volume Change Move
  !=============================================================================================================================
  subroutine VolumeChangeMove(simno,error)
    integer, intent(in)  :: simno
    logical, intent(Out) :: error

    integer :: sys,spc,mol,atm
    real(PR):: ratio, ratio3
    logical :: overlap
    integer ::N,nmoles,natoms
    real(PR):: arg, energy
    real(PR), dimension(3) :: vec,com, boxlength
    type(StorageInteractions) :: TrialInteractions,Interactions
    logical :: underflow

    error=.false.

    !** Select a system
    if(NumberOfSystems > 1)then
      call SelectSystem(sys)
    else
      sys=1
    end if
    VolumeChange(sys)%Attempts=VolumeChange(sys)%Attempts+1

    !--------------------------------------------------------------------------------------------------------------------------
    !** Create Trial Configuration
    !--------------------------------------------------------------------------------------------------------------------------
    !** Change box dimensions
    ratio=(1._PR+((2._PR*RandomNumber()-1._PR)*VolumeChange(sys)%arg(1))/CurrentSimulationCell(sys)%Volume)
    ratio3=ratio**(1._PR/3._PR)
    TrialSimulationCell(sys)=CurrentSimulationCell(sys)
    boxlength=CurrentSimulationCell(sys)%BoxLength*ratio3

    call SetSimulationCellBoxLengths(TrialSimulationCell(sys),boxlength)

    if(any(0.5_PR*TrialSimulationCell(sys)%BoxWidth < CutOffDistance))then
      error=.true.
      write(ErrorMessage,'(2a,i5,4x,a,i2,a)')__FILE__,':',__LINE__, &
        'CutOffDistance exceeds half the boxwidth for system ',sys,'. Increase system size or decrease cut-off distance'
      return
    end if

    !** Generate New Coordinates
    do spc=1,NumberOfSpecies
      nmoles=CurrentCoordinates(spc,sys)%NumberOfMolecules
      natoms=Molecule(spc)%NumberOfAtoms
      call ReAllocateMemoryForCoordinates(TrialCoordinates(spc,sys),nmoles,natoms,error)
      if(error)return
      TrialCoordinates(spc,sys)%NumberOfMolecules=nmoles
      do mol=1,nmoles
        !** Get the center of mass position of the reference coordinates
        com=CurrentCoordinates(spc,sys)%CenterOfMass(:,mol)
        TrialCoordinates(spc,sys)%CenterOfMass(:,mol)=CurrentCoordinates(spc,sys)%CenterOfMass(:,mol)*ratio3
        do atm=1,natoms
          vec=CurrentCoordinates(spc,sys)%Positions(1:3,mol,atm)-com
          TrialCoordinates(spc,sys)%Positions(1:3,mol,atm)=TrialCoordinates(spc,sys)%CenterOfMass(1:3,mol)+vec
        end do
      end do
    end do

    !--------------------------------------------------------------------------------------------------------------------------
    !** Calculate Change in Energy
    !--------------------------------------------------------------------------------------------------------------------------
    if(MainCellList(sys)%active)then
      call ResizeCellListArrays(CurrentCoordinates(:,sys)%NumberOfMolecules,TrialSimulationCell(sys),MainCellList(sys),error)
      if(error)return
      call CreateCellList(TrialCoordinates(:,sys),TrialSimulationCell(sys),MainCellList(sys))
    end if
    call TotalShortRangePairwiseInteraction(TrialCoordinates(:,sys),TrialSimulationCell(sys),MainCellList(sys),Interactions,overlap,error)
    if(overlap)then
      if(MainCellList(sys)%active)then
        call ResizeCellListArrays(CurrentCoordinates(:,sys)%NumberOfMolecules,CurrentSimulationCell(sys),MainCellList(sys),error)
        if(error)return
        call CreateCellList(CurrentCoordinates(:,sys),CurrentSimulationCell(sys),MainCellList(sys))
      end if
      return
    end if
    if(error)return
    TrialInteractions=Interactions
    call TotalLongRangePairwiseInteraction(TrialCoordinates(:,sys)%NumberOfMolecules,TrialSimulationCell(sys)%Volume,Interactions)
    TrialInteractions=TrialInteractions+Interactions

    if(COULOMB_INTERACTION)then
      call CalculateKfactors(sys,TrialSimulationCell(sys),TrialKfactor)
      call CalculateTotalStructureFactor(sys,TrialCoordinates(:,sys),TrialSimulationCell(sys),TrialStructureFactor,error)
      if(error)return
      call EwaldFourierInteraction(sys,TrialSimulationCell(sys)%Volume,TrialKfactor,TrialStructureFactor,Interactions)
      TrialInteractions=TrialInteractions+Interactions
      call TotalIntraMolecularFourierInteraction(sys,TrialCoordinates(:,sys),energy)
      TrialInteractions%Energy%IntraFourier=energy
    end if

    call UpdateStorageInteractions(TrialInteractions)
    !--------------------------------------------------------------------------------------------------------------------------
    !** Accept or Reject the move
    !--------------------------------------------------------------------------------------------------------------------------
    N=TotalNumberOfMolecules(sys)
    arg=(TrialInteractions%Energy%Total-CurrentInteractions(sys)%Energy%Total)+Pressure(simno,sys)*(TrialSimulationCell(sys)%Volume-CurrentSimulationCell(sys)%Volume)
    arg=ratio**real(N-1,PR)*exp(-beta(simno,sys)*arg)
!!$    call ieee_get_flag(ieee_underflow,underflow)
!!$    if(underflow)then
!!$      call ieee_set_flag(ieee_underflow,.false.)
!!$      arg=0._PR
!!$    end if
    if(RandomNumber() < arg)then
      VolumeChange(sys)%Success=VolumeChange(sys)%Success+1
      do spc=1,NumberOfSpecies
        CurrentCoordinates(spc,sys)=TrialCoordinates(spc,sys)
      end do
      CurrentSimulationCell(sys)=TrialSimulationCell(sys)
      if(COULOMB_INTERACTION)then
        CurrentKfactor(sys)=TrialKfactor
        CurrentStructureFactor(sys)=TrialStructureFactor
      end if
      CurrentInteractions(sys)=TrialInteractions
      call UpdateStorageInteractions(CurrentInteractions(sys))
    else
      if(MainCellList(sys)%active)then
        call ResizeCellListArrays(CurrentCoordinates(:,sys)%NumberOfMolecules,CurrentSimulationCell(sys),MainCellList(sys),error)
        if(error)return
        call CreateCellList(CurrentCoordinates(:,sys),CurrentSimulationCell(sys),MainCellList(sys))
      end if
    end if
    
  end subroutine VolumeChangeMove

  !=============================================================================================================================
  subroutine UpdateMoveParameters

    integer :: spc
    real(PR) :: ratio, MinBoxWidth
    integer::sys

    do sys=1,NumberOfSystems

      do spc=1,NumberOfSpecies
        !** Translation
        if(Translation(spc,sys)%Weight /= 0 .and. Translation(spc,sys)%Attempts /= 0)then
          ratio=real(Translation(spc,sys)%Success,PR)/real(Translation(spc,sys)%Attempts,PR)
          if(ratio < 0.5_PR)then
            Translation(spc,sys)%arg(1)=0.95_PR*Translation(spc,sys)%arg(1)
          else
            Translation(spc,sys)%arg(1)=1.05_PR*Translation(spc,sys)%arg(1)
            MinBoxWidth=minval(CurrentSimulationCell(sys)%Boxwidth)
            Translation(spc,sys)%arg(1)=min(Translation(spc,sys)%arg(1),0.5_PR*MinBoxWidth)
          end if
          Translation(spc,sys)%Success=0
          Translation(spc,sys)%Attempts=0
        end if
   
        !** Rotation
        if(Rotation(spc,sys)%Weight /= 0 .and. Rotation(spc,sys)%Attempts /= 0)then
          ratio=real(Rotation(spc,sys)%Success,PR)/real(Rotation(spc,sys)%Attempts,PR)
          if(ratio < 0.5_PR)then
            Rotation(spc,sys)%arg(1)=0.95_PR*Rotation(spc,sys)%arg(1)
          else
            Rotation(spc,sys)%arg(1)=1.05_PR*Rotation(spc,sys)%arg(1)
            Rotation(spc,sys)%arg(1)=min(Rotation(spc,sys)%arg(1),PI)
          end if
          Rotation(spc,sys)%Success=0
          Rotation(spc,sys)%Attempts=0
        end if
      end do
   
      if(VolumeChange(sys)%Weight /= 0 .and. VolumeChange(sys)%Attempts /= 0)then
        ratio=real(VolumeChange(sys)%Success,PR)/real(VolumeChange(sys)%Attempts,PR)
        if(ratio < 0.5_PR)then
          VolumeChange(sys)%arg(1)=0.95_PR*VolumeChange(sys)%arg(1)
        else
          VolumeChange(sys)%arg(1)=1.05_PR*VolumeChange(sys)%arg(1)
        end if
        VolumeChange(sys)%Success=0
        VolumeChange(sys)%Attempts=0
      end if
    end do
  
  end subroutine UpdateMoveParameters

  !=============================================================================================================================
  subroutine ResetMoveStatistics
    Translation%Success=0
    Translation%Attempts=0
    Rotation%Success=0
    Rotation%Attempts=0
    Insertion%Success=0
    Insertion%Attempts=0
    Deletion%Success=0
    Deletion%Attempts=0
    CutAndRegrow%Success=0
    CutAndRegrow%Attempts=0
    VolumeChange%Success=0
    VolumeChange%Attempts=0
  end subroutine ResetMoveStatistics

  !=============================================================================================================================
  subroutine SelectSystemAndSpecies(spc,sys)
    integer, intent(out) :: spc,sys

    if(NumberOfSystems > 1)then
      sys=min(int(RandomNumber()*NumberOfSystems)+1,NumberOfSystems)
    else
      sys=1
    end if

    do
      spc=min(int(RandomNumber()*NumberOfSpecies)+1,NumberOfSpecies)
      if(TotalWeights(spc,sys) /= 0)exit
    end do

  end subroutine SelectSystemAndSpecies

  !=================================================================================================
  subroutine SelectSystem(sys)
    integer, intent(out) :: sys

    do
      sys=min(int(RandomNumber()*NumberOfSystems)+1,NumberOfSystems)
      if(VolumeChange(sys)%Weight /= 0)exit
    end do

  end subroutine SelectSystem

  !=============================================================================================================================
  subroutine DisplayMoveStatistics(sys,unitno)
    integer, intent(In)   :: sys,unitno

    character(len=strlen) :: ratio,string1,string2
    integer               :: spc,spc2

    write(unitno,'(a)')'------------------------------------------------------------------------------'
    write(unitno,'(25x,a,i2)')'SYSTEM : ',sys
    write(unitno,'(a)')'------------------------------------------------------------------------------'
    do spc=1,NumberOfSpecies
      write(unitno,'(2a)')'Move Statitics for ',trim(Molecule(spc)%Name)
      if(Translation(spc,sys)%Weight /= 0)then
        write(ratio,'(f10.4)')real(Translation(spc,sys)%Success,PR)/real(Translation(spc,sys)%Attempts,PR)
        ratio=adjustl(ratio)
        write(string1,'(i20)')Translation(spc,sys)%Success
        string1=adjustl(string1)
        write(string2,'(i20)')Translation(spc,sys)%Attempts
        string2=adjustl(string2)
        write(unitno,'(a,t31,7a)')'Translation Move','Success Ratio : ',trim(ratio),' (',trim(string1),' out of ',trim(string2),')'
      end if
      if(Rotation(spc,sys)%Weight /= 0)then
        write(ratio,'(f10.4)')real(Rotation(spc,sys)%Success,PR)/real(Rotation(spc,sys)%Attempts,PR)
        ratio=adjustl(ratio)
        write(string1,'(i20)')Rotation(spc,sys)%Success
        string1=adjustl(string1)
        write(string2,'(i20)')Rotation(spc,sys)%Attempts
        string2=adjustl(string2)
        write(unitno,'(a,t31,7a)')'Rotation Move','Success Ratio : ',trim(ratio),' (',trim(string1),' out of ',trim(string2),')'
      end if
      if(Insertion(spc,sys)%Weight /= 0)then
        write(ratio,'(f10.4)')real(Insertion(spc,sys)%Success,PR)/real(Insertion(spc,sys)%Attempts,PR)
        ratio=adjustl(ratio)
        write(string1,'(i20)')Insertion(spc,sys)%Success
        string1=adjustl(string1)
        write(string2,'(i20)')Insertion(spc,sys)%Attempts
        string2=adjustl(string2)
        write(unitno,'(a,t31,7a)')'Insertion Move','Success Ratio : ',trim(ratio),' (',trim(string1),' out of ',trim(string2),')'
      end if
      if(Deletion(spc,sys)%Weight /= 0)then
        write(ratio,'(f10.4)')real(Deletion(spc,sys)%Success,PR)/real(Deletion(spc,sys)%Attempts,PR)
        ratio=adjustl(ratio)
        write(string1,'(i20)')Deletion(spc,sys)%Success
        string1=adjustl(string1)
        write(string2,'(i20)')Deletion(spc,sys)%Attempts
        string2=adjustl(string2)
        write(unitno,'(a,t31,7a)')'Deletion Move','Success Ratio : ',trim(ratio),' (',trim(string1),' out of ',trim(string2),')'
      end if
      if(CutAndRegrow(spc,sys)%Weight /= 0)then
        write(ratio,'(f10.4)')real(CutAndRegrow(spc,sys)%Success,PR)/real(CutAndRegrow(spc,sys)%Attempts,PR)
        ratio=adjustl(ratio)
        write(string1,'(i20)')CutAndRegrow(spc,sys)%Success
        string1=adjustl(string1)
        write(string2,'(i20)')CutAndRegrow(spc,sys)%Attempts
        string2=adjustl(string2)
        write(unitno,'(a,t31,7a)')'CutAndRegrow Move','Success Ratio : ',trim(ratio),' (',trim(string1),' out of ',trim(string2),')'
      end if
    end do
    write(unitno,'(a)')'------------------------------------------------------------------------------'
    if(VolumeChange(sys)%Weight /= 0)then
      write(ratio,'(f10.4)')real(VolumeChange(sys)%Success,PR)/real(VolumeChange(sys)%Attempts,PR)
      ratio=adjustl(ratio)
      write(string1,'(i20)')VolumeChange(sys)%Success
      string1=adjustl(string1)
      write(string2,'(i20)')VolumeChange(sys)%Attempts
      string2=adjustl(string2)
      write(unitno,'(a,t31,7a)')'Volume Change Move','Success Ratio : ',trim(ratio),' (',trim(string1),' out of ',trim(string2),')'
      write(unitno,'(a)')'------------------------------------------------------------------------------'
    end if

  end subroutine DisplayMoveStatistics

  !=============================================================================================================================
  subroutine UpdateStructureFactor(sys,trial,current,new,old)
    integer, intent(in) :: sys
    type(StructureFactorInfo), intent(in) :: current,new,old
    type(StructureFactorInfo), intent(inout) :: trial

    integer :: nk
    
    do nk=1,NumberOfKVectors(sys)
      trial%RhoK(nk)=current%RhoK(nk)+new%RhoK(nk)-old%RhoK(nk)
      trial%RhokV(nk)=current%RhokV(nk)+new%RhokV(nk)-old%RhokV(nk)
    end do
    trial%Dipole(1)=current%Dipole(1)+new%Dipole(1)-old%Dipole(1)
    trial%Dipole(2)=current%Dipole(2)+new%Dipole(2)-old%Dipole(2)
    trial%Dipole(3)=current%Dipole(3)+new%Dipole(3)-old%Dipole(3)
    trial%DipoleV(1)=current%DipoleV(1)+new%DipoleV(1)-old%DipoleV(1)
    trial%DipoleV(2)=current%DipoleV(2)+new%DipoleV(2)-old%DipoleV(2)
    trial%DipoleV(3)=current%DipoleV(3)+new%DipoleV(3)-old%DipoleV(3)

  end subroutine UpdateStructureFactor

  !=============================================================================================================================
  subroutine AddToTheStructureFactor(sys,trial,current,new)
    integer, intent(in) :: sys
    type(StructureFactorInfo), intent(in) :: current,new
    type(StructureFactorInfo), intent(inout) :: trial

    integer :: nk
    
    do nk=1,NumberOfKVectors(sys)
      trial%RhoK(nk)=current%RhoK(nk)+new%RhoK(nk)
      trial%RhokV(nk)=current%RhokV(nk)+new%RhokV(nk)
    end do
    trial%Dipole(1)=current%Dipole(1)+new%Dipole(1)
    trial%Dipole(2)=current%Dipole(2)+new%Dipole(2)
    trial%Dipole(3)=current%Dipole(3)+new%Dipole(3)
    trial%DipoleV(1)=current%DipoleV(1)+new%DipoleV(1)
    trial%DipoleV(2)=current%DipoleV(2)+new%DipoleV(2)
    trial%DipoleV(3)=current%DipoleV(3)+new%DipoleV(3)

  end subroutine AddToTheStructureFactor

  !=============================================================================================================================
  subroutine RemoveFromTheStructureFactor(sys,trial,current,old)
    integer, intent(in) :: sys
    type(StructureFactorInfo), intent(in) :: current,old
    type(StructureFactorInfo), intent(inout) :: trial

    integer :: nk
    
    do nk=1,NumberOfKVectors(sys)
      trial%RhoK(nk)=current%RhoK(nk)-old%RhoK(nk)
      trial%RhokV(nk)=current%RhokV(nk)-old%RhokV(nk)
    end do
    trial%Dipole(1)=current%Dipole(1)-old%Dipole(1)
    trial%Dipole(2)=current%Dipole(2)-old%Dipole(2)
    trial%Dipole(3)=current%Dipole(3)-old%Dipole(3)
    trial%DipoleV(1)=current%DipoleV(1)-old%DipoleV(1)
    trial%DipoleV(2)=current%DipoleV(2)-old%DipoleV(2)
    trial%DipoleV(3)=current%DipoleV(3)-old%DipoleV(3)

  end subroutine RemoveFromTheStructureFactor

end module mcmoves

