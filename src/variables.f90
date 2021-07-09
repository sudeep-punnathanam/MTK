module variables
  use consts, only: PR, PI, strlen, lstrlen, xlstrlen, MAX_NO_OF_SIMULATIONS, MAX_NO_OF_SYSTEMS, MAX_NO_OF_SPECIES, NO_OF_SUBSETS, &
    K_B, UNIT_PRESSURE, UNIT_DENSITY
  use utils, only: ErrorMessage, ReadStringFromFile, ReadString, GetSubsetNumber, swap, uppercase, FindStringInFile, GetFileUnit, lowercase
  use atoms_and_molecules, only: GetSpeciesNumber, Molecule
  use config, only: SpeciesCoordinates, AllocateMemoryForCoordinates
  use simcell, only: SimcellInfo, SetSimulationCellBoxLengthsAndBoxAngles
  use storage, only: StorageInteractions
  use ewaldsum, only: KfactorInfo, StructureFactorInfo
  use cell_list, only: CellListInfo, InitializeCellListArrays, CreateCellList

  implicit none
  public
  save

  !=============================================================================================================================
  !                           Simulation Variables
  !=============================================================================================================================
  logical  :: REDUCED_UNITS
  real(PR) :: CutOffDistance, CutOffDistanceSq
  integer  :: NumberOfSimulations=0
  integer  :: NumberOfEquilibrationCycles=0,NumberOfProductionCycles=0
  integer  :: RestartFrequency=1000,ConfigsFrequency=10
  integer  :: DisplayFrequency=1000, UpdateFrequency=20


  !=============================================================================================================================
  !                           Simulation Cell Variables
  !=============================================================================================================================
  integer  :: NumberOfSystems
  type(SimcellInfo), dimension(MAX_NO_OF_SYSTEMS), target :: CurrentSimulationCell, TrialSimulationCell, OldSimulationCell

  !=============================================================================================================================
  !                           Coordinates
  !=============================================================================================================================
  type(SpeciesCoordinates), dimension(MAX_NO_OF_SPECIES,MAX_NO_OF_SYSTEMS), target :: CurrentCoordinates, TrialCoordinates, OldCoordinates
  integer, dimension(MAX_NO_OF_SYSTEMS) :: TotalNumberOfMolecules=0
  real(PR), dimension(MAX_NO_OF_SYSTEMS) :: TotalMass=0._PR
!!$  integer, dimension(MAX_NO_OF_SPECIES,MAX_NO_OF_SYSTEMS) :: InitialNumberOfMolecules
  character(len=strlen), dimension(MAX_NO_OF_SPECIES,MAX_NO_OF_SYSTEMS) :: initmethod

  !=============================================================================================================================
  !                           State Variables
  !=============================================================================================================================
  real(PR), dimension(MAX_NO_OF_SIMULATIONS,MAX_NO_OF_SYSTEMS) :: Temperature=1._PR,beta=1._PR
  real(PR), dimension(MAX_NO_OF_SIMULATIONS,MAX_NO_OF_SYSTEMS) :: Pressure=1._PR
  real(PR), dimension(MAX_NO_OF_SPECIES,MAX_NO_OF_SIMULATIONS,MAX_NO_OF_SYSTEMS) :: Fugacity=1._PR
  real(PR), dimension(NO_OF_SUBSETS,MAX_NO_OF_SIMULATIONS,MAX_NO_OF_SYSTEMS) :: FugacityRatio=1._PR

  !=============================================================================================================================
  !                           Storage Interactions
  !=============================================================================================================================
  type(StorageInteractions), dimension(MAX_NO_OF_SYSTEMS) :: CurrentInteractions, OldInteractions
  type(KfactorInfo), dimension(MAX_NO_OF_SYSTEMS) :: CurrentKfactor, OldKfactor
  type(StructureFactorInfo), dimension(MAX_NO_OF_SYSTEMS) :: CurrentStructureFactor, OldStructureFactor

  !=============================================================================================================================
  type(CellListInfo), dimension(MAX_NO_OF_SYSTEMS) :: MainCellList

contains
  !=============================================================================================================================
  !**       Reduced units
  !=============================================================================================================================
  subroutine SetReducedUnits(SimulationInput)
    character(len=lstrlen), dimension(:), intent(in) :: SimulationInput

    character(len=lstrlen) :: string
    integer :: lineno

    call ReadStringFromFile(string,'reduced units',SimulationInput,'SimulationInput',lineno)
    if(lineno == 0)then
      REDUCED_UNITS=.false.
    else
      REDUCED_UNITS=.true.
    end if

  end subroutine SetReducedUnits

  !=============================================================================================================================
  !**       Cut-off Distance
  !=============================================================================================================================
  subroutine SetCutOffDistance(SimulationInput,error)
    character(len=lstrlen), dimension(:), intent(in) :: SimulationInput
    logical, intent(out) :: error

    character(len=lstrlen) :: string
    integer :: lineno

    error=.false.
    call ReadStringFromFile(string,'potential cutoff distance',SimulationInput,'SimulationInput',lineno,error=error)
    if(error)return
    read(string,*)CutOffDistance
    CutOffDistanceSq=CutOffDistance**2

  end subroutine SetCutOffDistance

  !=============================================================================================================================
  !**       Simulation Variables
  !=============================================================================================================================
  subroutine ReadSimulationVariables(SimulationInput,error)
    character(len=lstrlen), dimension(:), intent(in) :: SimulationInput
    logical, intent(out) :: error

    character(len=lstrlen) :: string='',filename='', line=''
    integer :: StartLine, StartLine_2, lineno, unitno
    integer :: pos, sys, spc, spc1, spc2
    integer :: simno,SubsetNo
    character(len=strlen) :: spcname='',spc1name='',spc2name='',CellShape=''
    real(pr), dimension(3) :: boxlength, boxangle
    real(pr) :: initial_mass, initial_density
    real(PR) :: rvalue
    integer :: nmoles, natoms
    logical :: lopen

    error=.false.

    !----------------------------------------------------------------------------------------------------
    !** Real or reduced units
    call ReadStringFromFile(string,'reduced units',SimulationInput,'SimulationInput',lineno)
    if(lineno == 0)then
      REDUCED_UNITS=.false.
    else
      REDUCED_UNITS=.true.
    end if

    !----------------------------------------------------------------------------------------------------
    !** Cut-off Distance
    call ReadStringFromFile(string,'potential cutoff distance',SimulationInput,'SimulationInput',lineno,error=error)
    if(error)return
    read(string,*)CutOffDistance
    CutOffDistanceSq=CutOffDistance**2

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
    !** Number Of Equilibration Cycles
    call ReadStringFromFile(string,'number of equilibration cycles',SimulationInput,'Simulation Input',lineno,error=error)
    if(error)return
    read(string,*)NumberOfEquilibrationCycles

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
    !** Restart Frequency
    call ReadStringFromFile(string,'restart frequency',SimulationInput,'Simulation Input',lineno,error=error)
    if(error)return
    read(string,*)RestartFrequency

    !----------------------------------------------------------------------------------------------------
    !** Configs Frequency
    call ReadStringFromFile(string,'configs frequency',SimulationInput,'Simulation Input',lineno,error=error)
    if(error)return
    read(string,*)ConfigsFrequency


    !----------------------------------------------------------------------------------------------------
    !**       Simulation Cells
    !----------------------------------------------------------------------------------------------------
    StartLine=1
    scell:do
      call ReadStringFromFile(string,'simulation cell',SimulationInput(StartLine:),'Simulation Input',lineno)
      if(lineno /= 0)NumberOfSystems=NumberOfSystems+1
      if(NumberOfSystems == 0)then
        write(ErrorMessage,'(2a,i5,4x,a)')__FILE__,':',__LINE__, &
          'Tag "Simulation Cell" missing from Simulation Input'
        error=.true.
        return
      else if(NumberOfSystems > MAX_NO_OF_SYSTEMS)then
        write(ErrorMessage,'(2a,i5,4x,a,i3,a,i3)')__FILE__,':',__LINE__, &
          'Specified number Of systems ',NumberOfSystems,' exceeds MAX_NO_OF_SYSTEMS :',MAX_NO_OF_SYSTEMS
        error=.true.
        return
      else if(lineno == 0)then
        exit scell
      end if
      StartLine=StartLine+lineno

      !----------------------------------------------------------------------------------------------------
      !** Read Simulation Cell Parameters
      sys=NumberOfSystems
      call ReadStringFromFile(string,'shape',SimulationInput(StartLine:),'Simulation Input',lineno,'simulation cell',error)
      if(error)return
      !** Box Shape
      CellShape=uppercase(trim(string))
      select case (trim(CellShape))
      case ('BOX')
        !** Box Lengths
        call ReadStringFromFile(string,'box lengths',SimulationInput(StartLine:),'Simulation Input',lineno,'simulation cell',error)
        if(error)return
        read(string,*)boxlength
        !** Box Angles
        call ReadStringFromFile(string,'box angles',SimulationInput(StartLine:),'Simulation Input',lineno,'simulation cell',error)
        if(error)return
        read(string,*)boxangle
        !** Orthorhombic or Non-orthorhombic
        if(any(abs(boxangle - 90._Pr) > 1.0e-6_PR))CurrentSimulationCell(sys)%NonOrthorhombic=.true.
        !** convert angles from degrees to radians
        boxangle=boxangle*PI/180._PR
      case('CUBIC')
        CurrentSimulationCell(sys)%NonOrthorhombic=.true.
        call ReadStringFromFile(string,'density',SimulationInput(StartLine:),'Simulation Input',lineno,'simulation cell',error)
        if(error)return
        read(string,*)initial_density
        if(.not. REDUCED_UNITS)initial_density=initial_density/UNIT_DENSITY

        !** Calculate initial mass and box lengths
        StartLine_2=StartLine
        initial_mass=0._PR
        do
          call ReadStringFromFile(string,'species',SimulationInput(StartLine_2:),'Simulation Input',lineno,'simulation cell')
          if(lineno == 0)exit
          pos=ReadString(string,'name',string)
          spcname=adjustl(string)
          spc=GetSpeciesNumber(spcname,error)
          StartLine_2=StartLine_2+lineno

          call ReadStringFromFile(string,'initial number of molecules',SimulationInput(StartLine_2:),'Simulation Input',lineno,'species',error)
          if(error)return
          read(string,*)nmoles
          initial_mass=initial_mass+nmoles*Molecule(spc)%MolecularWeight
        end do
        boxlength=(initial_mass/initial_density)**(1._PR/3._PR)
        boxangle=PI/2._PR
      case default
        error=.true.
        write(ErrorMessage,'(2a,i5,4x,2a)')__FILE__,':',__LINE__, &
          'Unknown cell shape :',trim(CellShape)
        return 
      end select
      call SetSimulationCellBoxLengthsAndBoxAngles(CurrentSimulationCell(sys),boxlength,boxangle)

    end do scell

    !----------------------------------------------------------------------------------------------------
    !**       Method for initialization of coordinates
    !----------------------------------------------------------------------------------------------------
    StartLine=1
    sys=0
    scell_2:do
      call ReadStringFromFile(string,'simulation cell',SimulationInput(StartLine:),'Simulation Input',lineno)
      if(lineno /= 0)then
        sys=sys+1
      else
        exit scell_2
      end if
      StartLine=StartLine+lineno

      StartLine_2=StartLine
      do
        call ReadStringFromFile(string,'species',SimulationInput(StartLine_2:),'Simulation Input',lineno,'simulation cell')
        if(lineno == 0)exit
        pos=ReadString(string,'name',string)
        spcname=adjustl(string)
        spc=GetSpeciesNumber(spcname,error)
        natoms=Molecule(spc)%NumberOfAtoms
        StartLine_2=StartLine_2+lineno

        call ReadStringFromFile(string,'initialization method',SimulationInput(StartLine_2:),'Simulation Input',lineno,'species',error)
        if(error)return
        initmethod(spc,sys)=trim(adjustl(string))
        read(initmethod(spc,sys),*)string
        string=lowercase(string)
        string=adjustl(string)
        select case(trim(string))
        case ('random')
          call ReadStringFromFile(string,'initial number of molecules',SimulationInput(StartLine_2:),'Simulation Input',lineno,'species',error)
          if(error)return
          read(string,*)nmoles
        case ('file')
          read(initmethod(spc,sys),*)string,filename
          filename=adjustl(filename)
          unitno=GetFileUnit(trim(filename),lopen,error)
          if(error)return
          if(.not. lopen)then
            open(unit=unitno,file=trim(filename))
          end if

          call FindStringInFile(trim(Molecule(spc)%Name),unitno,lineno,line,error)
          if(error)return
          if(lineno == 0)then
            error=.true.
            write(ErrorMessage,'(2a,i4,4x,4a)')__FILE__,':',__LINE__, &
              'Failed to find coordinates of species ',trim(Molecule(spc)%Name),' in the file ',trim(filename)
            return
          end if
          read(line,*)nmoles
          close(unit=unitno)
        case default
          error=.true.
          write(ErrorMessage,'(2a,i5,4x,2a)')__FILE__,':',__LINE__, &
            'Unknown initialization method : ',trim(string)
          return
        end select
        !** Allocate memory to coordinates
        call AllocateMemoryForCoordinates(CurrentCoordinates(spc,sys),nmoles,natoms,error)
        if(error)return
        call AllocateMemoryForCoordinates(TrialCoordinates(spc,sys),nmoles,natoms,error)
        if(error)return
        call AllocateMemoryForCoordinates(OldCoordinates(spc,sys),nmoles,natoms,error)
        if(error)return
        TotalNumberOfMolecules(sys)=TotalNumberOfMolecules(sys)+nmoles
        TotalMass(sys)=TotalMass(sys)+Molecule(spc)%MolecularWeight*real(nmoles,PR)
      end do
    end do scell_2

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
      else if(sys > NumberOfSystems)then
        write(ErrorMessage,'(2a,i5,4x,a,i2,a,i2)')__FILE__,':',__LINE__, &
          'System number, ',sys,', exceeds Number of Systems :',NumberOfSystems
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
      
    !** Read Pressure
    StartLine=1
    do
      call ReadStringFromFile(string,'pressure',SimulationInput(StartLine:),'Simulation Input',lineno)
      if(lineno == 0)exit
      read(string,*)simno,sys,rvalue
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
      end if
      Pressure(simno,sys)=rvalue
      if(.not. REDUCED_UNITS)Pressure(simno,sys)=Pressure(simno,sys)/UNIT_PRESSURE
      StartLine=StartLine+lineno
    end do

    !** Read Fugacity
    StartLine=1
    do
      call ReadStringFromFile(string,'fugacity',SimulationInput(StartLine:),'Simulation Input',lineno)
      if(lineno == 0)exit
      read(string,*)spcname,simno,sys,rvalue
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
      end if
      spc=GetSpeciesNumber(adjustl(spcname),error)
      if(error)return
      Fugacity(spc,simno,sys)=rvalue
      if(.not. REDUCED_UNITS)Fugacity(spc,simno,sys)=Fugacity(spc,simno,sys)/UNIT_PRESSURE
      StartLine=StartLine+lineno
    end do

    !** Read FugacityRatio
    StartLine=1
    do
      call ReadStringFromFile(string,'fugacityratio',SimulationInput(StartLine:),'Simulation Input',lineno)
      if(lineno == 0)exit
      read(string,*)spc1name,spc2name,simno,sys,rvalue
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
      end if
      spc1=GetSpeciesNumber(adjustl(spc1name),error)
      if(error)return
      spc2=GetSpeciesNumber(adjustl(spc2name),error)
      if(error)return
      if(spc1 > spc2)call swap(spc1,spc2)
      SubsetNo=GetSubsetNumber(spc1,spc2)
      FugacityRatio(SubsetNo,simno,sys)=rvalue
      StartLine=StartLine+lineno
    end do

  end subroutine ReadSimulationVariables

  !=============================================================================================================================
  subroutine ReadAndActivateCellLists(SimulationInput,error)
    character(len=lstrlen), dimension(:), intent(in) :: SimulationInput
    logical, intent(out) :: error

    integer :: sys, lineno
    character(len=lstrlen) :: string

    error=.false.
    call ReadStringFromFile(string,'cell list',SimulationInput,'Simulation Input',lineno)
    if(lineno /= 0)then
      do sys=1,NumberOfSystems
        MainCellList(sys)%active=.true.
        read(string,*)MainCellList(sys)%MinimumCellWidth
        call InitializeCellListArrays(CurrentCoordinates(:,sys)%NumberOfMolecules,CurrentSimulationCell(sys),MainCellList(sys),error)
        if(error)return
        call CreateCellList(CurrentCoordinates(:,sys),CurrentSimulationCell(sys),MainCellList(sys))
      end do
    end if
  end subroutine ReadAndActivateCellLists
end module variables
