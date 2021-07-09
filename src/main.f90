program main
  use consts, only: lstrlen, Dashed_Line
  use variables, only: SetReducedUnits, SetCutOffDistance, ReadSimulationVariables, CurrentSimulationCell, NumberOfSystems, ReadAndActivateCellLists
  use utils, only: ErrorMessage, GetFileUnit, ReadStringFromFile, uppercase
  use random, only: SetRandomNumberSeed
  use atoms_and_molecules, only: ReadAtomInfo, ReadMoleculeInfo
  use forcefield, only: ReadPairwiseInteractions, ReadIntraMolecularPotentials
  use rosenbluth, only: CalculateIdealGasRosenbluthWeight, ReadIdealGasRosenbluthWeight, ReadCBMC_Moves
  use inter, only: AllocateMemoryToSeparationVectors
  use ewaldsum, only: ReadAndConfigureEwaldSettings
  use mcmoves, only: ReadMoveInformation
  use initial, only: InitializeMoleculeCoordinates
  use averages, only: ReadPropertyCalculationInformation
  use montecarlo, only: MonteCarloSimulation

  implicit none
  
  character(len=lstrlen), dimension(1000) :: SimulationInput
  character(len=lstrlen) :: PROJECT_DIR
  integer :: NumberOfInputLines

  integer, parameter :: STANDARD_MONTE_CARLO=1, NON_INTERACTING_EINSTEIN=2, INTERACTING_EINSTEIN=3, &
    GIBBS_ENSEMBLE_MONTE_CARLO=4
  integer :: SimulationType

  integer :: LogUnitNo
  integer :: lineno, iseed, ierror
  logical :: lopen, error
  character(len=lstrlen) :: string

  error=.false.
  !------------------------------------------------------------------------------------------------------------------------
  !    Open the log file "SimulationLog" and assign a unit number
  !------------------------------------------------------------------------------------------------------------------------
  LogUnitNo=GetFileUnit('SimulationLog',lopen,error)
  if(error)then
    write(0,'(2a,i5,4x,a)')__FILE__,':',__LINE__, &
      'Failed to obtain unit number for log file "SimulationLog"'
    stop
  end if
  open(unit=LogUnitNo,file='SimulationLog',position='append',action='write',iostat=ierror)
  if(ierror /=0)then
    write(0,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
      'Failed to open log file "SimulationLog". IOSTAT =',ierror
    stop
  end if

  !------------------------------------------------------------------------------------------------------------------------
  !    Open the value of the environment variable PROJECT_DIR
  !------------------------------------------------------------------------------------------------------------------------
  call get_environment_variable(name='PROJECT_DIR',value=PROJECT_DIR,status=ierror)
  if(ierror /=0)then
    !** Environment variable PROJECT_DIR not defined. Set it to the default value, i.e. the current working directory
    PROJECT_DIR='.'
  end if
  PROJECT_DIR=adjustl(PROJECT_DIR)

  !------------------------------------------------------------------------------------------------------------------------
  !    Read Input File "SimulationInput"
  !------------------------------------------------------------------------------------------------------------------------
  call ReadInputFile(SimulationInput,NumberOfInputLines,LogUnitNo,error)
  if(error)call PrintErrorMessage(LogUnitNo)

  !------------------------------------------------------------------------------------------------------------------------
  !    Extract simulation input variables from SimulationInput
  !------------------------------------------------------------------------------------------------------------------------
  !** Simulation Type
  call ReadStringFromFile(string,'simulation type',SimulationInput,'SimulationInput',lineno,error=error)
  if(error)call PrintErrorMessage(LogUnitNo)
  string=uppercase(trim(string))
  select case(trim(string))
  case ('STANDARD MONTE CARLO')
    SimulationType=STANDARD_MONTE_CARLO
  case ('NON INTERACTING EINSTEIN')
    SimulationType=NON_INTERACTING_EINSTEIN
  case ('INTERACTING EINSTEIN')
    SimulationType=INTERACTING_EINSTEIN
  case ('GIBBS ENSEMBLE MONTE CARLO')
    SimulationType=GIBBS_ENSEMBLE_MONTE_CARLO
  case ('IDEAL GAS ROSENBLUTH FACTOR')
    call CalculateIdealGasRosenbluthWeight(PROJECT_DIR,SimulationInput,LogUnitNo,error)
    if(error) call PrintErrorMessage(LogUnitNo)
    stop
  case default
    write(ErrorMessage,'(2a,i5,4x,2a)')__FILE__,':',__LINE__, &
      'Unknown simulation type :',trim(string)
    call PrintErrorMessage(LogUnitNo)
  end select

  !------------------------------------------------------------------------------------------------------------------------
  !** Reduced Units
  call SetReducedUnits(SimulationInput(1:NumberOfInputLines))

  !------------------------------------------------------------------------------------------------------------------------
  !** Cut-off Distance
  call SetCutOffDistance(SimulationInput(1:NumberOfInputLines),error)
  if(error)call PrintErrorMessage(LogUnitNo)

  !---------------------------------------------------------------------------------------------------
  !** Set Random Number Seed
  call ReadStringFromFile(string,'random number seed',SimulationInput,'Simulation Input',lineno,error=error)
  if(error)call PrintErrorMessage(LogUnitNo)
  read(string,*)iseed
  call SetRandomNumberSeed(iseed)

  !------------------------------------------------------------------------------------------------------------------------
  !** Atoms
  call ReadAtomInfo(PROJECT_DIR,LogUnitNo,error)
  if(error)call PrintErrorMessage(LogUnitNo)

  !------------------------------------------------------------------------------------------------------------------------
  !** Molecules
  call ReadMoleculeInfo(PROJECT_DIR,SimulationInput(1:NumberOfInputLines),error)
  if(error)call PrintErrorMessage(LogUnitNo)

  !------------------------------------------------------------------------------------------------------------------------
  !** CBMC Moves
  call ReadCBMC_Moves(PROJECT_DIR,LogUnitNo,error)
  if(error)call PrintErrorMessage(LogUnitNo)

  !------------------------------------------------------------------------------------------------------------------------
  !** Forcefield
  call ReadPairwiseInteractions(PROJECT_DIR,LogUnitNo,error)
  if(error)call PrintErrorMessage(LogUnitNo)
  call ReadIntraMolecularPotentials(PROJECT_DIR,LogUnitNo,error)
  if(error)call PrintErrorMessage(LogUnitNo)
  call AllocateMemoryToSeparationVectors(error)
  if(error)call PrintErrorMessage(LogUnitNo)

  !------------------------------------------------------------------------------------------------------------------------
  !** Simulation Variables
  call ReadSimulationVariables(SimulationInput(1:NumberOfInputLines),error)
  if(error)call PrintErrorMessage(LogUnitNo)

  !------------------------------------------------------------------------------------------------------------------------
  !** Monte Carlo Moves
  call ReadMoveInformation(SimulationInput(1:NumberOfInputLines),error)
  if(error)call PrintErrorMessage(LogUnitNo)

  !------------------------------------------------------------------------------------------------------------------------
  !** Read Ideal Gas Rosenbluth Weight
  call ReadIdealGasRosenbluthWeight(SimulationInput(1:NumberOfInputLines),error)
  if(error)call PrintErrorMessage(LogUnitNo)

  !------------------------------------------------------------------------------------------------------------------------
  !** Ewald Summation
  call ReadAndConfigureEwaldSettings(SimulationInput(1:NumberOfInputLines),CurrentSimulationCell,NumberOfSystems,error)
  if(error)call PrintErrorMessage(LogUnitNo)

  !---------------------------------------------------------------------------------------------------
  !** Read Property Calculation Information
  call ReadPropertyCalculationInformation(SimulationInput(1:NumberOfInputLines),error)
  if(error)call PrintErrorMessage(LogUnitNo)

  !------------------------------------------------------------------------------------------------------------------------
  !** Initialize Coordinates
  call InitializeMoleculeCoordinates(error)
  if(error)call PrintErrorMessage(LogUnitNo)

  !------------------------------------------------------------------------------------------------------------------------
  !** Initialize Cell List
  call ReadAndActivateCellLists(SimulationInput(1:NumberOfInputLines),error)
  if(error)call PrintErrorMessage(LogUnitNo)

  !*****************************************************************************************
  !** Perform the simulation
  !*****************************************************************************************
  select case (SimulationType)
  case (STANDARD_MONTE_CARLO)
    call MonteCarloSimulation(LogUnitNo,error)
!!$  case (GIBBS_ENSEMBLE_MONTE_CARLO)
!!$    call GibbsEnsembleSimulation(LogUnitNo,RestartMode,error)
!!$  case (NON_INTERACTING_EINSTEIN)
!!$    call NonInteractingEinsteinSimulation(LogUnitNo,RestartMode,error)
!!$  case (INTERACTING_EINSTEIN)
!!$    call InteractingEinsteinSimulation(LogUnitNo,RestartMode,error)
!!$  case (IDEAL_FLEXIBLE_MOLECULE)
!!$    call idealflexmolrwfactor
  end select

  if(error)call PrintErrorMessage(LogUnitNo)

contains

  !*****************************************************************************************
  !**  Subroutine for reading the input file "SimulationInput"
  !*****************************************************************************************
  subroutine ReadInputFile(SimulationInput,NumberOfInputLines,LogUnitNo,error)
    character(len=lstrlen), dimension(:), intent(out) :: SimulationInput
    integer, intent(out)    :: NumberOfInputLines
    integer, intent(in)     :: LogUnitNo
    logical, intent(Out)    :: error

    integer                 :: lineno,ierror
    integer                 :: InputUnitNo
    character(len=lstrlen)  :: line
    logical                 :: lopen

    error=.false.
    !** Open the input file
    InputUnitNo=GetFileUnit('SimulationInput',lopen,error)
    if(error)return
    if(lopen)then
      write(ErrorMessage,'(2a,i5,4x,a)')__FILE__,':',__LINE__, &
        'Input file "SimulationInput" already open.'
      error=.true.
      return
    end if
    open(unit=InputUnitNo,file='SimulationInput',action='read',iostat=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
        'Failed to open file "SimulationInput". IOSTAT = ',ierror
      error=.true.
      return
    end if

    write(LogUnitNo,'(a)')trim(Dashed_Line)
    write(LogUnitNo,'(a)')'#                         Simulation Input'
    write(LogUnitNo,'(a)')trim(Dashed_Line)
    lineno=1
    do
      read(unit=InputUnitNo,fmt='(a)',iostat=ierror)line
      if(ierror == -1)exit
      if(ierror /= 0)then
        write(ErrorMessage,'(2a,i5,4x,a,i4,a,i3)')__FILE__,':',__LINE__, &
          'Failed to read line number ',lineno,' in file "SimulationInput". IOSTAT = ',ierror
        error=.true.
        return
      end if
      SimulationInput(lineno)=adjustl(line)
      write(LogUnitNo,'(a)')trim(adjustl(line))
      lineno=lineno+1
    end do
    write(LogUnitNo,*)
    NumberOfInputLines=lineno

    !** Close the input file
    close(unit=InputUnitNo,iostat=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to close file "SimulationInput". IOSTAT = ',ierror
      error=.true.
      return
    end if
    call flush(LogUnitNo)

  end subroutine ReadInputFile

  !*****************************************************************************************
  !** Print error messages if any to the log file "SimulationLog" and stop the program
  !*****************************************************************************************
  subroutine PrintErrorMessage(LogUnitNo)
    integer, intent(in)     :: LogUnitNo
    write(LogUnitNo,'(a)')trim(ErrorMessage)
    stop
  end subroutine PrintErrorMessage
end program main
