module averages
  use consts, only: PR, MAX_NO_OF_SYSTEMS, MAX_NO_OF_SPECIES, strlen, lstrlen, UNIT_ENERGY, UNIT_PRESSURE, UNIT_VOLUME, &
    UNIT_DENSITY, UNIT_CONCENTRATION, AVOGADRO_NUMBER
  use utils, only: ReadString, GetFileUnit, NumberFormat
  use variables, only: REDUCED_UNITS, NumberOfSystems, CurrentCoordinates, CurrentSimulationCell, CurrentInteractions, & 
    TotalNumberOfMolecules, TotalMass, beta, NumberOfProductionCycles
  use atoms_and_molecules, only: GetSpeciesNumber, NumberOfSpecies
  use storage, only: StorageInteractions, UpdateStorageInteractions
  use rosenbluth, only: CBMC_grow, CBMC_Move
  use inter, only: MoleculeSystemLongRangePairwiseInteraction
  implicit none
  private
  save
  
  !============================================================================================================
  !                           Properties
  !============================================================================================================
  !** Property Definition
  type :: Property
    character(strlen) :: Name=''
    logical          :: Calculate=.false.
    real(PR)         :: CurrentValue=0._PR
    real(PR)         :: Block=0._PR
    real(PR)         :: SquaredBlock=0._PR
    real(PR)         :: Total=0._PR
    real(PR)         :: SquaredTotal=0._PR
    real(PR)         :: BlockSquaredTotal=0._PR
    real(PR)         :: Average=0._PR
    real(PR)         :: Fluctuation=0._PR
    real(PR)         :: Error=0._PR
    integer          :: Sample=0
    integer          :: BlockLength=0
    integer          :: NumberOfBlocks=0
    integer          :: SamplingFrequency=0
  end type Property

  !============================================================================================================
  !** Histogram Definition
  type :: Histogram
    type(Property), dimension(:), allocatable :: DataPoint
    integer                                   :: NumberOfBins=0
    real(PR)                                  :: BinWidth=0._PR
  end type Histogram

  integer, parameter :: NumberOfTestInsertions=1000

  !============================================================================================================
  type(Property), dimension(MAX_NO_OF_SYSTEMS)                    :: Energy_c
  type(Property), dimension(MAX_NO_OF_SYSTEMS)                    :: Pressure_c
  type(Property), dimension(MAX_NO_OF_SYSTEMS)                    :: Volume_c
  type(Property), dimension(MAX_NO_OF_SYSTEMS)                    :: Density_c
  type(Property), dimension(MAX_NO_OF_SPECIES,MAX_NO_OF_SYSTEMS)  :: Number_c 
  type(Property), dimension(MAX_NO_OF_SPECIES,MAX_NO_OF_SYSTEMS)  :: Concentration_c 
  type(Property), dimension(MAX_NO_OF_SPECIES,MAX_NO_OF_SYSTEMS)  :: Fugacity_c
  type(Property), dimension(MAX_NO_OF_SYSTEMS)                    :: dUdKappa_c

  !============================================================================================================
  public :: Energy_c, Pressure_c, Volume_c, Density_c, Number_c, Concentration_c, Fugacity_c, dUdKappa_c
  public :: ReadPropertyCalculationInformation
  public :: UpdateStatistics, SampleStatistics, ComputeStatistics, ResetStatistics, PrintStatistics

contains
  !============================================================================================================
  subroutine ReadPropertyCalculationInformation(SimulationInput,error)
    character(len=lstrlen), dimension(:), intent(in) :: SimulationInput
    logical, intent(out) :: error

    integer :: i,spc,spc1,spc2,sys,j

    character(len=strlen)  :: spcname,spc1name,spc2name
    character(len=lstrlen) :: line, string

    error=.false.
    !** Read Property Calculation Information
    do i=1,size(SimulationInput)
      line=SimulationInput(i)
      line=adjustl(line)

      if(ReadString(line,'calculate',string) == 1)then
        line=trim(adjustl(string))
        if(ReadString(line,'pressure',string) == 1)then
          call ReadProperty(Pressure_c,'Pressure')
        else if(ReadString(line,'density',string) == 1)then
          call ReadProperty(Density_c,'Density')
        else if(ReadString(line,'volume',string) == 1)then
          call ReadProperty(Volume_c,'Volume')
        else if(ReadString(line,'energy',string) == 1)then
          call ReadProperty(Energy_c,'Energy')
        else if(ReadString(line,'number',string) == 1)then
          read(string,*)spcname
          spc=GetSpeciesNumber(spcname,error)
          if(error)return
          j=ReadString(string,trim(spcname),string)
          call ReadProperty(Number_c(spc,:),'Number_'//trim(spcname))
        else if(ReadString(line,'concentration',string) == 1)then
          read(string,*)spcname
          spc=GetSpeciesNumber(spcname,error)
          if(error)return
          j=ReadString(string,trim(spcname),string)
          call ReadProperty(Concentration_c(spc,:),'Concentration_'//trim(spcname))
        else if(ReadString(line,'fugacity',string) == 1)then
          read(string,*)spcname
          spc=GetSpeciesNumber(spcname,error)
          j=ReadString(string,trim(spcname),string)
          if(error)return
          call ReadProperty(Fugacity_c(spc,:),'Fugacity_'//trim(spcname))
        else if(ReadString(line,'dudkappa',string) == 1)then
          call ReadProperty(dUdKappa_c,'dUdKappa')
         end if
      end if
    end do

    contains
      subroutine ReadProperty(obj,objname)
        type(Property), dimension(:), intent(inout) :: obj
        character(len=*), intent(in)                :: objname
        obj(1)%Calculate=.true.
        obj(1)%Name=trim(adjustl(objname))
        read(string,*)obj(1)%SamplingFrequency
        obj(1)%BlockLength=NumberOfProductionCycles/obj(1)%SamplingFrequency/5
        obj(2:)=obj(1)
      end subroutine ReadProperty
  end subroutine ReadPropertyCalculationInformation

  subroutine UpdateStatistics(obj)
    type(Property), intent(INOUT) :: obj

    real(PR) :: BlockAverage
    
    obj%Sample=obj%Sample+1
    obj%Block=obj%Block+obj%CurrentValue
    obj%SquaredBlock=obj%SquaredBlock+(obj%CurrentValue)**2

    !** Block Statistics
    if(obj%Sample == obj%BlockLength)then
      obj%NumberOfBlocks=obj%NumberOfBlocks+1
      BlockAverage=obj%Block/real(obj%BlockLength,PR)

      obj%Total=obj%Total+obj%Block
      obj%SquaredTotal=obj%SquaredTotal+obj%SquaredBlock
      obj%BlockSquaredTotal=obj%BlockSquaredTotal+(BlockAverage)**2

      obj%Sample=0
      obj%Block=0._PR
      obj%SquaredBlock=0._PR
    end if
  end subroutine UpdateStatistics

 subroutine SampleStatistics(cycleno,simno,error)
   integer, intent(in) :: cycleno,simno
   logical, intent(out) :: error

   integer :: i,j

   error=.false.
   call CalculateEnergy(cycleno)
   call CalculatePressure(cycleno,simno)
   call CalculateVolume(cycleno)
   call CalculateDensity(cycleno)
   call CalculatedUdKappa(cycleno)
   if(error)return
   do i=1,NumberOfSpecies
     call CalculateNumber(i,cycleno)
     call CalculateConcentration(i,cycleno)
     call CalculateFugacity(i,cycleno,simno,error)
     if(error)return
   end do

 end subroutine SampleStatistics

  subroutine ComputeStatistics(simno)
    integer, intent(in) :: simno

    integer :: i,j, sys
    real(PR) :: kT
        
    do sys=1,NumberOfSystems
      kT=1._PR/beta(simno,sys)
      call CalculateAverages(Energy_c(sys))
      call CalculateAverages(Pressure_c(sys))
      call CalculateAverages(Volume_c(sys))
      call CalculateAverages(Density_c(sys))
      call CalculateAverages(dUdKappa_c(sys))

      do i=1,NumberOfSpecies
        call CalculateAverages(Number_c(i,sys))
        call CalculateAverages(Concentration_c(i,sys))
        call CalculateAverages(Fugacity_c(i,sys))
        if(Fugacity_c(i,sys)%Calculate)then
          Fugacity_c(i,sys)%Error=abs(Fugacity_c(i,sys)%Error/Fugacity_c(i,sys)%Average**2)
          Fugacity_c(i,sys)%Fluctuation=abs(Fugacity_c(i,sys)%Fluctuation/Fugacity_c(i,sys)%Average**2)
          Fugacity_c(i,sys)%Average=kT/Fugacity_c(i,sys)%Average
        end if
      end do
    end do

  contains
    subroutine CalculateAverages(obj)

      type(Property), intent(INOUT) :: obj
      if(.not. obj%Calculate)return
      obj%Average=obj%Total/real(obj%NumberOfBlocks*obj%BlockLength,PR)

      obj%Fluctuation=obj%SquaredTotal/real(obj%NumberOfBlocks*obj%BlockLength,PR)-(obj%Average)**2
      obj%Fluctuation=sqrt(obj%Fluctuation)

      obj%Error=(obj%BlockSquaredTotal-real(obj%NumberOfBlocks,PR)*(obj%Average)**2)/real(obj%NumberOfBlocks-1,PR)
      obj%Error=obj%Error/real(obj%NumberOfBlocks,PR)
      !** 90% Confidence interval
      obj%Error=2.13184679_PR*sqrt(obj%Error)  

    end subroutine CalculateAverages
  end subroutine ComputeStatistics

  subroutine ResetStatistics
    integer :: i,j, sys
   
    do sys=1,NumberOfSystems
      call ResetAverages(Energy_c(sys))
      call ResetAverages(Pressure_c(sys))
      call ResetAverages(Volume_c(sys))
      call ResetAverages(Density_c(sys))
      call ResetAverages(dUdKappa_c(sys))
      do i=1,NumberOfSpecies
        call ResetAverages(Number_c(i,sys))
        call ResetAverages(Concentration_c(i,sys))
        call ResetAverages(Fugacity_c(i,sys))
      end do
    end do

  contains
    subroutine ResetAverages(obj)

      type(Property), intent(INOUT) :: obj

      if(.not. obj%Calculate)return
      obj%Block=0._PR
      obj%SquaredBlock=0._PR
      obj%Total=0._PR
      obj%SquaredTotal=0._PR
      obj%BlockSquaredTotal=0._PR
      obj%Sample=0
      obj%NumberOfBlocks=0
    end subroutine ResetAverages
  end subroutine ResetStatistics

  subroutine PrintStatistics(simno,error)
    integer, intent(in) :: simno
    logical, intent(out) :: error
    
    integer :: i,j, sys
    error=.false.

    do sys=1,NumberOfSystems
      call PrintAverages(Energy_c(sys),UNIT_ENERGY*AVOGADRO_NUMBER,'J/mol',error)
      if(error)return
      call PrintAverages(Pressure_c(sys),UNIT_PRESSURE/1000._PR,'kPa',error)
      if(error)return
      call PrintAverages(Volume_c(sys),UNIT_VOLUME,'m3',error)
      if(error)return
      call PrintAverages(Density_c(sys),UNIT_DENSITY,'kg per m3',error)
      if(error)return
      call PrintAverages(dUdKappa_c(sys),UNIT_ENERGY*AVOGADRO_NUMBER,'J/mol',error)
      if(error)return
      do i=1,NumberOfSpecies
        call PrintAverages(Number_c(i,sys),1._PR,' ',error)
        if(error)return
        call PrintAverages(Concentration_c(i,sys),UNIT_CONCENTRATION,'per m3',error)
        if(error)return
        call PrintAverages(Fugacity_c(i,sys),UNIT_PRESSURE/1000._PR,'kPa',error)
        if(error)return
      end do
    end do
  contains

    subroutine PrintAverages(obj,unitvalue,unitstr,error)
      type(Property), intent(INOUT) :: obj
      real(PR), intent(IN)          :: unitvalue
      character(len=*), intent(IN)  :: unitstr
      logical, intent(out) :: error

      integer               :: unitno
      character(len=strlen) :: fmtstr
      character(len=strlen) :: filename
      character(len=strlen) :: avg,err,fluc
      logical               :: lopen,lexist

      if(.not. obj%Calculate)return

      error=.false.
      filename=trim(adjustl(obj%Name)) // '.out'
      unitno=GetFileUnit(trim(adjustl(filename)),lopen,error)
      if(error)return
      inquire(file=trim(adjustl(filename)),exist=lexist)
      if(.not. lexist)then
        open(unit=unitno,file=trim(adjustl(filename)))
        if(REDUCED_UNITS)then
          write(unitno,'(a,t20,a,t40,a,t60,a)')'Number','Value','Error','Fluctuation'
        else
          write(unitno,'(a,t20,3a,t40,3a,t60,3a)')'Number','Value (',trim(adjustl(unitstr)),')','Error (',trim(adjustl(unitstr)),&
            ')','Fluctuation (',trim(adjustl(unitstr)),')'
        end if
        write(unitno,'(a)')'================================================================================='
      else
        open(unit=unitno,file=trim(adjustl(filename)),position='append')
      end if

      if(.not. REDUCED_UNITS)then
        obj%Average=obj%Average*unitvalue
        obj%Error=obj%Error*unitvalue
        obj%Fluctuation=obj%Fluctuation*unitvalue
      end if
      fmtstr=NumberFormat(obj%Average)
      write(avg,fmtstr)obj%Average
      write(err,fmtstr)obj%Error
      write(fluc,fmtstr)obj%Fluctuation
      write(unitno,'(i3,t20,a,t40,a,t60,a)')simno,trim(adjustl(avg)),trim(adjustl(err)),trim(adjustl(fluc))

      close(unit=unitno)

    end subroutine PrintAverages
  end subroutine PrintStatistics

  subroutine CalculateEnergy(cycleno)
    integer, intent(in) :: cycleno
    integer                    :: nmoles,sys

    do sys=1,NumberOfSystems
      if(.not. Energy_c(sys)%Calculate)cycle
      if(mod(cycleno,Energy_c(sys)%SamplingFrequency) /= 0)cycle
      Energy_c(sys)%CurrentValue=CurrentInteractions(sys)%Energy%Total/real(TotalNumberOfMolecules(sys),PR)
      call UpdateStatistics(Energy_c(sys))
    end do

  end subroutine CalculateEnergy

  subroutine CalculatePressure(cycleno,simno)
    integer, intent(in) :: cycleno
    integer, intent(in) :: simno

    real(PR)            :: kT
    integer             :: nmoles
    integer             ::sys

    do sys=1,NumberOfSystems
      if(.not. Pressure_c(sys)%Calculate)cycle
      if(mod(cycleno,Pressure_c(sys)%SamplingFrequency) /= 0)cycle
   
      nmoles=TotalNumberOfMolecules(sys)
      kT=1._PR/beta(simno,sys)
      Pressure_c(sys)%CurrentValue=(kT*real(nmoles,PR) + CurrentInteractions(sys)%Virial%Total)/CurrentSimulationCell(sys)%Volume

      call UpdateStatistics(Pressure_c(sys))
    end do

  end subroutine CalculatePressure

  subroutine CalculateVolume(cycleno)
    integer, intent(in) :: cycleno
    integer             ::sys

    do sys=1,NumberOfSystems
      if(.not. Volume_c(sys)%Calculate)cycle
      if(mod(cycleno,Volume_c(sys)%SamplingFrequency) /= 0)cycle

      Volume_c(sys)%CurrentValue=CurrentSimulationCell(sys)%Volume

      call UpdateStatistics(Volume_c(sys))
    end do

  end subroutine CalculateVolume

  subroutine CalculateDensity(cycleno)
    integer, intent(in) :: cycleno
    integer             ::sys

    do sys=1,NumberOfSystems
      if(.not. Density_c(sys)%Calculate)cycle
      if(mod(cycleno,Density_c(sys)%SamplingFrequency) /= 0)cycle

      if(REDUCED_UNITS)then
          Density_c(sys)%CurrentValue=real(TotalNumberOfMolecules(sys),PR)/CurrentSimulationCell(sys)%Volume
      else
        Density_c(sys)%CurrentValue=TotalMass(sys)/CurrentSimulationCell(sys)%Volume
      end if

      call UpdateStatistics(Density_c(sys))
    end do

  end subroutine CalculateDensity

  subroutine CalculateNumber(spc,cycleno)
    integer, intent(IN) :: spc
    integer, intent(in) :: cycleno
    integer::sys

    do sys=1,NumberOfSystems
      if(.not. Number_c(spc,sys)%Calculate)cycle
      if(mod(cycleno,Number_c(spc,sys)%SamplingFrequency) /= 0)cycle

      Number_c(spc,sys)%CurrentValue=real(CurrentCoordinates(spc,sys)%NumberOfMolecules,PR)
      
      call UpdateStatistics(Number_c(spc,sys))
    end do

  end subroutine CalculateNumber
 

  subroutine CalculateConcentration(spc,cycleno)
    integer, intent(IN) :: spc
    integer, intent(in) :: cycleno
    integer             :: nmoles
    integer             ::sys

    do sys=1,NumberOfSystems
      if(.not. Concentration_c(spc,sys)%Calculate)cycle
      if(mod(cycleno,Concentration_c(spc,sys)%SamplingFrequency) /= 0)cycle

      nmoles=CurrentCoordinates(spc,sys)%NumberOfMolecules
      Concentration_c(spc,sys)%CurrentValue=real(nmoles,PR)/CurrentSimulationCell(sys)%Volume

      call UpdateStatistics(Concentration_c(spc,sys))
    end do

  end subroutine CalculateConcentration

  subroutine CalculateFugacity(spc,cycleno,simno,error)
    integer, intent(in)  :: spc,cycleno,simno
    logical, intent(out) :: error

    integer :: sys,i,nmoles,SequenceNumber
    real(PR) :: AverageRosenbluthWeight, RosenbluthWeight, IdealRwFactor
    logical :: overlap
    type(StorageInteractions) :: Interactions, LongRangeInteractions

    error=.false.
    do sys=1,NumberOfSystems
      if(.not. Fugacity_c(spc,sys)%Calculate)cycle
      if(mod(cycleno,Fugacity_c(spc,sys)%SamplingFrequency) /= 0)cycle

      nmoles=CurrentCoordinates(spc,sys)%NumberOfMolecules+1

      AverageRosenbluthWeight=0._pr
      do i=1,NumberOfTestInsertions
        SequenceNumber=1
        RosenbluthWeight=CBMC_grow(SequenceNumber,sys,spc,nmoles,beta(simno,sys),.false.,overlap)
        if(overlap)cycle
        AverageRosenbluthWeight=AverageRosenbluthWeight+RosenbluthWeight
      end do
      AverageRosenbluthWeight=AverageRosenbluthWeight/real(NumberOfTestInsertions,PR)

      !** Calculate Long Range van der Waals and Ewald Fourier Interactions
      call MoleculeSystemLongRangePairwiseInteraction(spc,nmoles+1,CurrentCoordinates(:,sys)%NumberOfMolecules,CurrentSimulationCell(sys)%Volume,LongRangeInteractions)
      call UpdateStorageInteractions(LongRangeInteractions)
      AverageRosenbluthWeight=AverageRosenbluthWeight*exp(-beta(simno,sys)*LongRangeInteractions%Energy%Total)

      IdealRwFactor=CBMC_Move(spc)%GrowthSequence(SequenceNumber)%IdealRosenbluthWeight
      Fugacity_c(spc,sys)%CurrentValue=(AverageRosenbluthWeight/IdealRwFactor)*CurrentSimulationCell(sys)%Volume/real(nmoles+1,PR)
      call UpdateStatistics(Fugacity_c(spc,sys))
    end do

  end subroutine CalculateFugacity

  subroutine CalculatedUdKappa(cycleno)
    integer, intent(in) :: cycleno
    integer             :: nmoles,sys

    do sys=1,NumberOfSystems
      if(.not. dUdKappa_c(sys)%Calculate)cycle
      if(mod(cycleno,dUdKappa_c(sys)%SamplingFrequency) /= 0)cycle
      dUdKappa_c(sys)%CurrentValue=CurrentInteractions(sys)%dUdKappa%Total/real(TotalNumberOfMolecules(sys),PR)
      call UpdateStatistics(dUdKappa_c(sys))
    end do

  end subroutine CalculatedUdKappa

end module averages
