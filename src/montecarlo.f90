module montecarlo
  use consts, only: PR, PI, strlen, MIN_NO_OF_PARTICLE_MOVES, UNIT_DENSITY, UNIT_ENERGY, AVOGADRO_NUMBER
  use variables, only: NumberOfEquilibrationCycles, NumberOfProductionCycles, NumberOfSimulations, &
    CurrentInteractions, CurrentCoordinates, CurrentSimulationCell, UpdateFrequency, TotalNumberOfMolecules, &
    DisplayFrequency, UpdateFrequency, RestartFrequency, REDUCED_UNITS, TotalNumberOfMolecules, TotalMass, &
    CurrentKfactor, CurrentStructureFactor, MainCellList, beta
  use mcmoves, only: ParticleMove, VolumeChangeMove, UpdateMoveParameters, ResetMoveStatistics, DisplayMoveStatistics, VolumeChange
  use averages, only: SampleStatistics, ComputeStatistics, ResetStatistics, PrintStatistics
  use inter, only: TotalShortRangePairwiseInteraction, TotalLongRangePairwiseInteraction
  use utils, only: ErrorMessage, RealToString, OpenFile, CloseFile
  use random, only: RandomNumber
  use storage, only: StorageInteractions, UpdateStorageInteractions, DisplayStorage, operator(+)
  use atoms_and_molecules, only: Molecule, NumberOfSpecies
  use restart, only: ReadRestartFile, WriteRestartFile
  use ewaldsum, only: CalculateKfactors, CalculateTotalStructureFactor, EwaldFourierInteraction, TotalIntraMolecularFourierInteraction, COULOMB_INTERACTION
  implicit none
  private
  save

  public :: MonteCarloSimulation
contains
  !=================================================================================================
  subroutine MonteCarloSimulation(LogUnitNo,error)
    integer, intent(in)   :: LogUnitNo
    logical, intent(Out)  :: error

    integer :: spc
    integer :: MoveType, CycleLength, ParticleMoveWeight
    logical :: overlap, EquilibrationStage, RestartMode
    character(len=strlen) :: filename, intstr
    type(StorageInteractions) :: Interactions
    real(PR) :: energy

    integer :: simno, cycleno, iterno, sys

    error=.false.
    sys=1

    simno=1
    do while(simno <= NumberOfSimulations)
      !** Read Restart File if any
      call ReadRestartFile(simno,cycleno,iterno,EquilibrationStage,RestartMode,error)
      if(error)return
      if(RestartMode)then
        RestartMode=.false.
        if(EquilibrationStage)go to 10
        go to 20
      end if

      !** Reset Statistics
      call ResetStatistics

      !** Save initial configuration to file
      write(intstr,*)simno
      filename='InitialConfig.' // trim(adjustl(intstr))
      call WriteConfigsToFile(filename,error)
      if(error)return

      !** Calculate initial energy
      call TotalShortRangePairwiseInteraction(CurrentCoordinates(:,sys),CurrentSimulationCell(sys),MainCellList(sys),CurrentInteractions(sys),overlap,error)
      if(error)return
      if(overlap)then
        write(ErrorMessage,'(2a,i5,4x,a)')__FILE__,':',__LINE__, &
          'Failed to calculate energy of system initially. Possible overlap between atoms'
        error=.true.
        return
      end if

      call TotalLongRangePairwiseInteraction(CurrentCoordinates(:,sys)%NumberOfMolecules,CurrentSimulationCell(sys)%Volume, &
        Interactions)
      CurrentInteractions(sys)=CurrentInteractions(sys)+Interactions
      if(COULOMB_INTERACTION)then
        call CalculateKfactors(sys,CurrentSimulationCell(sys),CurrentKfactor(sys))
        call CalculateTotalStructureFactor(sys,CurrentCoordinates(:,sys),CurrentSimulationCell(sys),CurrentStructureFactor(sys),error)
        if(error)return
        call EwaldFourierInteraction(sys,CurrentSimulationCell(sys)%Volume,CurrentKfactor(sys),CurrentStructureFactor(sys),Interactions)
        CurrentInteractions(sys)=CurrentInteractions(sys)+Interactions
        call TotalIntraMolecularFourierInteraction(sys,CurrentCoordinates(:,sys),energy)
        CurrentInteractions(sys)%Energy%IntraFourier=energy
      end if
      call UpdateStorageInteractions(CurrentInteractions(sys))
      call DisplayStorage(CurrentInteractions(sys)%Energy,'INITIAL ENERGIES',LogUnitNo,REDUCED_UNITS)
      call DisplayStorage(CurrentInteractions(sys)%Virial,'INITIAL VIRIALS',LogUnitNo,REDUCED_UNITS)
      call flush(LogUnitNo)

      !** Equilibration Stage **
      EquilibrationStage=.true.
      cycleno=1
10    do while(cycleno <= NumberOfEquilibrationCycles)
        if(mod(cycleno,RestartFrequency) == 0)then
          call WriteRestartFile(simno,cycleno,iterno,EquilibrationStage,error)
          if(error)return
        end if
        ParticleMoveWeight=max(TotalNumberOfMolecules(sys),MIN_NO_OF_PARTICLE_MOVES)
        CycleLength=ParticleMoveWeight+VolumeChange(sys)%Weight
        do iterno=1,CycleLength
          MoveType=min(int(RandomNumber()*CycleLength)+1,CycleLength)
          if(MoveType <= ParticleMoveWeight)then
            call ParticleMove(simno,error)
            if(error)return
          else if(MoveType <= ParticleMoveWeight+VolumeChange(sys)%Weight)then
            call VolumeChangeMove(simno,error)
            if(error)return
          end if
        end do
        if(mod(cycleno,DisplayFrequency)==0)then
          call DisplaySimulationStatistics
          call flush(LogUnitNo)
        end if
        if(mod(cycleno,UpdateFrequency) == 0)call UpdateMoveParameters
        cycleno=cycleno+1
      end do
      call ResetMoveStatistics

      !** Production Stage **
      EquilibrationStage=.false.
      cycleno=1
20    do while(cycleno <=NumberOfProductionCycles)
        if(mod(cycleno,RestartFrequency) == 0)then
          call WriteRestartFile(simno,cycleno,iterno,EquilibrationStage,error)
          if(error)return
        end if

        ParticleMoveWeight=max(TotalNumberOfMolecules(sys),MIN_NO_OF_PARTICLE_MOVES)
        CycleLength=ParticleMoveWeight+VolumeChange(sys)%Weight
        do iterno=1,CycleLength
          MoveType=min(int(RandomNumber()*CycleLength)+1,CycleLength)
          if(MoveType <= ParticleMoveWeight)then
            call ParticleMove(simno,error)
            if(error)return
          else if(MoveType <= ParticleMoveWeight+VolumeChange(sys)%Weight)then
            call VolumeChangeMove(simno,error)
            if(error)return
          end if
        end do
        !** Sample Averages
        call SampleStatistics(cycleno,simno,error)
        if(error)return
        if(mod(cycleno,DisplayFrequency)==0)then
          call DisplaySimulationStatistics
          call flush(LogUnitNo)
        end if
        cycleno=cycleno+1
      end do

      call ComputeStatistics(simno)
      call PrintStatistics(simno,error)
      if(error)return

      !** Display Current Interactions
      call DisplayStorage(CurrentInteractions(sys)%Energy,'STORED ENERGIES',LogUnitNo,REDUCED_UNITS)
      call DisplayStorage(CurrentInteractions(sys)%Virial,'STORED VIRIALS',LogUnitNo,REDUCED_UNITS)

      !** Make a fresh estimation of Current Interactions for consistency check
      call TotalShortRangePairwiseInteraction(CurrentCoordinates(:,sys),CurrentSimulationCell(sys),MainCellList(sys),CurrentInteractions(sys),overlap,error)
      if(error)return
      if(overlap)then
        write(ErrorMessage,'(2a,i5,4x,a)')__FILE__,':',__LINE__, &
          'Failed to calculate energy of system initially. Possible overlap between atoms'
        error=.true.
        return
      end if
      call TotalLongRangePairwiseInteraction(CurrentCoordinates(:,sys)%NumberOfMolecules,CurrentSimulationCell(sys)%Volume, &
        Interactions)
      CurrentInteractions(sys)=CurrentInteractions(sys)+Interactions
      if(COULOMB_INTERACTION)then
        call CalculateKfactors(sys,CurrentSimulationCell(sys),CurrentKfactor(sys))
        call CalculateTotalStructureFactor(sys,CurrentCoordinates(:,sys),CurrentSimulationCell(sys),CurrentStructureFactor(sys),error)
        if(error)return
        call EwaldFourierInteraction(sys,CurrentSimulationCell(sys)%Volume,CurrentKfactor(sys),CurrentStructureFactor(sys),Interactions)
        CurrentInteractions(sys)=CurrentInteractions(sys)+Interactions
        call TotalIntraMolecularFourierInteraction(sys,CurrentCoordinates(:,sys),energy)
        CurrentInteractions(sys)%Energy%IntraFourier=energy
      end if
      call UpdateStorageInteractions(CurrentInteractions(sys))
      call DisplayStorage(CurrentInteractions(sys)%Energy,'CALCULATED ENERGIES',LogUnitNo,REDUCED_UNITS)
      call DisplayStorage(CurrentInteractions(sys)%Virial,'CALCULATED VIRIALS',LogUnitNo,REDUCED_UNITS)
      call flush(LogUnitNo)

      !** Save final configuration to file
      write(intstr,*)simno
      filename='FinalConfig.' // trim(adjustl(intstr))
      call WriteConfigsToFile(filename,error)
      if(error)return

      simno=simno+1
    end do

  contains
    !=============================================================================================================================
    subroutine DisplaySimulationStatistics
      character(len=strlen) :: string1,string2

      write(string1,*)simno
      string1=adjustl(string1)
      if(EquilibrationStage)then
        write(LogUnitNo,'(3a)')'SIMULATION NUMBER ',trim(string1),' :  EQUILIBRATION STAGE'
        write(string2,*)NumberOfEquilibrationCycles
        string2=adjustl(string2)
      else
        write(LogUnitNo,'(3a)')'SIMULATION NUMBER ',trim(string1),' :  PRODUCTION STAGE'
        write(string2,*)NumberOfProductionCycles
        string2=adjustl(string2)
      end if
      write(string1,*)cycleno
      string1=adjustl(string1)
      write(LogUnitNo,'(4a)')'Cycle Number ',trim(string1),' out of ',trim(string2)

      call DisplayMoveStatistics(sys,LogUnitNo)
      call DisplayProperties

      write(LogUnitNo,'(a)')'=============================================================================='
      write(LogUnitNo,*)

    end subroutine DisplaySimulationStatistics

    !=============================================================================================================================
    subroutine DisplayProperties

      integer             :: spc

      character(len=strlen), dimension(10) :: value

      do spc=1,NumberOfSpecies
        write(LogUnitNo,'(2a,t31,a,i5)')'Number of molecules of ',trim(Molecule(spc)%Name),':',CurrentCoordinates(spc,sys)%NumberOfMolecules
      end do
      if(REDUCED_UNITS)then
        write(LogUnitNo,'(a,t31,a,3f15.4)')'Box Length',':',CurrentSimulationCell(sys)%BoxLength(:)
        write(LogUnitNo,'(a,t31,a,3f15.4)')'Box Angles',':',CurrentSimulationCell(sys)%BoxAngle(:)*180._PR/PI
        write(LogUnitNo,'(a,t31,a,f15.4)')'SimulationVolume',':',CurrentSimulationCell(sys)%Volume
        write(LogUnitNo,'(a,t31,a,f15.4)')'Number Density',':',real(TotalNumberOfMolecules(sys),PR)/CurrentSimulationCell(sys)%Volume
        value(1)=RealToString(CurrentInteractions(sys)%Energy%Total)
        value(2)=RealToString(CurrentInteractions(sys)%Energy%vdW)
        value(3)=RealToString(CurrentInteractions(sys)%Energy%CoulReal)
        value(4)=RealToString(CurrentInteractions(sys)%Energy%CoulFourier)
        value(5)=RealToString(CurrentInteractions(sys)%Energy%Einstein)
        write(LogUnitNo,'(a,t31,a,a15,a,4a15,a)')'Total Energy ',':',trim(value(1)),'  ( ',trim(value(2)), trim(value(3)), trim(value(4)) , trim(value(5)) ,' )'
     else
        value(1)=RealToString(CurrentSimulationCell(sys)%BoxLength(1))
        value(2)=RealToString(CurrentSimulationCell(sys)%BoxLength(2))
        value(3)=RealToString(CurrentSimulationCell(sys)%BoxLength(3))
        write(LogUnitNo,'(a,t31,a,3a15)')'Box Length (Ang)',':',trim(value(1)),trim(value(2)),trim(value(3))
        value(1)=RealToString(CurrentSimulationCell(sys)%BoxAngle(1)*180._PR/PI)
        value(2)=RealToString(CurrentSimulationCell(sys)%BoxAngle(2)*180._PR/PI)
        value(3)=RealToString(CurrentSimulationCell(sys)%BoxAngle(3)*180._PR/PI)
        write(LogUnitNo,'(a,t31,a,3a15)')'Box Angles',':',trim(value(1)),trim(value(2)),trim(value(3))
        value(1)=RealToString(CurrentSimulationCell(sys)%Volume)
        write(LogUnitNo,'(a,t31,a,a15)')'Simulation Volume (Ang^3)',':',trim(value(1))
        value(1)=RealToString(TotalMass(sys)/CurrentSimulationCell(sys)%Volume*UNIT_DENSITY)
        write(LogUnitNo,'(a,t31,a,a15)')'Density (kg/m^3)',':',trim(value(1))
        value(1)=RealToString(CurrentInteractions(sys)%Energy%Total*UNIT_ENERGY*AVOGADRO_NUMBER)
        value(2)=RealToString(CurrentInteractions(sys)%Energy%vdW*UNIT_ENERGY*AVOGADRO_NUMBER)
        value(3)=RealToString(CurrentInteractions(sys)%Energy%CoulReal*UNIT_ENERGY*AVOGADRO_NUMBER)
        value(4)=RealToString(CurrentInteractions(sys)%Energy%CoulFourier*UNIT_ENERGY*AVOGADRO_NUMBER)
        value(5)=RealToString(CurrentInteractions(sys)%Energy%Einstein*UNIT_ENERGY*AVOGADRO_NUMBER)
        write(LogUnitNo,'(a,t31,a,a15,a,4a15,a)')'Total Energy (J/mol)',':',trim(value(1)),'  ( ',trim(value(2)), trim(value(3)), trim(value(4)) , trim(value(5)) ,' )'
     end if
    end subroutine DisplayProperties

    !=========================================================================================================
    subroutine WriteConfigsToFile(filename,error)
      character(len=*), intent(in) :: filename
      logical, intent(out) :: error

      integer :: nmoles,natoms
      integer :: spc,mol,atm
      integer :: unitno
      logical :: lopen

      error=.false.
      call OpenFile(filename,'write','append',unitno,error)
      if(error)return

      do spc=1,NumberOfSpecies
        nmoles=CurrentCoordinates(spc,sys)%NumberOfMolecules
        natoms=Molecule(spc)%NumberOfAtoms
        write(unitno,'(a,4x,i5)')trim(adjustl(Molecule(spc)%Name)),nmoles
        do mol=1,nmoles
          do atm=1,natoms
            write(unitno,100)mol,atm,CurrentCoordinates(spc,sys)%Positions(:,mol,atm)
          end do
        end do
      end do

      write(unitno,*)
      call CloseFile(unitno,error)
      if(error)return

100   format(i5,2x,i5,3(4x,f15.4))
    end subroutine WriteConfigsToFile

  end subroutine MonteCarloSimulation

  !=================================================================================================

end module montecarlo

