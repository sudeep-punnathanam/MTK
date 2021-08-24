module restart
  use consts, only: PR
  use variables, only: NumberOfSystems, CurrentCoordinates, TrialCoordinates, OldCoordinates, CurrentSimulationCell, CurrentInteractions, &
    TotalNumberOfMolecules, TotalMass
  use random, only: mtprng_state, SetRandomNumberState, GetRandomNumberState, RandomNumber
  use mcmoves, only: Translation, Rotation, Insertion, Deletion, CutAndRegrow, VolumeChange
  use averages, only: Energy_c, Pressure_c, Volume_c, Density_c, Concentration_c, Number_c, Fugacity_c, dUdKappa_c
  use atoms_and_molecules, only: Molecule, NumberOfSpecies
  use config, only: ReAllocateMemoryForCoordinates
  use utils, only: ErrorMessage, GetFileUnit
  implicit none
  public
  save
contains
  !=========================================================================================================
  subroutine WriteRestartFile(simno,cycleno,iterno,EquilibrationStage,error)
    integer, intent(in) :: simno,cycleno,iterno
    logical, intent(in) :: EquilibrationStage
    logical, intent(out) :: error

    integer :: unitno,spc,nmoles
    type(mtprng_state) :: rstate,rstate2
    logical :: lopen
    integer ::sys
    error=.false.

    unitno=GetFileUnit('RestartFile',lopen,error)
    if(error)return
    if(.not. lopen)then
      open(unit=unitno,file='RestartFile',form='unformatted')
    else
      error=.true.
      write(ErrorMessage,'(2a,i4,4x,a)')__FILE__,':',__LINE__,'Restart file already open'
      return
    end if

    !** Save Simulation Stage Information
    write(unitno)EquilibrationStage
    write(unitno)simno,cycleno,iterno

    !** Save state of the random number generator
    call GetRandomNumberState(rstate)
    write(unitno)rstate

    do sys=1,NumberOfSystems
      !** Save Simulation Cell Information
      write(unitno)CurrentSimulationCell(sys)

      !** Save Coordinates
      do spc=1,NumberOfSpecies
        nmoles=CurrentCoordinates(spc,sys)%NumberOfMolecules
        write(unitno)nmoles
        write(unitno)CurrentCoordinates(spc,sys)%Positions(:,1:nmoles,:)
        write(unitno)CurrentCoordinates(spc,sys)%CenterOfMass(:,1:nmoles)
        write(unitno)CurrentCoordinates(spc,sys)%BodyAxes(:,:,1:nmoles)
        write(unitno)CurrentCoordinates(spc,sys)%LatticePositions(:,1:nmoles)
      end do
    end do

    !** Save Move Information Statistics
    write(unitno)Translation
    write(unitno)Rotation
    write(unitno)Insertion
    write(unitno)Deletion
    write(unitno)CutAndRegrow
    write(unitno)VolumeChange

    !** Save Property Accumulators
    write(unitno)Energy_c
    write(unitno)Pressure_c
    write(unitno)Volume_c
    write(unitno)Density_c
    write(unitno)Number_c
    write(unitno)Concentration_c
    write(unitno)Fugacity_c
    write(unitno)dUdKappa_c

    close(unit=unitno)

  end subroutine WriteRestartFile
  !=========================================================================================================

  subroutine ReadRestartFile(simno,cycleno,iterno,EquilibrationStage,lexist,error)
    integer, intent(out) :: simno,cycleno,iterno
    logical, intent(out) :: EquilibrationStage
    logical, intent(out) :: lexist
    logical, intent(out) :: error

    integer :: unitno,spc,nmoles, natoms
    type(mtprng_state) :: rstate,rstate2
    logical :: lopen
    integer :: sys

    error=.false.
    inquire(file='RestartFile',exist=lexist)
    if(.not. lexist)return

    unitno=GetFileUnit('RestartFile',lopen,error)
    if(error)return
    if(.not. lopen)then
      open(unit=unitno,file='RestartFile',form='unformatted')
    else
      error=.true.
      write(ErrorMessage,'(2a,i4,4x,a)')__FILE__,':',__LINE__,'Restart file already open'
      return
    end if
   
    !** Read Simulation Stage
    read(unitno)EquilibrationStage
    read(unitno)simno,cycleno,iterno

    !** Read state of the random number generator
    read(unitno)rstate
    call SetRandomNumberState(rstate)

    do sys=1,NumberOfSystems
      !** Read Simulation Cell Information
      read(unitno)CurrentSimulationCell(sys)

      !** Read Coordinates
      TotalNumberOfMolecules(sys)=0
      TotalMass(sys)=0._PR
      do spc=1,NumberOfSpecies
        natoms=Molecule(spc)%NumberOfAtoms
        read(unitno)nmoles
        CurrentCoordinates(spc,sys)%NumberOfMolecules=nmoles
        !** Allocate memory to coordinates
        call ReAllocateMemoryForCoordinates(CurrentCoordinates(spc,sys),nmoles,natoms,error)
        if(error)return
        call ReAllocateMemoryForCoordinates(TrialCoordinates(spc,sys),nmoles,natoms,error)
        if(error)return
        call ReAllocateMemoryForCoordinates(OldCoordinates(spc,sys),nmoles,natoms,error)
        if(error)return
        TotalNumberOfMolecules(sys)=TotalNumberOfMolecules(sys)+nmoles
        TotalMass(sys)=TotalMass(sys)+Molecule(spc)%MolecularWeight*real(nmoles,PR)
        read(unitno)CurrentCoordinates(spc,sys)%Positions(:,1:nmoles,:)
        read(unitno)CurrentCoordinates(spc,sys)%CenterOfMass(:,1:nmoles)
        read(unitno)CurrentCoordinates(spc,sys)%BodyAxes(:,:,1:nmoles)
        read(unitno)CurrentCoordinates(spc,sys)%LatticePositions(:,1:nmoles)
      end do
    end do

    !** Read Move Information Statistics
    read(unitno)Translation
    read(unitno)Rotation
    read(unitno)Insertion
    read(unitno)Deletion
    read(unitno)CutAndRegrow
    read(unitno)VolumeChange

    !** Read Property Accumulators
    read(unitno)Energy_c
    read(unitno)Pressure_c
    read(unitno)Volume_c
    read(unitno)Density_c
    read(unitno)Number_c
    read(unitno)Concentration_c
    read(unitno)Fugacity_c
    read(unitno)dUdKappa_c

    close(unit=unitno)

  end subroutine ReadRestartFile
  !=========================================================================================================
end module restart
