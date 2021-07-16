module initial
  use consts, only: PR, xlstrlen, strlen, INCR_SIZE
  use utils, only: ErrorMessage, FindStringInFile, lowercase, RotateVectors, OpenFile, CloseFile
  use atoms_and_molecules, only: NumberOfSpecies, Molecule, GetCenterOfMass
  use variables, only: NumberOfSystems, CurrentCoordinates, TrialCoordinates, initmethod, &
    CurrentSimulationCell, MainCellList, beta
  use config, only: SetCenterOfMass
  use storage, only: StorageInteractions
  use inter, only: MoleculeSystemShortRangePairwiseInteraction
  use simcell, only: RandomPointInCell
  use random, only: UniformRandomRotationMatrix
  use rosenbluth, only: CBMC_grow
  implicit none
  private
  save

  public :: InitializeMoleculeCoordinates

contains

  !============================================================================================================
  subroutine InitializeMoleculeCoordinates(error)
    logical, intent(Out)  :: error

    character(len=strlen) :: string, method, filename
    integer :: sys,spc, mol, nmoles
    logical :: overlap

    error=.false.
    do sys=1,NumberOfSystems
      do spc=1,NumberOfSpecies
        !** Initialize coordinates
        read(initmethod(spc,sys),*)string
        method=lowercase(string)
        method=adjustl(method)
        select case(trim(method))
        case ('random')
          !** do nothing at present
        case ('file')
          read(initmethod(spc,sys),*)string,filename
          filename=adjustl(filename)
          call ReadConfigsFromFile(sys,spc,filename,error)
          if(error)return
        case default
          error=.true.
          write(ErrorMessage,'(2a,i5,4x,2a)')__FILE__,':',__LINE__, &
            'Unknown initialization method : ',trim(method)
          return
        end select
      end do

      !** Random Initialization
      do spc=1,NumberOfSpecies
        read(initmethod(spc,sys),*)string
        method=lowercase(string)
        method=adjustl(method)
        if(trim(method) == 'random')then
          CurrentCoordinates(spc,sys)%NumberOfMolecules=0
          mol=1
          nmoles=size(CurrentCoordinates(spc,sys)%Positions,2)
          do while(mol <= nmoles)
            call InsertMolecule(sys,spc,overlap,error)
            if(error)return
            if(overlap)cycle
            CurrentCoordinates(spc,sys)%NumberOfMolecules=mol
            mol=mol+1
          end do
        end if
      end do
    end do
    
  end subroutine InitializeMoleculeCoordinates

  !=========================================================================================================
  subroutine ReadConfigsFromFile(sys,spc,filename,error)
    integer, intent(in) :: sys,spc
    character(len=*), intent(in) :: filename
    logical, intent(out) :: error

    integer :: nmoles,natoms
    integer :: mol,atm
    integer :: unitno,lineno
    real(PR) :: vec(3)
    logical :: lopen
    character(len=xlstrlen) :: line
    integer :: im,ia
    
    error=.false.
    call OpenFile(filename,'read','rewind',unitno,error)
    if(error)return

    call FindStringInFile(trim(Molecule(spc)%Name),unitno,lineno,line,error)
    if(error)return
    if(lineno == 0)then
      error=.true.
      write(ErrorMessage,'(2a,i4,4x,4a)')__FILE__,':',__LINE__, &
        'Failed to find coordinates of species ',trim(Molecule(spc)%Name),' in the file ',trim(filename)
      return
    end if
    read(line,*)nmoles
    natoms=Molecule(spc)%NumberOfAtoms
    do mol=1,nmoles
      do atm=1,natoms
        read(unitno,100)im,ia,CurrentCoordinates(spc,sys)%Positions(:,mol,atm)
      end do
      call SetCenterOfMass(CurrentCoordinates(spc,sys),spc,mol)
    end do
    call CloseFile(unitno,error)
    if(error)return

!!$    !** Einstein Crystal
!!$    if(SimulationType == NON_INTERACTING_EINSTEIN .or. SimulationType == INTERACTING_EINSTEIN)then
!!$      do mol=1,nmoles
!!$        CurrentCoordinates(spc,sys)%LatticePosition(:,mol)=CurrentCoordinates(spc,sys)%CenterOfMass(:,mol)
!!$        select case(trim(Molecule(spc)%Shape))
!!$        case ('spherical')
!!$          !** do nothing
!!$        case ('linear')
!!$          vec=CurrentCoordinates(spc,sys)%Positions(:,mol,1)-CurrentCoordinates(spc,sys)%Positions(:,mol,2)
!!$          vec=vec/sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
!!$          InitialCoordinates(spc,sys)%BodyAxes(:,1,mol)=vec
!!$          CurrentCoordinates(spc,sys)%BodyAxes(:,1,mol)=vec
!!$        case ('non-linear')
!!$          vec=2._PR*CurrentCoordinates(spc,sys)%Positions(:,mol,1)- &
!!$            (CurrentCoordinates(spc,sys)%Positions(:,mol,2)+CurrentCoordinates(spc,sys)%Positions(:,mol,3))
!!$          vec=vec/sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
!!$          InitialCoordinates(spc,sys)%BodyAxes(:,2,mol)=vec
!!$          CurrentCoordinates(spc,sys)%BodyAxes(:,2,mol)=vec
!!$          InitialCoordinates(spc,sys)%BodyAxes(1,1,mol)=vec(2)/sqrt(vec(1)**2+vec(2)**2)
!!$          InitialCoordinates(spc,sys)%BodyAxes(2,1,mol)=-vec(1)/sqrt(vec(1)**2+vec(2)**2)
!!$          InitialCoordinates(spc,sys)%BodyAxes(3,1,mol)=0._PR
!!$          CurrentCoordinates(spc,sys)%BodyAxes(:,1,mol)=InitialCoordinates(spc,sys)%BodyAxes(:,1,mol)
!!$        case default
!!$          error=.true.
!!$          write(ErrorMessage,'(2a,i4,4x,2a)')__FILE__,':',__LINE__, &
!!$            'Unknown Molecule shape: ',trim(Molecule(spc)%Shape)
!!$          return
!!$        end select
!!$      end do
!!$    end if

100 format(i5,2x,i5,3(4x,f15.4))
  end subroutine ReadConfigsFromFile

  !========================================================================================================================================================================
  subroutine InsertMolecule(sys,spc,overlap,error)
    integer, intent(in)      :: sys,spc
    logical, intent(out)     :: overlap
    logical, intent(out)     :: error

    type(StorageInteractions) :: Interactions

    real(PR), dimension(3) :: vec
    real(PR), dimension(3,3) :: RotationMatrix
    real(PR) :: RosenbluthWeight
    integer :: mol,natoms, SequenceNumber
    logical :: existing_molecule

    error=.false.
    overlap=.false.

    mol=CurrentCoordinates(spc,sys)%NumberOfMolecules+1
    natoms=Molecule(spc)%NumberOfAtoms

    !--------------------------------------------------------------------------------------------------------------------------
    !** Create Trial Configuration
    !--------------------------------------------------------------------------------------------------------------------------
    !** Calculate Change in Energy and Rosenbluth Weight
    existing_molecule=.false.
    SequenceNumber=1
    RosenbluthWeight=CBMC_grow(SequenceNumber,sys,spc,mol,beta(1,sys),existing_molecule,overlap)
    if(overlap)return
    CurrentCoordinates(spc,sys)%Positions(:,mol,:)=TrialCoordinates(spc,sys)%Positions(:,mol,:)
    CurrentCoordinates(spc,sys)%CenterOfMass(:,mol)=GetCenterOfMass(spc,TrialCoordinates(spc,sys)%Positions(:,mol,:))

  end subroutine InsertMolecule

  !============================================================================================================
end module initial
