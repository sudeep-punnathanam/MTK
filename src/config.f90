module config
  use consts, only: PR, INCR_SIZE
  use utils, only: ErrorMessage
  use atoms_and_molecules, only: Molecule, Atom
  implicit none
  public
  save

  !============================================================================================================
  !** Atomic Coordinates, Center of Mass, etc.
  type :: SpeciesCoordinates
    integer                                         :: NumberOfMolecules=0
    real(PR), dimension(:,:,:), allocatable         :: Positions
    real(PR), dimension(:,:), allocatable           :: CenterOfMass
    real(PR), dimension(:,:), allocatable           :: LatticePositions
    real(PR), dimension(:,:,:), allocatable         :: BodyAxes
  end type SpeciesCoordinates
  !============================================================================================================

contains
  !=============================================================================
  subroutine AllocateMemoryForCoordinates(coords,nmoles,natoms,error)
    type(SpeciesCoordinates), intent(inout) :: coords
    integer, intent(in) :: nmoles, natoms
    logical, intent(out) :: error

    integer :: ierror

    error=.false.
    coords%NumberOfMolecules=nmoles
    allocate( &
      coords%Positions(3,nmoles,natoms), &
      coords%CenterOfMass(3,nmoles), &
      coords%BodyAxes(3,2,nmoles), &
      coords%LatticePositions(3,nmoles), &
      STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to allocate memory to Arrays. STAT = ',ierror
      error=.true.
      return
    end if

  end subroutine AllocateMemoryForCoordinates

  !=============================================================================
  subroutine ReAllocateMemoryForCoordinates(coords,ArraySize,natoms,error)
    type(SpeciesCoordinates), intent(inout) :: coords
    integer, intent(IN)  :: ArraySize,natoms
    logical, intent(OUT) :: error

    real(PR), dimension(:,:,:), allocatable :: temp
    integer :: nmoles
    integer :: NewArraySize
    integer :: ierror
    integer :: sys

    error=.false.
    if(ArraySize <= size(coords%Positions,2))return

    nmoles=coords%NumberOfMolecules
    NewArraySize=((ArraySize-1)/INCR_SIZE+1)*INCR_SIZE

    !** Create temp
    allocate(temp(3,nmoles,natoms),STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to allocate memory to variable temp. STAT = ',ierror
      error=.true.
      return
    end if

    !** Positions
    temp=coords%Positions(:,1:nmoles,1:natoms)
    deallocate(coords%Positions,STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to deallocate memory from variable Positions. STAT = ',ierror
      error=.true.
      return
    end if
    allocate(coords%Positions(3,NewArraySize,natoms),STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to allocate memory to variable Positions. STAT = ',ierror
      error=.true.
      return
    end if
    coords%Positions(:,1:nmoles,1:natoms)=temp

    !** Center of Mass
    temp(:,1:nmoles,1)=coords%CenterOfMass(:,1:nmoles)
    deallocate(coords%CenterOfMass,STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to deallocate memory from variable CenterOfMass. STAT = ',ierror
      error=.true.
      return
    end if
    allocate(coords%CenterOfMass(3,NewArraySize),STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to allocate memory to variable CenterOfMass. STAT = ',ierror
      error=.true.
      return
    end if
    coords%CenterOfMass(:,1:nmoles)=temp(:,:,1)

    !** Lattice Positions
    temp(:,1:nmoles,1)=coords%LatticePositions(:,1:nmoles)
    deallocate(coords%LatticePositions,STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to deallocate memory from variable LatticePositions. STAT = ',ierror
      error=.true.
      return
    end if
    allocate(coords%LatticePositions(3,NewArraySize),STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to allocate memory to variable LatticePositions. STAT = ',ierror
      error=.true.
      return
    end if
    coords%LatticePositions(:,1:nmoles)=temp(:,:,1)

    !** Resize temp
    deallocate(temp,STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to deallocate memory from variable temp. STAT = ',ierror
      error=.true.
      return
    end if
    allocate(temp(3,2,nmoles),STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to deallocate memory from variable temp. STAT = ',ierror
      error=.true.
      return
    end if

    !** Body Axes
    temp=coords%BodyAxes(:,:,1:nmoles)
    deallocate(coords%BodyAxes,STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to deallocate memory from variable BodyAxes. STAT = ',ierror
      error=.true.
      return
    end if
    allocate(coords%BodyAxes(3,2,NewArraySize),STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to allocate memory to variable BodyAxes. STAT = ',ierror
      error=.true.
      return
    end if
    coords%BodyAxes(:,:,1:nmoles)=temp

    !** deallocate temp
    deallocate(temp,STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to deallocate memory from variable temp. STAT = ',ierror
      error=.true.
      return
    end if
  end subroutine ReAllocateMemoryForCoordinates

  !=============================================================================
  subroutine SetCenterOfMass(coords,spc,mol)
    type(SpeciesCoordinates), intent(inout) :: coords
    integer, intent(IN)                     :: spc,mol

    integer :: atm,atype
    real(PR) :: mass

    coords%CenterOfMass(:,mol)=0._PR
    mass=0._PR
    do atm=1,Molecule(spc)%NumberOfAtoms
      atype=Molecule(spc)%AtomType(atm)
      coords%CenterOfMass(:,mol)=coords%CenterOfMass(:,mol)+coords%Positions(:,mol,atm)*Atom(atype)%Mass
      mass=mass+Atom(atype)%Mass
    end do
    coords%CenterOfMass(:,mol)=coords%CenterOfMass(:,mol)/mass

  end subroutine SetCenterOfMass

end module config
