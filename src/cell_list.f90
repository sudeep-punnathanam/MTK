module cell_list
  use consts, only: PR, INCR_SIZE, MAX_NO_OF_SPECIES
  use atoms_and_molecules, only: Molecule, NumberOfSpecies
  use simcell, only: SimcellInfo, ApplyBoundaryCondition
  use config, only: SpeciesCoordinates
  use utils, only: ErrorMessage, Matrix_x_Vector

  implicit none
  public
  save

  type :: SpeciesCellListInfo
    integer, dimension(:,:), allocatable :: HeadOfCell, LinkedList
  end type SpeciesCellListInfo

  type :: CellListInfo
    logical :: active=.false.
    integer, dimension(3) :: NumberOfCells=0
    integer :: TotalNumberOfCells=0
    real(PR), dimension(3) :: CellWidth, CellLength
    real(PR) :: MinimumCellWidth
    integer, dimension(:,:), allocatable :: CellNeighbors
    type(SpeciesCellListInfo), dimension(MAX_NO_OF_SPECIES) :: SpeciesCellList
  end type CellListInfo
  
contains
  subroutine InitializeCellListArrays(NumberOfMolecules,SimulationCell,CellList,error)
    integer, dimension(:), intent(in) :: NumberOfMolecules
    type(SimcellInfo), intent(in)     :: SimulationCell
    type(CellListInfo), intent(inout) :: CellList
    logical, intent(Out)              :: error

    
    integer :: ic,ix,iy,iz
    integer :: nc,nx,ny,nz
    integer :: inx,iny,inz
    integer :: spc, nmoles, natoms, ncells
    integer :: ierror

    error=.false.

    CellList%NumberOfCells=int(SimulationCell%BoxWidth/CellList%MinimumCellWidth)
    if(any(CellList%NumberOfCells <= 3))then
      CellList%active=.false.
      return
    end if
    CellList%TotalNumberOfCells=product(CellList%NumberOfCells)

    ncells=CellList%TotalNumberOfCells
    do spc=1,NumberOfSpecies
      nmoles=NumberOfMolecules(spc)
      natoms=Molecule(spc)%NumberOfAtoms

      allocate(CellList%SpeciesCellList(spc)%HeadOfCell(ncells,natoms),STAT=ierror)
      if(ierror /= 0)then
        write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
          'Failed to allocate memory to variable HeadOfCell. STAT = ',ierror
        error=.true.
        return
      end if
      allocate(CellList%SpeciesCellList(spc)%LinkedList(nmoles,natoms),STAT=ierror)
      if(ierror /= 0)then
        write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
          'Failed to allocate memory to variable LinkedList. STAT = ',ierror
        error=.true.
        return
      end if
    end do

    allocate(CellList%CellNeighbors(27,ncells), STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to allocate memory to variable CellNeighbors. STAT = ',ierror
      error=.true.
      return
    end if

    do ic=1,CellList%TotalNumberOfCells
      iz=(ic-1)/CellList%NumberOfCells(1)/CellList%NumberOfCells(2)
      iy=((ic-1)-iz*CellList%NumberOfCells(1)*CellList%NumberOfCells(2))/CellList%NumberOfCells(1)
      ix=(ic-1)-iz*CellList%NumberOfCells(1)*CellList%NumberOfCells(2)-iy*CellList%NumberOfCells(1)
      nc=1
      do nz=-1,1
        inz=iz+nz
        if(inz == -1)inz=CellList%NumberOfCells(3)-1
        if(inz == CellList%NumberOfCells(3))inz=0
        do ny=-1,1
          iny=iy+ny
          if(iny == -1)iny=CellList%NumberOfCells(2)-1
          if(iny == CellList%NumberOfCells(2))iny=0
          do nx=-1,1
            inx=ix+nx
            if(inx == -1)inx=CellList%NumberOfCells(1)-1
            if(inx == CellList%NumberOfCells(1))inx=0
            CellList%CellNeighbors(nc,ic)=1+inx+iny*CellList%NumberOfCells(1)+inz*CellList%NumberOfCells(1)*CellList%NumberOfCells(2)
            nc=nc+1
          end do
        end do
      end do
    end do

  end subroutine InitializeCellListArrays

  subroutine DeallocateCellListArrays(CellList,error)
    type(CellListInfo), intent(inout) :: CellList
    logical, intent(Out)              :: error

    integer :: spc, ierror

    error=.false.

    do spc=1,NumberOfSpecies
      deallocate(CellList%SpeciesCellList(spc)%HeadOfCell,STAT=ierror)
      if(ierror /= 0)then
        write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
          'Failed to deallocate memory from variable HeadOfCell. STAT = ',ierror
        error=.true.
        return
      end if
      deallocate(CellList%SpeciesCellList(spc)%LinkedList,STAT=ierror)
      if(ierror /= 0)then
        write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
          'Failed to deallocate memory from variable LinkedList. STAT = ',ierror
        error=.true.
        return
      end if
    end do

    deallocate(CellList%CellNeighbors, STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to deallocate memory from variable CellNeighbors. STAT = ',ierror
      error=.true.
      return
    end if

  end subroutine DeallocateCellListArrays

  subroutine ResizeCellListArrays(NumberOfMolecules,SimulationCell,CellList,error)
    integer, dimension(:), intent(in) :: NumberOfMolecules
    type(SimcellInfo), intent(in)     :: SimulationCell
    type(CellListInfo), intent(inout) :: CellList
    logical, intent(Out) :: error

    integer :: spc,natoms, TotalNumberOfCells
    integer :: ierror,sys

    TotalNumberOfCells=product(int(SimulationCell%BoxWidth/CellList%MinimumCellWidth))
    if(CellList%TotalNumberOfCells == TotalNumberOfCells)return

    call DeallocateCellListArrays(CellList,error)
    if(error)return

    call InitializeCellListArrays(NumberOfMolecules,SimulationCell,CellList,error)
    if(error)return

  end subroutine ResizeCellListArrays

  subroutine CheckAndIncreaseCellListSize(spc,nmoles,LinkedList,error)
    integer, intent(In)                      :: spc
    integer, dimension(:,:), allocatable, intent(inout)   :: LinkedList
    logical, intent(Out)                     :: error

    integer, dimension(:,:), allocatable :: itemp
    integer :: nmoles,natoms
    integer :: ListSize,NewListSize
    integer :: ierror

    error=.false.

    natoms=Molecule(spc)%NumberOfAtoms
    ListSize=size(LinkedList,1)

    if(nmoles <= ListSize)return

    NewListSize=ListSize+INCR_SIZE

    allocate(itemp(ListSize,natoms),STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to allocate memory to variable itemp. STAT = ',ierror
      error=.true.
      return
    end if
    itemp=LinkedList

    deallocate(LinkedList,STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to deallocate memory from variable LinkedList. STAT = ',ierror
      error=.true.
      return
    end if
    allocate(LinkedList(NewListSize,natoms),STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to allocate memory to variable LinkedList. STAT = ',ierror
      error=.true.
      return
    end if
    LinkedList(1:ListSize,1:natoms)=itemp

    deallocate(itemp,STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to deallocate memory from variable itemp. STAT = ',ierror
      error=.true.
      return
    end if

  end subroutine CheckAndIncreaseCellListSize

  subroutine CreateCellList(Coordinates,SimulationCell,CellList)
    type(SpeciesCoordinates), dimension(:), intent(in) :: Coordinates
    type(SimcellInfo), intent(in)                      :: SimulationCell
    type(CellListInfo), intent(inout)                  :: CellList

    integer :: spc,atm,mol
    integer :: icel,i,sys


    CellList%CellWidth=SimulationCell%BoxWidth/real(CellList%NumberOfCells,PR)
    CellList%CellLength=CellList%CellWidth*sqrt(SimulationCell%RLatVecs(:,1)**2+SimulationCell%RLatVecs(:,2)**2+SimulationCell%RLatVecs(:,3)**2)

    do spc=1,NumberOfSpecies
      CellList%SpeciesCellList(spc)%HeadOfCell=0
      do atm=1,Molecule(spc)%NumberOfAtoms
        do mol=1,Coordinates(spc)%NumberOfMolecules
          icel=GetCellNumber(Coordinates(spc)%Positions(:,mol,atm),SimulationCell,CellList)
          CellList%SpeciesCellList(spc)%LinkedList(mol,atm)=CellList%SpeciesCellList(spc)%HeadOfCell(icel,atm)
          CellList%SpeciesCellList(spc)%HeadOfCell(icel,atm)=mol
        end do
      end do
    end do
  end subroutine CreateCellList

  subroutine UpdateCellList(spc,mol,OldPosition,NewPosition,SimulationCell,CellList)
    integer, intent(in)                  :: spc,mol
    real(PR), dimension(:,:), intent(in) :: OldPosition, NewPosition
    type(SimcellInfo), intent(in)        :: SimulationCell
    type(CellListInfo), intent(inout)    :: CellList

    integer :: i,j,k,atm
    integer :: icelold, icelnew

    i=mol
    do atm=1,Molecule(spc)%NumberOfAtoms
      icelold=GetCellNumber(OldPosition(:,atm),SimulationCell,CellList)
      icelnew=GetCellNumber(NewPosition(:,atm),SimulationCell,CellList)

      if(icelold /= icelnew)then
        j=CellList%SpeciesCellList(spc)%HeadOfCell(icelold,atm)
        if(j==i)then
          CellList%SpeciesCellList(spc)%HeadOfCell(icelold,atm)=CellList%SpeciesCellList(spc)%LinkedList(i,atm)
        else
          k=CellList%SpeciesCellList(spc)%LinkedList(j,atm)
          do while(k /= i)
            j=k
            k=CellList%SpeciesCellList(spc)%LinkedList(j,atm)
          end do
          CellList%SpeciesCellList(spc)%LinkedList(j,atm)=CellList%SpeciesCellList(spc)%LinkedList(i,atm)
        end if
        CellList%SpeciesCellList(spc)%LinkedList(i,atm)=CellList%SpeciesCellList(spc)%HeadOfCell(icelnew,atm)
        CellList%SpeciesCellList(spc)%HeadOfCell(icelnew,atm)=i
      end if
    end do

  end subroutine UpdateCellList

  subroutine InsertMoleculeToCellList(spc,mol,NewPosition,SimulationCell,CellList,error)
    integer, intent(in)                  :: spc,mol
    real(PR), dimension(:,:), intent(in) :: NewPosition
    type(SimcellInfo), intent(in)        :: SimulationCell
    type(CellListInfo), intent(inout)    :: CellList
    logical, intent(Out)                 :: error

    integer :: atm,icel

    error=.false.
    call CheckAndIncreaseCellListSize(spc,mol,CellList%SpeciesCellList(spc)%LinkedList,error)
    if(error)return

    do atm=1,Molecule(spc)%NumberOfAtoms
      icel=GetCellNumber(NewPosition(:,atm),SimulationCell,CellList)
      CellList%SpeciesCellList(spc)%LinkedList(mol,atm)=CellList%SpeciesCellList(spc)%HeadOfCell(icel,atm)
      CellList%SpeciesCellList(spc)%HeadOfCell(icel,atm)=mol
    end do

  end subroutine InsertMoleculeToCellList

  subroutine DeleteMoleculeFromCellList(spc,mol,nmoles,OldPosition,LastPosition,SimulationCell,CellList)
    integer, intent(in)                  :: spc,mol,nmoles
    real(PR), dimension(:,:), intent(in) :: OldPosition,LastPosition
    type(SimcellInfo), intent(in)        :: SimulationCell
    type(CellListInfo), intent(inout)    :: CellList

    integer :: atm,icel
    integer :: j,k

   
    do atm=1,Molecule(spc)%NumberOfAtoms
      !** Remove CurrentMolecule from Cell List
      icel=GetCellNumber(OldPosition(:,atm),SimulationCell,CellList)
      j=CellList%SpeciesCellList(spc)%HeadOfCell(icel,atm)
      if(j == mol)then
        CellList%SpeciesCellList(spc)%HeadOfCell(icel,atm)=CellList%SpeciesCellList(spc)%LinkedList(mol,atm)
      else
        k=CellList%SpeciesCellList(spc)%LinkedList(j,atm)
        do while(k /= mol)
          j=k
          k=CellList%SpeciesCellList(spc)%LinkedList(j,atm)
        end do
        CellList%SpeciesCellList(spc)%LinkedList(j,atm)=CellList%SpeciesCellList(spc)%LinkedList(mol,atm)
      end if

      if(mol == nmoles)cycle
      !** Put the last molecule in place of Current Molecule
      icel=GetCellNumber(LastPosition(:,atm),SimulationCell,CellList)
      j=CellList%SpeciesCellList(spc)%HeadOfCell(icel,atm)
      if(j == nmoles)then
        CellList%SpeciesCellList(spc)%HeadOfCell(icel,atm)=mol
      else
        k=CellList%SpeciesCellList(spc)%LinkedList(j,atm)
        do while(k /= nmoles)
          j=k
          k=CellList%SpeciesCellList(spc)%LinkedList(j,atm)
        end do
        CellList%SpeciesCellList(spc)%LinkedList(j,atm)=mol
      end if
      CellList%SpeciesCellList(spc)%LinkedList(mol,atm)=CellList%SpeciesCellList(spc)%LinkedList(nmoles,atm)
    end do

  end subroutine DeleteMoleculeFromCellList

  function GetCellNumber(rvec,SimulationCell,CellList) result(icel)
    real(PR), dimension(3), intent(IN) :: rvec
    type(SimcellInfo), intent(in)      :: SimulationCell
    type(CellListInfo), intent(in)     :: CellList
    integer :: icel

    integer, dimension(3) :: nvec
    real(PR), dimension(3) :: vec

    vec=rvec
    call ApplyBoundaryCondition(SimulationCell,vec)
    if(SimulationCell%NonOrthorhombic)vec=Matrix_x_Vector(SimulationCell%RLatVecs,vec)

    nvec=min(int((vec+SimulationCell%HalfBoxLength)/CellList%CellLength),CellList%NumberOfCells-1)
    icel=1+nvec(1)
    icel=icel+nvec(2)*CellList%NumberOfCells(1)
    icel=icel+nvec(3)*CellList%NumberOfCells(1)*CellList%NumberOfCells(2)

  end function GetCellNumber
end module cell_list
