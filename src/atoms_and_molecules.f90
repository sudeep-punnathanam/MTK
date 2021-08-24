module atoms_and_molecules
  use consts, only: PR, strlen, lstrlen, Dashed_Line, MAX_NO_OF_ATOMTYPES, MAX_NO_OF_SPECIES
  use utils, only: ErrorMessage, ReadString, FindStringInFile, lowercase, split, OpenFile, CloseFile

  implicit none
  private
  save
  
  !=============================================================================================================
  !** Atom Definition
  type :: AtomInfo
    character(len=strlen)    :: Name=''
    real(PR)                 :: Mass=0._PR
    real(PR)                 :: Charge=0._PR
    logical                  :: has_charge=.false. 
  end type AtomInfo

  !=============================================================================================================
  !** Group Definition
  type :: GroupInfo
    integer                                        :: NumberOfAtoms=0
    integer, dimension(:), allocatable             :: AtomNumber
    real(PR), dimension(:,:), allocatable          :: ReferencePosition
  end type GroupInfo

  !=============================================================================================================
  !** Molecule Definition
  type :: MoleculeInfo
    character(len=strlen)                          :: Name='',Shape='spherical'
    logical                                        :: HasPartialCharges=.false.
    real(PR)                                       :: MolecularWeight=0._PR

    integer                                        :: NumberOfAtoms=0
    integer, dimension(:), allocatable             :: AtomType

    integer                                        :: NumberOfGroups=0
    type(GroupInfo), dimension(:), allocatable     :: Group
  end type MoleculeInfo

  !=============================================================================================================
  integer                                          :: NumberOfAtomTypes=0
  type(AtomInfo), dimension(MAX_NO_OF_ATOMTYPES)   :: Atom

  integer                                          :: NumberOfSpecies=0
  type(MoleculeInfo), dimension(MAX_NO_OF_SPECIES) :: Molecule

  !=============================================================================================================
  public :: Atom, NumberOfAtomTypes
  public :: ReadAtomInfo, GetAtomType
  public :: Molecule, NumberOfSpecies
  public :: ReadMoleculeInfo, GetSpeciesNumber, GetCenterOfMass

contains

  !=============================================================================================================
  subroutine ReadAtomInfo(PROJECT_DIR,LogUnitNo,error)
    character(len=lstrlen), intent(in) :: PROJECT_DIR
    integer, intent(in)                :: LogUnitNo
    logical, intent(out)               :: error

    character(len=lstrlen) :: filename,line
    character(len=strlen), dimension(20) :: fields
    integer :: unitno, nfields
    integer :: atype
    integer :: ierror
    logical :: lopen

    error=.false.
    filename=trim(PROJECT_DIR)//'/Atoms.def'
    call OpenFile(filename,'read','rewind',unitno,error)
    if(error)return

    write(LogUnitNo,'(a)')trim(Dashed_Line)
    write(LogUnitNo,'(a,20x,a)')'#','Atoms'
    write(LogUnitNo,'(a)')trim(Dashed_Line)
    read(unitno,*)
    atype=0
    do
      read(unit=unitno,fmt='(a)',iostat=ierror)line
      if(ierror == -1)exit
      if(len_trim(line) == 0)cycle
      write(LogUnitNo,'(a)')trim(line)

      atype=atype+1
      if(atype > MAX_NO_OF_ATOMTYPES)then
        error=.true.
        write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__,  &
          'Number of atom types exceeds MAX_NO_OF_ATOMTYPES: ',MAX_NO_OF_ATOMTYPES
        return
      end if

      nfields=split(line,fields)
      read(fields(1),*)Atom(atype)%Name
      read(fields(2),*)Atom(atype)%Mass
      if(nfields==3)then
        Atom(atype)%has_charge=.true.
        read(fields(3),*)Atom(atype)%Charge
      end if
    end do
    write(LogUnitNo,*)
    NumberOfAtomTypes=atype

    call CloseFile(unitno,error)
    if(error)return
  end subroutine ReadAtomInfo

  !=============================================================================================================
  subroutine ReadMoleculeInfo(PROJECT_DIR,SimulationInput,error)
    character(len=lstrlen), intent(in)               :: PROJECT_DIR
    character(len=lstrlen), dimension(:), intent(in) :: SimulationInput
    logical, intent(out)                             :: error

    integer :: NumberOfInputLines
    character(len=lstrlen) :: line,string,spcname,filename, atomname
    integer :: l, pos, spc, atm, atype
    integer :: unitno, ierror, lineno, atomno, icount, indx
    logical :: lopen
    real(PR), dimension(3) :: com, vec
    real(PR) :: mass

    error=.false.
    NumberOfInputLines=size(SimulationInput)
    do l=1,NumberOfInputLines
      if(ReadString(SimulationInput(l),'species',string) /= 1)cycle
      pos=ReadString(string,'name',string)
      spcname=adjustl(string)
      spc=GetSpeciesNumber(spcname,error)
      if(spc /= 0)cycle

      NumberOfSpecies=NumberOfSpecies+1
      if(NumberOfSpecies > MAX_NO_OF_SPECIES)then
        error=.true.
        write(ErrorMessage,'(2a,i5,4x,a,i5)')__FILE__,':',__LINE__, &
          'Number of species exceeds MAX_NO_OF_SPECIES: ',MAX_NO_OF_SPECIES
        return
      end if

      spc=NumberOfSpecies
      Molecule(spc)%Name=trim(spcname)
      filename=trim(PROJECT_DIR)//'/'//adjustl(trim(Molecule(spc)%Name))//'.mol'
      call OpenFile(filename,'read','rewind',unitno,error)
      if(error)return

      !** Determine Number of groups and allocate space
      icount=0
      do
        read(unit=unitno,fmt='(a)',iostat=ierror)line
        if(ierror == -1)exit
        line=adjustl(line)
        indx=index(line,'Group')
        if(indx == 1)then
          icount=icount+1
        end if
      end do
      Molecule(spc)%NumberOfGroups=icount
      allocate(Molecule(spc)%Group(icount),STAT=ierror)
      if(ierror /=0)then
        error=.true.
        write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
          'Failed to allocate memory. STAT = ',ierror
        return
      end if

      !** Determine Number of atoms and allocate space
      rewind(unitno)
      icount=0
      do
        read(unit=unitno,fmt='(a)',iostat=ierror)line
        if(ierror == -1)exit
        line=adjustl(line)
        indx=index(line,'Group')
        if(indx == 1)then
          icount=icount+1
          read(line,*)string,Molecule(spc)%Group(icount)%NumberOfAtoms
          Molecule(spc)%NumberOfAtoms=Molecule(spc)%NumberOfAtoms+Molecule(spc)%Group(icount)%NumberOfAtoms

          allocate(Molecule(spc)%Group(icount)%AtomNumber(Molecule(spc)%Group(icount)%NumberOfAtoms), &
            Molecule(spc)%Group(icount)%ReferencePosition(3,Molecule(spc)%Group(icount)%NumberOfAtoms), &
            STAT=ierror)
          if(ierror /=0)then
            error=.true.
            write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
              'Failed to allocate memory. STAT = ',ierror
            return
          end if
        end if
      end do
      allocate(Molecule(spc)%AtomType(Molecule(spc)%NumberOfAtoms),STAT=ierror)
      if(ierror /=0)then
        error=.true.
        write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
          'Failed to allocate memory. STAT = ',ierror
        return
      end if
      
      !** Read Atom Types and Reference Position for Rigid Groups
      rewind(unitno)
      icount=0
      do
        read(unit=unitno,fmt='(a)',iostat=ierror)line
        if(ierror == -1)exit
        line=adjustl(line)
        indx=index(line,'Group')
        if(indx == 1)then
          icount=icount+1
          if(Molecule(spc)%Group(icount)%NumberOfAtoms == 1)then
            read(unitno,*)atomno,atomname
            Molecule(spc)%Group(icount)%AtomNumber(1)=atomno
            Molecule(spc)%Group(icount)%ReferencePosition(:,1)=0._PR
            atype=GetAtomType(atomname,error)
            if(error)return
            Molecule(spc)%AtomType(atomno)=atype
            Molecule(spc)%MolecularWeight=Molecule(spc)%MolecularWeight+Atom(atype)%Mass
            if(Atom(atype)%has_charge)Molecule(spc)%HasPartialCharges=.true.
          else
            read(unitno,'(a)')line
            line=lowercase(line)
            if(trim(adjustl(line)) == 'rigid')then
              com=0._PR
              mass=0._PR
              do atm=1,Molecule(spc)%Group(icount)%NumberOfAtoms
                read(unitno,*)atomno,atomname,vec
                Molecule(spc)%Group(icount)%AtomNumber(atm)=atomno
                Molecule(spc)%Group(icount)%ReferencePosition(:,atm)=vec
                atype=GetAtomType(atomname,error)
                if(error)return
                Molecule(spc)%AtomType(atomno)=atype
                Molecule(spc)%MolecularWeight=Molecule(spc)%MolecularWeight+Atom(atype)%Mass
                if(Atom(atype)%has_charge)Molecule(spc)%HasPartialCharges=.true.

                !** Calculate center of mass of the group
                com=com+Molecule(spc)%Group(icount)%ReferencePosition(:,atm)*Atom(Molecule(spc)%AtomType(atomno))%Mass
                mass=mass+Atom(Molecule(spc)%AtomType(atomno))%Mass
              end do
              !** Adjust Reference positions such that center of mass is at origin
              com=com/mass
              Molecule(spc)%Group(icount)%ReferencePosition=Molecule(spc)%Group(icount)%ReferencePosition- &
                spread(com,2,Molecule(spc)%Group(icount)%NumberOfAtoms)
            else
              !** Only rigid groups implemented
              error=.true.
              write(ErrorMessage,'(2a,i5,4x,a)')__FILE__,':',__LINE__, &
                'Only rigid groups implemented'
              return
            end if
          end if
        end if
      end do

      call CloseFile(unitno,error)
      if(error)return
    end do

  end subroutine ReadMoleculeInfo

  !=============================================================================================================
  function GetAtomType(atomname,error)
    character(len=strlen), intent(In) :: atomname
    logical, intent(Out)              :: error
    integer                           :: GetAtomType

    integer :: atm
    error=.false.
    do atm=1,NumberOfAtomTypes
      if(trim(atomname) == trim(Atom(atm)%Name))then
        GetAtomType=atm
        return
      end if
    end do
    GetAtomType=0
    error=.true.
    write(ErrorMessage,'(2a,i5,4x,3a)')__FILE__,':',__LINE__,'Atom ',trim(atomname),' does not exist.'

  end function GetAtomType

  !=============================================================================================================
  function GetSpeciesNumber(spcname,error)
    character(len=strlen), intent(In) :: spcname
    logical, intent(Out)              :: error
    integer                           :: GetSpeciesNumber

    integer :: spc

    error=.false.
    do spc=1,NumberOfSpecies
      if(trim(spcname) == trim(Molecule(spc)%Name))then
        GetSpeciesNumber=spc
        return
      end if
    end do
    GetSpeciesNumber=0
    error=.true.
    write(ErrorMessage,'(2a,i5,4x,3a)')__FILE__,':',__LINE__,'Species ',trim(spcname),' does not exist'

  end function GetSpeciesNumber

  !=============================================================================
  function GetCenterOfMass(spc,atvecs) result(com)
    integer, intent(IN)                  :: spc
    real(PR), dimension(:,:), intent(IN) :: atvecs
    real(PR), dimension(3)               :: com

    integer :: atm,atype
    real(PR) :: mass

    com=0._PR
    mass=0._PR
    do atm=1,Molecule(spc)%NumberOfAtoms
      atype=Molecule(spc)%AtomType(atm)
      com=com+atvecs(:,atm)*Atom(atype)%Mass
      mass=mass+Atom(atype)%Mass
    end do
    com=com/mass

  end function GetCenterOfMass

end module atoms_and_molecules
