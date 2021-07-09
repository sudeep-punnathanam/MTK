module forcefield
  use consts, only: PR, PI, TWOPI, strlen, lstrlen, Dashed_Line, MAX_NO_OF_ATOMTYPES, K_B, COULOMBIC_CONVERSION_FACTOR, &
    MAX_NO_OF_SPECIES
  use variables, only: REDUCED_UNITS, CutOffDistance, CutOffDistanceSq
  use utils, only: ErrorMessage, GetFileUnit, lowercase, ErrorFunctionComplement, split
  use atoms_and_molecules, only: Atom, NumberOfAtomTypes, GetAtomType, NumberOfSpecies, Molecule
  use ewaldsum, only: Ewald_alpha, Ewald_alpha_by_sqrt_pi, Ewald_kappa, Ewald_kappa_by_two_alpha, PairCoulombInteraction, COULOMB_SCREENING
  use storage, only: InteractionAccumulator
  use random, only: GaussianRandomNumber
  use utils, only: CROSS_PRODUCT, FindStringInFile
  implicit none
  public
  save

  !=============================================================================
  !** No interaction
  integer, parameter :: OFF=0, ON=1
  !** Non-Bonded Two-Body Potentials
  integer, parameter :: HARD_SPHERE=1,LENNARD_JONES=2,TRUNCATED_LENNARD_JONES=3,SHIFTED_LENNARD_JONES=4,YUKAWA=5
  !** Bond Interactions
  integer, parameter :: RIGID_BOND=1,HARMONIC_BOND=2
  !** Angle Interactions
  integer, parameter :: RIGID_ANGLE=1,HARMONIC_ANGLE=2,HARMONIC_COSINE=3
  !** Dihedral Interactions
  integer, parameter :: COSINE_SERIES_1=1

  !=============================================================================
  type :: PotentialModel_info
    character(len=strlen)   :: Name=''
    character(len=strlen)   :: Atom1,Atom2,Atom3,Atom4
    integer                 :: ptype=OFF
    real(PR), dimension(10) :: arg=0._PR
  end type PotentialModel_info
  !=============================================================================

  type(PotentialModel_info), dimension(MAX_NO_OF_ATOMTYPES,MAX_NO_OF_ATOMTYPES) :: PairPotential

  !** Intramolecular forcefield 
  type :: IntraMolecule_Info
    integer :: NumberOfBonds=0, NumberOfAngles=0, NumberOfTorsions=0, NumberOfIntraPairs=0, NumberOfIntraCouls=0
    integer, dimension(:,:), allocatable :: BondAtom, AngleAtom, TorsionAtom, IntraPairAtom, IntraCoulAtom
    type(PotentialModel_info), dimension(:), allocatable :: Bond, Angle, Torsion, IntraPair, IntraCoul
  end type IntraMolecule_Info
  type(IntraMolecule_Info), dimension(MAX_NO_OF_SPECIES) :: Molecule_intra

contains

  !================================================================================
  subroutine ReadPairwiseInteractions(PROJECT_DIR,LogUnitNo,error)
    character(len=lstrlen), intent(in) :: PROJECT_DIR
    integer, intent(in)                :: LogUnitNo
    logical, intent(out)               :: error

    character(len=lstrlen) :: filename,line, potname, rulename
    character(len=strlen) :: AtomName1, AtomName2
    integer :: type1, type2
    integer :: unitno, indx
    integer :: ierror
    logical :: lopen, found

    error=.false.
    filename=trim(PROJECT_DIR)//'/TwoBodyInteractions.def'
    unitno=GetFileUnit(trim(filename),lopen,error)
    if(error)return
    if(lopen)then
      error=.true.
      write(ErrorMessage,'(2a,i5,4x,2a)')__FILE__,':',__LINE__,trim(filename),' already open'
      return
    end if

    open(unit=unitno,file=trim(filename),action='read',iostat=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to open file ,',trim(filename),'. IOSTAT = ',ierror
      error=.true.
      return
    end if

    write(LogUnitNo,'(a)')trim(Dashed_Line)
    write(LogUnitNo,'(a,20x,a)')'#','Pairwise Interactions'
    write(LogUnitNo,'(a)')trim(Dashed_Line)
    !** Like Interactions
    do type1=1,NumberOfAtomTypes
      PairPotential(type1,type1)%Atom1=Atom(type1)%Name
      PairPotential(type1,type1)%Atom2=Atom(type1)%Name
      rewind(unitno)
      read(unitno,*)
      found=.false.
      inner_loop:do
        read(unit=unitno,fmt='(a)',iostat=ierror)line
        if(ierror == -1)exit inner_loop
        if(len_trim(line) == 0)cycle inner_loop
        read(line,*)AtomName1,AtomName2,potname
        if(trim(AtomName1) == Atom(type1)%Name .and. trim(AtomName2) == trim(Atom(type1)%Name))then
          write(LogUnitNo,'(a)')trim(line)
          indx=index(line,trim(potname))
          line=line(indx:)
          call ReadPairPotential(PairPotential(type1,type1),line,error)
          if(error)return
          found=.true.
          exit inner_loop
        end if
      end do inner_loop
      if(.not. found)then
        error=.true.
        write(ErrorMessage,'(2a,i5,4x,5a)')__FILE__,':',__LINE__, &
          'Pair potential between atoms ',trim(AtomName1),' and ',trim(AtomName2),' not found.'
        return
      end if
    end do

    !** Un-Like Interactions
    rewind(unitno)
    read(unitno,*)
    do
      read(unit=unitno,fmt='(a)',iostat=ierror)line
      if(ierror == -1)exit
      if(len_trim(line) == 0)cycle
      read(line,*)AtomName1,AtomName2,rulename
      rulename=adjustl(rulename)
      type1=GetAtomType(AtomName1,error)
      if(error)return
      type2=GetAtomType(AtomName2,error)
      if(error)return
      if(type1 /= type2)then
        write(LogUnitNo,'(a)')trim(line)
        indx=index(line,trim(rulename))
        line=line(indx:)
        PairPotential(type1,type2)%Atom1=Atom(type1)%Name
        PairPotential(type1,type2)%Atom2=Atom(type2)%Name
        call ReadMixingRules(type1,type2,line,error)
        if(error)return
      end if
    end do

    close(unit=unitno,iostat=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
        'Failed to close file ,',trim(filename),'. IOSTAT = ',ierror
      error=.true.
      return
    end if

    !** Coulomb Interactions
    do type1=1,NumberOfAtomTypes
      do type2=type1,NumberOfAtomTypes
        if(Atom(type1)%has_charge .and. Atom(type2)%has_charge)then
          PairCoulombInteraction(type1,type2)=.true.
          PairCoulombInteraction(type2,type1)=.true.
        else
          PairCoulombInteraction(type1,type2)=.false.
          PairCoulombInteraction(type2,type1)=.false.
        end if
      end do
    end do

  end subroutine ReadPairwiseInteractions

  !================================================================================
  ! Read and initialize Mixing Rules
  !================================================================================
  subroutine ReadMixingRules(type1,type2,line,error)
    integer, intent(In)                 :: type1,type2
    character(len=lstrlen), intent(In)  :: line
    logical, intent(Out)                :: error

    real(pr) :: sigma,epsilon,locut,hicut
    character(len=strlen) :: field,rulename

    error=.false.

    read(line,*)field
    rulename=lowercase(field)
    select case (trim(rulename))
    case ('off')
      PairPotential(type1,type2)%Name='OFF'
      PairPotential(type1,type2)%ptype=OFF
    case ('lorentz-berthelot')
      if(PairPotential(type1,type1)%ptype /= PairPotential(type2,type2)%ptype)then
        write(ErrorMessage,'(2a,i5,4x,5a)')__FILE__,':',__LINE__, &
          'Lorentz-Berthelot rule cannot be applied to atoms ',trim(Atom(type1)%Name), &
          ' and ',trim(Atom(type2)%Name),' since their potentials differ.'
        error=.true.
        return
      end if
      PairPotential(type1,type2)%Name=PairPotential(type1,type1)%Name
      PairPotential(type1,type2)%ptype=PairPotential(type1,type1)%ptype
      select case (PairPotential(type1,type2)%ptype)
      case (LENNARD_JONES)
        sigma=(PairPotential(type1,type1)%arg(1)+PairPotential(type2,type2)%arg(1))/2._PR
        epsilon=sqrt(PairPotential(type1,type1)%arg(2)*PairPotential(type2,type2)%arg(2))
        locut=(PairPotential(type1,type1)%arg(3)+PairPotential(type2,type2)%arg(3))/2._PR
        PairPotential(type1,type2)%arg(1)=sigma
        PairPotential(type1,type2)%arg(2)=epsilon
        PairPotential(type1,type2)%arg(3)=locut
        PairPotential(type1,type2)%arg(6)=sigma**2
        PairPotential(type1,type2)%arg(7)=locut**2
      case (TRUNCATED_LENNARD_JONES,SHIFTED_LENNARD_JONES)
        sigma=(PairPotential(type1,type1)%arg(1)+PairPotential(type2,type2)%arg(1))/2._PR
        epsilon=sqrt(PairPotential(type1,type1)%arg(2)*PairPotential(type2,type2)%arg(2))
        locut=(PairPotential(type1,type1)%arg(3)+PairPotential(type2,type2)%arg(3))/2._PR
        hicut=PairPotential(type1,type1)%arg(4)
        PairPotential(type1,type2)%arg(1)=sigma
        PairPotential(type1,type2)%arg(2)=epsilon
        PairPotential(type1,type2)%arg(3)=locut
        PairPotential(type1,type2)%arg(4)=hicut
        if(PairPotential(type1,type2)%ptype == SHIFTED_LENNARD_JONES)then
          PairPotential(type1,type2)%arg(5)=4._PR*epsilon*((sigma/hicut)**12-(sigma/hicut)**6)
        else
          PairPotential(type1,type2)%arg(5)=0._pr
        end if
        PairPotential(type1,type2)%arg(6)=sigma**2
        PairPotential(type1,type2)%arg(7)=locut**2
        PairPotential(type1,type2)%arg(8)=hicut**2
      case Default
        write(ErrorMessage,'(2a,i5,4x,5a)')__FILE__,':',__LINE__,&
          'Lorentz-Berthelot mixing rule for pair potential between atoms ',trim(Atom(type1)%Name),&
          ' and ',trim(Atom(type2)%Name),' not defined.'
        error=.true.
        return
      end select
    case Default
      call ReadPairPotential(PairPotential(type1,type1),line,error)
      if(error)return
    end select
    PairPotential(type2,type1)=PairPotential(type1,type2)

  end subroutine ReadMixingRules

  !================================================================================
  !**  Read intra-molecular potentials from the molecule files
  !================================================================================
  subroutine ReadIntraMolecularPotentials(PROJECT_DIR,LogUnitNo,error)
    character(len=lstrlen), intent(in) :: PROJECT_DIR
    integer, intent(in)                :: LogUnitNo
    logical, intent(out)               :: error

    character(len=lstrlen) :: filename,line, potname, string
    character(len=strlen), dimension(20) :: fields
    integer :: type1, type2, type3, type4
    integer :: unitno, indx, nfields, lineno
    integer :: ierror
    logical :: lopen, found
    integer :: spc, i

    error=.false.
    do spc=1,NumberOfSpecies
      filename=trim(PROJECT_DIR)//'/'//adjustl(trim(Molecule(spc)%Name))//'.mol'
      unitno=GetFileUnit(trim(filename),lopen,error)
      if(error)return
      if(lopen)then
        error=.true.
        write(ErrorMessage,'(2a,i5,4x,3a)')__FILE__,':',__LINE__, &
          'Molecule file ',trim(filename),' already open'
        return
      end if
      open(unit=unitno,file=trim(filename),action='read',iostat=ierror)
      if(ierror /= 0)then
        error=.true.
        write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
          'Failed to open molecule file: ',trim(filename),'. iostat = ',ierror
        return
      end if

      !** Bond Interaction
      call FindStringInFile('Number of Bonds',unitno,lineno,string,error)
      if(error)return
      if(lineno /= 0)then
        read(string,*)Molecule_intra(spc)%NumberOfBonds
        allocate(Molecule_intra(spc)%BondAtom(Molecule_intra(spc)%NumberOfBonds,2), &
                 Molecule_intra(spc)%Bond(Molecule_intra(spc)%NumberOfBonds), &
                 STAT=ierror)
        if(ierror /=0)then
          error=.true.
          write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
            'Failed to allocate memory. STAT = ',ierror
          return
        end if
        do i=1,Molecule_intra(spc)%NumberOfBonds
          read(unitno,'(a)')line
          read(line,*)Molecule_intra(spc)%BondAtom(i,1),Molecule_intra(spc)%BondAtom(i,2),potname
          indx=index(line,trim(potname))
          line=line(indx:)
          type1=Molecule(spc)%AtomType(Molecule_intra(spc)%BondAtom(i,1))
          type2=Molecule(spc)%AtomType(Molecule_intra(spc)%BondAtom(i,2))
          Molecule_intra(spc)%Bond(i)%Atom1=Atom(type1)%Name
          Molecule_intra(spc)%Bond(i)%Atom2=Atom(type2)%Name
          call ReadBondPotential(Molecule_intra(spc)%Bond(i),line,error)
          if(error)return
        end do
      else
        Molecule_intra(spc)%NumberOfBonds=0
      end if

      !** Angle Interaction
      call FindStringInFile('Number of Angles',unitno,lineno,string,error)
      if(error)return
      if(lineno /= 0)then
        read(string,*)Molecule_intra(spc)%NumberOfAngles
        allocate(Molecule_intra(spc)%AngleAtom(Molecule_intra(spc)%NumberOfAngles,3), &
                 Molecule_intra(spc)%Angle(Molecule_intra(spc)%NumberOfAngles), &
                 STAT=ierror)
        if(ierror /=0)then
          error=.true.
          write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
            'Failed to allocate memory. STAT = ',ierror
          return
        end if
        do i=1,Molecule_intra(spc)%NumberOfAngles
          read(unitno,'(a)')line
          read(line,*)Molecule_intra(spc)%AngleAtom(i,1),Molecule_intra(spc)%AngleAtom(i,2), &
            Molecule_intra(spc)%AngleAtom(i,3),potname
          indx=index(line,trim(potname))
          line=line(indx:)
          type1=Molecule(spc)%AtomType(Molecule_intra(spc)%AngleAtom(i,1))
          type2=Molecule(spc)%AtomType(Molecule_intra(spc)%AngleAtom(i,2))
          type3=Molecule(spc)%AtomType(Molecule_intra(spc)%AngleAtom(i,3))
          Molecule_intra(spc)%Bond(i)%Atom1=Atom(type1)%Name
          Molecule_intra(spc)%Bond(i)%Atom2=Atom(type2)%Name
          Molecule_intra(spc)%Bond(i)%Atom3=Atom(type3)%Name
          call ReadAnglePotential(Molecule_intra(spc)%Angle(i),line,error)
          if(error)return
        end do
      else
        Molecule_intra(spc)%NumberOfAngles=0
      end if

      !** Torsion Interaction
      call FindStringInFile('Number of Torsions',unitno,lineno,string,error)
      if(error)return
      if(lineno /= 0)then
        read(string,*)Molecule_intra(spc)%NumberOfTorsions
        allocate(Molecule_intra(spc)%TorsionAtom(Molecule_intra(spc)%NumberOfTorsions,4), &
                 Molecule_intra(spc)%Torsion(Molecule_intra(spc)%NumberOfTorsions), &
                 STAT=ierror)
        if(ierror /=0)then
          error=.true.
          write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
            'Failed to allocate memory. STAT = ',ierror
          return
        end if
        do i=1,Molecule_intra(spc)%NumberOfTorsions
          read(unitno,'(a)')line
          read(line,*)Molecule_intra(spc)%TorsionAtom(i,1),Molecule_intra(spc)%TorsionAtom(i,2), &
            Molecule_intra(spc)%TorsionAtom(i,3),Molecule_intra(spc)%TorsionAtom(i,4),potname
          indx=index(line,trim(potname))
          line=line(indx:)
          type1=Molecule(spc)%AtomType(Molecule_intra(spc)%TorsionAtom(i,1))
          type2=Molecule(spc)%AtomType(Molecule_intra(spc)%TorsionAtom(i,2))
          type3=Molecule(spc)%AtomType(Molecule_intra(spc)%TorsionAtom(i,3))
          type4=Molecule(spc)%AtomType(Molecule_intra(spc)%TorsionAtom(i,4))
          Molecule_intra(spc)%Bond(i)%Atom1=Atom(type1)%Name
          Molecule_intra(spc)%Bond(i)%Atom2=Atom(type2)%Name
          Molecule_intra(spc)%Bond(i)%Atom3=Atom(type3)%Name
          Molecule_intra(spc)%Bond(i)%Atom4=Atom(type4)%Name
          call ReadTorsionPotential(Molecule_intra(spc)%Torsion(i),line,error)
          if(error)return
        end do
      else
        Molecule_intra(spc)%NumberOfTorsions=0
      end if

      !** Intra Pair Interaction
      call FindStringInFile('Number of IntraPairs',unitno,lineno,string,error)
      if(error)return
      if(lineno /= 0)then
        read(string,*)Molecule_intra(spc)%NumberOfIntraPairs
        allocate(Molecule_intra(spc)%IntraPairAtom(Molecule_intra(spc)%NumberOfIntraPairs,2), &
                 Molecule_intra(spc)%IntraPair(Molecule_intra(spc)%NumberOfIntraPairs), &
                 STAT=ierror)
        if(ierror /=0)then
          error=.true.
          write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
            'Failed to allocate memory. STAT = ',ierror
          return
        end if
        do i=1,Molecule_intra(spc)%NumberOfIntraPairs
          read(unitno,'(a)')line
          read(line,*)Molecule_intra(spc)%IntraPairAtom(i,1),Molecule_intra(spc)%IntraPairAtom(i,2),potname
          indx=index(line,trim(potname))
          line=line(indx:)
          type1=Molecule(spc)%AtomType(Molecule_intra(spc)%IntraPairAtom(i,1))
          type2=Molecule(spc)%AtomType(Molecule_intra(spc)%IntraPairAtom(i,2))
          Molecule_intra(spc)%IntraPair(i)%Atom1=Atom(type1)%Name
          Molecule_intra(spc)%IntraPair(i)%Atom2=Atom(type2)%Name
          call ReadPairPotential(Molecule_intra(spc)%IntraPair(i),line,error)
          if(error)return
        end do
      else
        Molecule_intra(spc)%NumberOfIntraPairs=0
      end if

      !** Intra Coul Interaction
      call FindStringInFile('Number of IntraCouls',unitno,lineno,string,error)
      if(error)return
      if(lineno /= 0)then
        read(string,*)Molecule_intra(spc)%NumberOfIntraCouls
        allocate(Molecule_intra(spc)%IntraCoulAtom(Molecule_intra(spc)%NumberOfIntraCouls,2), &
                 Molecule_intra(spc)%IntraCoul(Molecule_intra(spc)%NumberOfIntraCouls), &
                 STAT=ierror)
        if(ierror /=0)then
          error=.true.
          write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
            'Failed to allocate memory. STAT = ',ierror
          return
        end if
        do i=1,Molecule_intra(spc)%NumberOfIntraCouls
          read(unitno,'(a)')line
          read(line,*)Molecule_intra(spc)%IntraCoulAtom(i,1),Molecule_intra(spc)%IntraCoulAtom(i,2), &
            Molecule_intra(spc)%IntraCoul(i)%arg(1)
          type1=Molecule(spc)%AtomType(Molecule_intra(spc)%IntraCoulAtom(i,1))
          type2=Molecule(spc)%AtomType(Molecule_intra(spc)%IntraCoulAtom(i,2))
          Molecule_intra(spc)%IntraCoul(i)%Atom1=Atom(type1)%Name
          Molecule_intra(spc)%IntraCoul(i)%Atom2=Atom(type2)%Name
        end do
      else
        Molecule_intra(spc)%NumberOfIntraCouls=0
      end if

      close(unit=unitno,iostat=ierror)
      if(ierror /= 0)then
        write(ErrorMessage,'(2a,i5,4x,3a,i4)')__FILE__,':',__LINE__, &
          'Failed to close file ,',trim(filename),'. IOSTAT = ',ierror
        error=.true.
        return
      end if
    end do

  end subroutine ReadIntraMolecularPotentials

  !================================================================================
  ! Read and initialize Pair Potential
  !================================================================================
  subroutine ReadPairPotential(PotentialModel,line,error)
    type(PotentialModel_info), intent(inout) :: PotentialModel
    character(len=lstrlen), intent(In)       :: line
    logical, intent(Out)                     :: error

    real(Kind=PR)   :: sigma,epsilon,locut,hicut,kappa
    character(len=strlen) :: potname

    error=.false.

    read(line,*)PotentialModel%Name
    potname=lowercase(PotentialModel%Name)
    select case (trim(potname))
    case ('off')
      PotentialModel%ptype=OFF
    case ('hard-sphere')
      read(line,*)potname,sigma
      PotentialModel%ptype=HARD_SPHERE
      PotentialModel%arg(1)=sigma
      PotentialModel%arg(2)=sigma**2
    case ('lennard-jones')
      read(line,*)potname,sigma,epsilon,locut
      if(.not. REDUCED_UNITS)epsilon=epsilon*K_B
      PotentialModel%ptype=LENNARD_JONES
      PotentialModel%arg(1)=sigma
      PotentialModel%arg(2)=epsilon
      PotentialModel%arg(3)=locut
      PotentialModel%arg(5)=0._PR
      PotentialModel%arg(6)=sigma**2
      PotentialModel%arg(7)=locut**2
    case ('truncated-lennard-jones')
      read(line,*)potname,sigma,epsilon,locut,hicut
      if(.not. REDUCED_UNITS)epsilon=epsilon*K_B
      PotentialModel%ptype=TRUNCATED_LENNARD_JONES
      PotentialModel%arg(1)=sigma
      PotentialModel%arg(2)=epsilon
      PotentialModel%arg(3)=locut
      PotentialModel%arg(4)=hicut
      PotentialModel%arg(5)=0._PR
      PotentialModel%arg(6)=sigma**2
      PotentialModel%arg(7)=locut**2
      PotentialModel%arg(8)=hicut**2
    case ('shifted-lennard-jones')
      read(line,*)potname,sigma,epsilon,locut,hicut
      if(.not. REDUCED_UNITS)epsilon=epsilon*K_B
      PotentialModel%ptype=SHIFTED_LENNARD_JONES
      PotentialModel%arg(1)=sigma
      PotentialModel%arg(2)=epsilon
      PotentialModel%arg(3)=locut
      PotentialModel%arg(4)=hicut
      PotentialModel%arg(5)=4._PR*epsilon*((sigma/hicut)**12-(sigma/hicut)**6)
      PotentialModel%arg(6)=sigma**2
      PotentialModel%arg(7)=locut**2
      PotentialModel%arg(8)=hicut**2
    case ('yukawa')
      read(line,*)potname,sigma,epsilon,kappa
      if(.not. REDUCED_UNITS)epsilon=epsilon*K_B
      PotentialModel%ptype=YUKAWA
      PotentialModel%arg(1)=sigma
      PotentialModel%arg(2)=epsilon
      PotentialModel%arg(3)=kappa
      PotentialModel%arg(4)=sigma**2
    case Default
      write(ErrorMessage,'(2a,i5,4x,6a)')__FILE__,':',__LINE__, &
        'Unknown Pair Potential ',trim(potname),' between atoms ', &
        trim(PotentialModel%Atom1),' and ',trim(PotentialModel%Atom2)
      error=.true.
      return
    end select

  end subroutine ReadPairPotential

  !================================================================================
  ! Read and initialize Bond Potential
  !================================================================================
  subroutine ReadBondPotential(PotentialModel,line,error)
    type(PotentialModel_info), intent(inout) :: PotentialModel
    character(len=lstrlen), intent(In)       :: line
    logical, intent(Out)                     :: error

    real(Kind=PR)   :: r0,k
    character(len=strlen) :: potname

    error=.false.

    read(line,*)PotentialModel%Name
    potname=lowercase(PotentialModel%Name)
    select case (trim(potname))
    case ('rigid-bond')
      read(line,*)potname,r0
      PotentialModel%ptype=RIGID_BOND
      PotentialModel%arg(1)=r0
    case ('harmonic-bond')
      read(line,*)potname,r0,k
      if(.not. REDUCED_UNITS)k=k*K_B
      PotentialModel%ptype=HARMONIC_BOND
      PotentialModel%arg(1)=r0
      PotentialModel%arg(2)=k
    case Default
      write(ErrorMessage,'(2a,i5,4x,6a)')__FILE__,':',__LINE__, &
        'Unknown Bond Potential ',trim(potname),' between atoms ', &
        trim(PotentialModel%Atom1),' and ',trim(PotentialModel%Atom2)
      error=.true.
      return
    end select

  end subroutine ReadBondPotential

  !================================================================================
  ! Read and initialize Angle Potential
  !================================================================================
  subroutine ReadAnglePotential(PotentialModel,line,error)
    type(PotentialModel_info), intent(inout) :: PotentialModel
    character(len=lstrlen), intent(In)       :: line
    logical, intent(Out)                     :: error

    real(Kind=PR)   :: theta0,k
    character(len=strlen) :: potname

    error=.false.

    read(line,*)PotentialModel%Name
    potname=lowercase(PotentialModel%Name)
    select case (trim(potname))
    case ('rigid-angle')
      read(line,*)potname,theta0
      PotentialModel%ptype=RIGID_ANGLE
      PotentialModel%arg(1)=theta0
    case ('harmonic-angle')
      read(line,*)potname,theta0,k
      if(.not. REDUCED_UNITS)k=k*K_B
      theta0=theta0*PI/180._PR
      PotentialModel%ptype=HARMONIC_ANGLE
      PotentialModel%arg(1)=theta0
      PotentialModel%arg(2)=k
    case ('harmonic-cosine')
      read(line,*)potname,theta0,k
      if(.not. REDUCED_UNITS)k=k*K_B
      theta0=theta0*PI/180._PR
      PotentialModel%ptype=HARMONIC_COSINE
      PotentialModel%arg(1)=cos(theta0)
      PotentialModel%arg(2)=k
    case Default
      write(ErrorMessage,'(2a,i5,4x,8a)')__FILE__,':',__LINE__, &
        'Unknown Angle Potential ',trim(potname),' between atoms ', &
        trim(PotentialModel%Atom1),', ',trim(PotentialModel%Atom2),' and ',trim(PotentialModel%Atom3)
      error=.true.
      return
    end select

  end subroutine ReadAnglePotential

  !================================================================================
  ! Read and initialize Torsion Potential
  !================================================================================
  subroutine ReadTorsionPotential(PotentialModel,line,error)
    type(PotentialModel_info), intent(inout) :: PotentialModel
    character(len=lstrlen), intent(In)       :: line
    logical, intent(Out)                     :: error

    real(Kind=PR), dimension(3)   :: c
    character(len=strlen) :: potname

    error=.false.

    read(line,*)PotentialModel%Name
    potname=lowercase(PotentialModel%Name)
    select case (trim(potname))
    case ('cosine-series-1')
      read(line,*)potname,c
      if(.not. REDUCED_UNITS)c=c*K_B
      PotentialModel%ptype=COSINE_SERIES_1
      PotentialModel%arg(1:3)=c
    case Default
      write(ErrorMessage,'(2a,i5,4x,10a)')__FILE__,':',__LINE__, &
        'Unknown Torsion Potential ',trim(potname),' between atoms ', &
        trim(PotentialModel%Atom1),', ',trim(PotentialModel%Atom2), &
        ', ',trim(PotentialModel%Atom3),' and ',trim(PotentialModel%Atom4)
      error=.true.
      return
    end select

  end subroutine ReadTorsionPotential

  !================================================================================
  !      van der Waals interactions
  !================================================================================
  pure function TwoBodyNonBondedInteraction(TwoBodyNonBonded,rijsq,vfactor) result(Accumulator)
    type(PotentialModel_info), intent(in)     :: TwoBodyNonBonded
    real(PR), dimension(:), intent(in)        :: rijsq,vfactor
    type(InteractionAccumulator)              :: Accumulator

    integer :: n,i
    real(PR) :: atmvirial

    real(PR) :: sigmasq,epsilon,locutsq,ucut, hicutsq
    real(PR) :: sr2,sr6,sr12
    real(PR) :: sigma,kappa,rij,nrg

    n=size(rijsq)
    select case (TwoBodyNonBonded%ptype)
    case (OFF)
    case (HARD_SPHERE)
      do i=1,n
        if(rijsq(i) < TwoBodyNonBonded%arg(2))then
          Accumulator%overlap=.true.
          return
        end if
      end do
    case (LENNARD_JONES)
      sigmasq=TwoBodyNonBonded%arg(6)
      epsilon=TwoBodyNonBonded%arg(2)
      locutsq=TwoBodyNonBonded%arg(7)
      do i=1,n
        if(rijsq(i) < locutsq)then
          Accumulator%overlap=.true.
          return
        else if(rijsq(i) < CutOffDistanceSq)then
          sr2=sigmasq/rijsq(i)
          sr6=sr2**3
          sr12=sr6**2
          Accumulator%Energy=Accumulator%Energy+4._PR*epsilon*(sr12-sr6)
          atmvirial=8._PR*epsilon*(2._PR*sr12-sr6)
          Accumulator%Virial=Accumulator%Virial+vfactor(i)*atmvirial
        end if
      end do
    case (TRUNCATED_LENNARD_JONES,SHIFTED_LENNARD_JONES)
      sigmasq=TwoBodyNonBonded%arg(6)
      epsilon=TwoBodyNonBonded%arg(2)
      locutsq=TwoBodyNonBonded%arg(7)
      hicutsq=TwoBodyNonBonded%arg(8)
      ucut=TwoBodyNonBonded%arg(5)
      do i=1,n
        if(rijsq(i) < locutsq)then
          Accumulator%overlap=.true.
          return
        else if(rijsq(i) < hicutsq)then
          sr2=sigmasq/rijsq(i)
          sr6=sr2**3
          sr12=sr6**2
          Accumulator%Energy=Accumulator%Energy+4._PR*epsilon*(sr12-sr6)-ucut
          atmvirial=8._PR*epsilon*(2._PR*sr12-sr6)
          Accumulator%Virial=Accumulator%Virial+vfactor(i)*atmvirial
        end if
      end do
    case (YUKAWA)
      sigma=TwoBodyNonBonded%arg(1)
      epsilon=TwoBodyNonBonded%arg(2)
      kappa=TwoBodyNonBonded%arg(3)
      locutsq=TwoBodyNonBonded%arg(4)
      do i=1,n
        if(rijsq(i) < locutsq)then
          Accumulator%overlap=.true.
          return
        else if(rijsq(i) < CutOffDistanceSq)then
          rij=sqrt(rijsq(i))
          nrg=epsilon*(sigma/rij)*exp(kappa*(sigma-rij))
          Accumulator%Energy=Accumulator%Energy+nrg
          atmvirial=nrg*(1._PR+kappa*rij)/3._PR
          Accumulator%Virial=Accumulator%Virial+vfactor(i)*atmvirial
          Accumulator%dUdKappa=Accumulator%dUdKappa+nrg*(1._PR-rij/sigma)
        end if
      end do
    end select
  end function TwoBodyNonBondedInteraction

  !================================================================================
  !      Ewald Short Range Interactions
  !================================================================================
  pure function EwaldRealSpaceInteraction(type1,type2,rijsq,vfactor) result(Accumulator)
    integer, intent(IN)                      :: type1,type2
    real(PR), dimension(:), intent(in)       :: rijsq,vfactor
    type(InteractionAccumulator)             :: Accumulator

    integer :: i,n
    real(PR) :: charge1, charge2
    real(PR) :: rij,erfc_kr
    real(PR) :: A,B,C,D,E,F
    real(PR) :: atmvirial

    n=size(rijsq)
    if(PairCoulombInteraction(type1,type2))then
      charge1=Atom(type1)%Charge
      charge2=Atom(type2)%Charge
      if(COULOMB_SCREENING)then
        do i=1,n
          if (rijsq(i) < CutOffDistanceSq)then
            rij=sqrt(rijsq(i))
            A=Ewald_alpha*rij+Ewald_kappa_by_two_alpha
            B=Ewald_alpha*rij-Ewald_kappa_by_two_alpha
            C=Ewald_kappa*rij
            D=ErrorFunctionComplement(A)
            E=ErrorFunctionComplement(B)
            F=exp(C)
            Accumulator%Energy=Accumulator%Energy+(D*F+E/F)/rij/2._PR
            atmvirial=(D*F*(1._PR-C)+E/F*(1._PR+C))/rij/2._PR+ Ewald_alpha_by_sqrt_pi*(exp(-A**2)*F+exp(-B**2)/F)
            Accumulator%Virial=Accumulator%Virial+atmvirial*vfactor(i)
          end if
        end do
      else
        do i=1,n
          if (rijsq(i) < CutOffDistanceSq)then
            rij=sqrt(rijsq(i))
            erfc_kr=ErrorFunctionComplement(Ewald_alpha*rij)
            Accumulator%Energy=Accumulator%Energy+erfc_kr/rij
            atmvirial=(erfc_kr/rij +2._PR*Ewald_alpha_by_sqrt_pi*exp(-(Ewald_alpha*rij)**2))
            Accumulator%Virial=Accumulator%Virial+atmvirial*vfactor(i)
          end if
        end do
      end if
      Accumulator%Energy=Accumulator%Energy*charge1*charge2
      Accumulator%Virial=Accumulator%Virial*charge1*charge2/3._PR
    end if

    !** Convert energy units
    if(.not. REDUCED_UNITS)then
      Accumulator%Energy=Accumulator%Energy*COULOMBIC_CONVERSION_FACTOR
      Accumulator%Virial=Accumulator%Virial*COULOMBIC_CONVERSION_FACTOR
    end if

  end function EwaldRealSpaceInteraction

  !================================================================================
  !      Tail corrections for van der Waals interactions
  !================================================================================
  pure function TailCorrections(TwoBodyNonBonded) result(Accumulator)
    type(PotentialModel_info), intent(in)     :: TwoBodyNonBonded
    type(InteractionAccumulator)              :: Accumulator

    real(PR) :: sigma, epsilon, hicut, sr3, sr9
    real(PR) :: kappa, factor1, factor2

    select case (TwoBodyNonBonded%ptype)
    case (OFF,TRUNCATED_LENNARD_JONES,SHIFTED_LENNARD_JONES)
    case (LENNARD_JONES)
      sigma=TwoBodyNonBonded%arg(1)
      epsilon=TwoBodyNonBonded%arg(2)
      hicut=CutOffDistance
      sr3=(sigma/hicut)**3
      sr9=sr3**3
      Accumulator%Energy=8._PR/3._PR*PI*(sigma**3)*epsilon*(sr9/3._PR-sr3)
      Accumulator%Virial=16._PR/3._PR*PI*(sigma**3)*epsilon*(2._PR*sr9/3._PR-sr3)
    case (YUKAWA)
      sigma=TwoBodyNonBonded%arg(1)
      epsilon=TwoBodyNonBonded%arg(2)
      kappa=TwoBodyNonBonded%arg(3)
      hicut=CutOffDistance
      Accumulator%Energy=TWOPI*epsilon*sigma/kappa**2*exp(kappa*(sigma-hicut))*(1._PR+kappa*hicut)
      Accumulator%Virial=TWOPI*epsilon*sigma/kappa**2*exp(kappa*(sigma-hicut))*(1._PR+kappa*hicut+(kappa*hicut)**2/3._PR)
      Accumulator%dUdKappa=TWOPI*epsilon/kappa**3*exp(kappa*(sigma-hicut))*(sigma*kappa+sigma*hicut*kappa**2-(kappa*hicut)**2-2._pr*kappa*hicut-2._PR)
    end select
  end function TailCorrections

  !====================================================================================
  !         Bond interaction
  !====================================================================================
  pure function BondInteraction(Bond,r) result(energy)
    type(PotentialModel_info), intent(in)   :: Bond
    real(PR), intent(in)                    :: r
    real(PR)                                :: energy

    real(PR) :: k,r0

    select case (Bond%ptype)
    case (OFF)
      energy=0._PR
    case (RIGID_BOND)
      energy=0._PR
    case (HARMONIC_BOND)
      r0=Bond%arg(1)
      k=Bond%arg(2)
      energy=0.5_PR*k*(r-r0)**2
    end select
  end function BondInteraction

  !====================================================================================
  !         Sample Bond Length for CBMC move
  !====================================================================================
  function GetBondLength(Bond,beta) result(BondLength)
    type(PotentialModel_info), intent(In)   :: Bond
    real(PR), intent(In)                    :: beta
    real(PR)                                :: BondLength

    real(PR) :: k,r0,sigma

    select case (Bond%ptype)
    case (RIGID_BOND)
      BondLength=Bond%arg(1)
    case (HARMONIC_BOND)
      r0=Bond%arg(1)
      k=Bond%arg(2)
      sigma=1._PR/SQRT(k*beta)
      BondLength=r0+sigma*GaussianRandomNumber()
    end select
  end function GetBondLength

  !====================================================================================
  !         Get Rosenbluth Weight for Bond Generation
  !====================================================================================
  function GetRosenbluthWeightForBond(Bond,beta,BondLength) result(RosenbluthWeight)
    type(PotentialModel_info), intent(in)   :: Bond
    real(PR), intent(in)                    :: beta
    real(PR), intent(in)                    :: BondLength
    real(PR)                                :: RosenbluthWeight

    select case (Bond%ptype)
    case (RIGID_BOND)
      RosenbluthWeight=1._PR
    case (HARMONIC_BOND)
      RosenbluthWeight=BondLength**2
    end select
  end function GetRosenbluthWeightForBond

  !====================================================================================
  !         Angle interaction
  !====================================================================================
  pure function AngleInteraction(Angle,theta) result(energy)
    type(PotentialModel_info), intent(In)    :: Angle
    real(PR), intent(In)                     :: theta
    real(PR)                                 :: energy

    real(PR) :: k,theta0,costheta0

    select case (Angle%ptype)
    case (OFF)
      energy=0._PR
    case (RIGID_ANGLE)
      energy=0._PR
    case (HARMONIC_ANGLE)
      theta0=Angle%arg(1)
      k=Angle%arg(2)
      energy=0.5_PR*k*(theta-theta0)**2
    case (HARMONIC_COSINE)
      costheta0=Angle%arg(1)
      k=Angle%arg(2)
      energy=0.5_PR*k*(cos(theta)-costheta0)**2
    end select
  end function AngleInteraction

  !====================================================================================
  !         Get Bend Angle for CBMC move
  !====================================================================================
  function GetBendAngle(Angle,beta) result(BendAngle)
    type(PotentialModel_info), intent(in)    :: Angle
    real(PR), intent(in)                     :: beta
    real(PR)                                 :: BendAngle

    real(PR) :: k,theta0,costheta0,sigma

    select case (Angle%ptype)
    case (RIGID_ANGLE)
      BendAngle=Angle%arg(1)
    case (HARMONIC_ANGLE)
      theta0=Angle%arg(1)
      k=Angle%arg(2)
      sigma=1._PR/SQRT(k*beta)
      BendAngle=theta0+sigma*GaussianRandomNumber()
    case (HARMONIC_COSINE)
      costheta0=Angle%arg(1)
      k=Angle%arg(2)
      sigma=1._PR/SQRT(k*beta)
      BendAngle=acos(costheta0+sigma*GaussianRandomNumber())
    end select
  end function GetBendAngle

  !====================================================================================
  !         Get Rosenbluth Weight for Angle Generation
  !====================================================================================
  function GetRosenbluthWeightForAngle(Angle,beta,BendAngle) result(RosenbluthWeight)
    type(PotentialModel_info), intent(in)   :: Angle
    real(PR), intent(in)                    :: beta
    real(PR), intent(in)                    :: BendAngle
    real(PR)                                :: RosenbluthWeight

    select case (Angle%ptype)
    case (RIGID_ANGLE)
      RosenbluthWeight=1._PR
    case (HARMONIC_ANGLE)
      RosenbluthWeight=sin(BendAngle)
    case (HARMONIC_COSINE)
      RosenbluthWeight=1._PR
    end select

  end function GetRosenbluthWeightForAngle

  !====================================================================================
  !         Dihedral interaction
  !====================================================================================
  pure function TorsionInteraction(Dihedral,phi) result(energy)
    type(PotentialModel_info), intent(In)    :: Dihedral
    real(PR), intent(In)                     :: phi
    real(PR)                                 :: energy

    real(PR) :: cosphi, cos2phi, cos3phi

    select case (Dihedral%ptype)
    case (COSINE_SERIES_1)
      cosphi=cos(phi)
      cos2phi=2._PR*cosphi**2-1._PR
      cos3phi=4._PR*cosphi**3-3._PR*cosphi
      energy=Dihedral%arg(1)*(1._PR+cosphi)+Dihedral%arg(2)*(1._PR-cos2phi)+Dihedral%arg(3)*(1._PR+cos3phi)
    end select
  end function TorsionInteraction

  !================================================================================
  !      Intra Pair interactions
  !================================================================================
  function IntraPairInteraction(IntraPair,rijsq,overlap) result(energy)
    type(PotentialModel_info), intent(in)    :: IntraPair
    real(PR), intent(in)                     :: rijsq
    logical, intent(out)                     :: overlap
    real(PR)                                 :: energy

    real(PR) :: sigmasq,epsilon,locutsq,ucut,hicutsq
    real(PR) :: sr2,sr6,sr12
    real(PR) :: sigma,kappa,rij

    overlap=.false.
    select case (IntraPair%ptype)
    case (HARD_SPHERE)
      energy=0._PR
      if(rijsq < IntraPair%arg(2))overlap=.true.
    case (LENNARD_JONES)
      sigmasq=IntraPair%arg(6)
      epsilon=IntraPair%arg(2)
      locutsq=IntraPair%arg(7)
      if(rijsq < locutsq)then
        overlap=.true.
      else
        sr2=sigmasq/rijsq
        sr6=sr2**3
        sr12=sr6**2
        energy=4._PR*epsilon*(sr12-sr6)
      end if
    case (TRUNCATED_LENNARD_JONES,SHIFTED_LENNARD_JONES)
      sigmasq=IntraPair%arg(6)
      epsilon=IntraPair%arg(2)
      locutsq=IntraPair%arg(7)
      hicutsq=IntraPair%arg(8)
      ucut=IntraPair%arg(5)
      energy=0._PR
      if(rijsq < locutsq)then
        overlap=.true.
      else if(rijsq < hicutsq)then
        sr2=sigmasq/rijsq
        sr6=sr2**3
        sr12=sr6**2
        energy=4._PR*epsilon*(sr12-sr6)-ucut
      end if
    case (YUKAWA)
      sigma=IntraPair%arg(1)
      epsilon=IntraPair%arg(2)
      kappa=IntraPair%arg(3)
      locutsq=IntraPair%arg(4)
      if(rijsq < locutsq)then
        overlap=.true.
      else
        rij=sqrt(rijsq)
        energy=epsilon*(sigma/rij)*exp(kappa*(sigma-rij))
      end if
    end select
  end function IntraPairInteraction

  !================================================================================
  !      Ewald Short Range Interactions
  !================================================================================
  pure function IntraCoulInteraction(IntraCoul,rijsq) result(energy)
    type(PotentialModel_info), intent(in)    :: IntraCoul
    real(PR), intent(in)                     :: rijsq
    real(PR)                                 :: energy

    real(PR) :: rij

    rij=sqrt(rijsq)
    energy=IntraCoul%arg(1)/rij

    !** Convert energy units
    if(.not. REDUCED_UNITS)then
      energy=energy*COULOMBIC_CONVERSION_FACTOR
    end if

  end function IntraCoulInteraction

  !================================================================================
end module forcefield

