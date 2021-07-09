module ewaldsum
  use consts, only: PR, PI, TWOPI, MAX_NO_OF_SYSTEMS, MAX_NO_OF_ATOMS, MAX_NO_OF_ATOMTYPES, lstrlen, COULOMBIC_CONVERSION_FACTOR
  use config, only: SpeciesCoordinates
  use atoms_and_molecules, only: Molecule, NumberOfSpecies, Atom, NumberOfAtomTypes
  use simcell, only: SimcellInfo
  use storage, only: StorageInteractions
  use utils, only: ErrorMessage, Matrix_x_Vector, Vector_x_Matrix, ReadStringFromFile, ErrorFunctionComplement
  implicit none
  private
  save

  !=============================================================================================================================
  !** Constants
  integer, parameter :: MAX_NO_OF_KVECTORS=20000, MAX_KMAX=15
  !** Variables
  logical  :: REDUCED_UNITS
  logical   :: TIN_FOIL_BC=.false.
  integer :: EwaldPrecision=5
  logical :: COULOMB_INTERACTION=.false., COULOMB_SCREENING=.false.
  logical, dimension(MAX_NO_OF_ATOMTYPES,MAX_NO_OF_ATOMTYPES) :: PairCoulombInteraction=.false.

  !** Width of the Gaussian
  real(PR) :: Ewald_alpha, Ewald_alpha_by_sqrt_pi

  !** Screened Coulomb Interactions
  real(PR) :: Ewald_kappa=0._PR,Ewald_kappasq=0._PR, Ewald_kappa_by_two_alpha=0._PR, Ewald_kappasq_by_four_alphasq=0._PR

  !** KVectors information
  integer, dimension(MAX_NO_OF_SYSTEMS) :: NumberOfKVectors
  integer, dimension(3,MAX_NO_OF_SYSTEMS) :: kmax
  type :: FourierMaskInfo
    logical, dimension(:,:,:), allocatable :: mask
  end type FourierMaskInfo
  type(FourierMaskInfo), dimension(MAX_NO_OF_SYSTEMS) :: FourierMask

  !** Ewald Kfactor, Structure Factor and Dipoles
  type :: KfactorInfo
    real(PR), dimension(MAX_NO_OF_KVECTORS) :: Kfactor=0._PR, KfactorV=0._PR
  end type KfactorInfo

  type :: StructureFactorInfo
    complex(PR), dimension(MAX_NO_OF_KVECTORS) :: RhoK=0._PR, RhokV=0._PR
    real(PR), dimension(3)                     :: Dipole=0._PR, DipoleV=0._PR
  end type StructureFactorInfo

  !** Temporary Vectors
  integer :: nvecs
  real(PR), dimension(3,MAX_NO_OF_ATOMS)  :: atvecs
  real(PR), dimension(3,MAX_NO_OF_ATOMS)  :: comvecs
  real(PR), dimension(MAX_NO_OF_ATOMS)    :: charges
  real(PR), dimension(MAX_NO_OF_ATOMS,0:MAX_KMAX) :: kxrx2
  real(PR), dimension(MAX_NO_OF_ATOMS,-MAX_KMAX:MAX_KMAX) :: kyry2,kzrz2
  complex(PR), dimension(MAX_NO_OF_ATOMS,0:MAX_KMAX) :: eikx
  complex(PR), dimension(MAX_NO_OF_ATOMS,-MAX_KMAX:MAX_KMAX) :: eiky,eikz
  complex(PR), dimension(MAX_NO_OF_ATOMS) :: eikxky,eikr

  public :: Ewald_alpha, Ewald_alpha_by_sqrt_pi, Ewald_kappa, Ewald_kappa_by_two_alpha, PairCoulombInteraction
  public :: NumberOfKVectors
  public :: COULOMB_INTERACTION, COULOMB_SCREENING
  public :: KfactorInfo, StructureFactorInfo
  public :: ReadAndConfigureEwaldSettings
  public :: CalculateKfactors, CalculateStructureFactor
  public :: CalculateTotalStructureFactor, CalculateMoleculeStructureFactor
  public :: EwaldFourierInteraction, TotalIntraMolecularFourierInteraction, IntramolecularFourier

contains

  !=============================================================================================================================
  !** Parameters chosen from guidelines given in David Fincham, Molecular Simulation, v 13, p 1-9, (1994)
  !=============================================================================================================================
  subroutine ReadAndConfigureEwaldSettings(SimulationInput,SimulationCell,NumberOfSystems,error)
    character(len=lstrlen), dimension(:), intent(in) :: SimulationInput
    type(SimcellInfo), dimension(:), intent(in) :: SimulationCell
    integer, intent(in) :: NumberOfSystems
    logical, intent(Out)   :: error

    integer                :: type1, type2, lineno, sys
    integer                :: spc, atm1
    integer                :: ierror
    integer                :: nk,kx,ky,kz,ksq
    real(PR)               :: rksqmax, rkx,rky,rkz,rksq
    real(PR)               :: prec, charge1, charge2
    real(PR)               :: CutOffDistance
    real(PR), dimension(3) :: boxlength
    character(len=lstrlen) :: line

    error=.false.

    call ReadStringFromFile(line,'reduced units',SimulationInput,'SimulationInput',lineno)
    if(lineno == 0)then
      REDUCED_UNITS=.false.
    else
      REDUCED_UNITS=.true.
    end if

    call ReadStringFromFile(line,'potential cutoff distance',SimulationInput,'SimulationInput',lineno,error=error)
    if(error)return
    read(line,*)CutOffDistance

    COULOMB_INTERACTION=any(Molecule%HasPartialCharges)
    if(.not. COULOMB_INTERACTION)return

    call ReadStringFromFile(line,'ewald precision',SimulationInput,'Simulation Input',lineno)
    if(lineno /= 0)read(line,*)EwaldPrecision
    prec=real(EwaldPrecision,PR)*log(10._pr)

    call ReadStringFromFile(line,'coulomb screening',SimulationInput,'Simulation Input',lineno)
    if(lineno /= 0)then
      COULOMB_SCREENING=.true.
      read(line,*)Ewald_kappa
      Ewald_kappasq=Ewald_kappa**2
    end if

    call ReadStringFromFile(line,'tin foil bc',SimulationInput,'Simulation Input',lineno)
    if(lineno /= 0)TIN_FOIL_BC=.true.

    !** Determine alpha
    Ewald_alpha=sqrt(prec)/CutOffDistance
    Ewald_alpha_by_sqrt_pi=Ewald_alpha/sqrt(PI)
    Ewald_kappa_by_two_alpha=Ewald_kappa/2._PR/Ewald_alpha
    Ewald_kappasq_by_four_alphasq=Ewald_kappa_by_two_alpha**2

    !** Determine kmax and the number of k vectors
    do sys=1,NumberOfSystems
      boxlength=SimulationCell(sys)%BoxLength
      kmax(:,sys)=int(prec*boxlength/CutOffDistance/PI)+1
      if(any(kmax(:,sys) > MAX_KMAX))then
        error=.true.
        write(ErrorMessage,'(2a,i5,4x,a,3i4,a,i4)')__FILE__,':',__LINE__, &
          'Values of kmax(',kmax(1:3,sys),') exceeds MAX_KMAX ', MAX_KMAX
        return
      end if
    
      !** allocate memory
      allocate(FourierMask(sys)%mask(0:kmax(1,sys),-kmax(2,sys):kmax(2,sys),-kmax(3,sys):kmax(3,sys)),STAT=ierror)
      if(ierror /= 0)then
        error=.true.
        write(ErrorMessage,'(2a,i5,4x,a,i5)')__FILE__,':',__LINE__, &
          'Memory Allocation Failure. ierror = ',ierror
        return
      end if
      FourierMask(sys)%mask=.true.

      rksqmax=maxval(TWOPI*(real(kmax(:,sys),PR)/SimulationCell(sys)%BoxWidth))**2
      nk=0
      do kx=0,kmax(1,sys)
        rkx=TWOPI*real(kx,PR)/SimulationCell(sys)%BoxLength(1)
        do ky=-kmax(2,sys),kmax(2,sys)
        rky=TWOPI*real(ky,PR)/SimulationCell(sys)%BoxLength(2)
          do kz=-kmax(3,sys),kmax(3,sys)
          rkz=TWOPI*real(kz,PR)/SimulationCell(sys)%BoxLength(3)
            ksq=kx**2+ky**2+kz**2
            rksq=rkx**2+rky**2+rkz**2
            if(rksq > rksqmax)cycle
            if (ksq == 0 .and. (.not. COULOMB_SCREENING))cycle
            FourierMask(sys)%mask(kx,ky,kz)=.false.
            nk = nk + 1
          end do
        end do
      end do
      NumberOfKVectors(sys)=nk
    end do

  end subroutine ReadAndConfigureEwaldSettings

  !=============================================================================================================================
  subroutine CalculateKfactors(sys,SimulationCell,Kfactor)
    integer, intent(in)              :: sys
    type(SimcellInfo), intent(in)    :: SimulationCell
    type(KfactorInfo), intent(inout) :: Kfactor
    
    integer :: nk,kx,ky,kz
    real(PR) :: boxlength(3),kvec(3)
    real(PR) :: rksq, b, prefactor
    integer :: ierror

    b=1._PR/4._PR/Ewald_alpha**2
    boxlength=SimulationCell%BoxLength
    nk=0
    do kx=0,kmax(1,sys)
      if(kx == 0)then
        prefactor=1._PR
      else
        prefactor=2._PR
      end if
      kvec(1)=TWOPI*real(kx,PR)/boxlength(1)
      do ky=-kmax(2,sys),kmax(2,sys)
        kvec(2)=TWOPI*real(ky,PR)/boxlength(2)
        do kz=-kmax(3,sys),kmax(3,sys)
          if(FourierMask(sys)%mask(kx,ky,kz))cycle
          nk = nk + 1
          kvec(3)=TWOPI*real(kz,PR)/boxlength(3)
          if(SimulationCell%NonOrthorhombic)kvec=Vector_x_Matrix(kvec,SimulationCell%RLatVecs)
          rksq=kvec(1)**2+kvec(2)**2+kvec(3)**2
          Kfactor%Kfactor(nk)=prefactor*TWOPI*exp(-(rksq+Ewald_kappasq)*b)/(rksq+Ewald_kappasq)/SimulationCell%Volume
          Kfactor%KfactorV(nk)=Kfactor%Kfactor(nk)*((rksq+3._pr*Ewald_kappasq)/(rksq+Ewald_kappasq)-2._PR*rksq*b)
        end do
      end do
    end do

  end subroutine CalculateKfactors

  !=============================================================================================================================
  subroutine CalculateTotalStructureFactor(sys,Coordinates,SimulationCell,StructureFactor,error)
    integer, intent(in)                                :: sys
    type(SpeciesCoordinates), dimension(:), intent(in) :: Coordinates
    type(SimcellInfo), intent(in)                      :: SimulationCell
    type(StructureFactorInfo), intent(inout)           :: StructureFactor
    logical, intent(out)                               :: error

    integer :: spc,mol,atm,atype

    error=.false.

    !** Prepare Temporary Arrays
    nvecs=0
    do spc=1,NumberOfSpecies
      do atm=1,Molecule(spc)%NumberOfAtoms
        atype=Molecule(spc)%AtomType(atm)
        if(Atom(atype)%has_charge)then
          do mol=1,Coordinates(spc)%NumberOfMolecules
            nvecs=nvecs+1
            atvecs(1,nvecs)=Coordinates(spc)%Positions(1,mol,atm)
            atvecs(2,nvecs)=Coordinates(spc)%Positions(2,mol,atm)
            atvecs(3,nvecs)=Coordinates(spc)%Positions(3,mol,atm)
            comvecs(1,nvecs)=Coordinates(spc)%Positions(1,mol,atm)-Coordinates(spc)%CenterOfMass(1,mol)
            comvecs(2,nvecs)=Coordinates(spc)%Positions(2,mol,atm)-Coordinates(spc)%CenterOfMass(2,mol)
            comvecs(3,nvecs)=Coordinates(spc)%Positions(3,mol,atm)-Coordinates(spc)%CenterOfMass(3,mol)
            charges(nvecs)=Atom(atype)%Charge
          end do
        end if
      end do
    end do
    if(nvecs > MAX_NO_OF_ATOMS)then
      write(ErrorMessage,'(2a,i5,4x,a,i10)')__FILE__,':',__LINE__, &
        'Number of atomic vectors exceeds MAX_NO_OF_ATOMS:',MAX_NO_OF_ATOMS
      error=.true.
      return
    end if
    if(nvecs == 0)return

    call CalculateStructureFactor(sys,SimulationCell,StructureFactor)
  end subroutine CalculateTotalStructureFactor

  !=============================================================================================================================
  subroutine CalculateMoleculeStructureFactor(sys,spc,MoleculePosition,CenterOfMass,SimulationCell,StructureFactor)
    integer, intent(in)                       :: sys,spc
    real(PR), dimension(:,:), intent(in)      :: MoleculePosition
    real(PR), dimension(3), intent(in)        :: CenterOfMass
    type(SimcellInfo), intent(in)             :: SimulationCell
    type(StructureFactorInfo), intent(inout)  :: StructureFactor

    integer :: atm,atype

    !** Prepare Temporary Arrays
    nvecs=0
    do atm=1,Molecule(spc)%NumberOfAtoms
      atype=Molecule(spc)%AtomType(atm)
      if(Atom(atype)%has_charge)then
        nvecs=nvecs+1
        atvecs(1,nvecs)=MoleculePosition(1,atm)
        atvecs(2,nvecs)=MoleculePosition(2,atm)
        atvecs(3,nvecs)=MoleculePosition(3,atm)
        comvecs(1,nvecs)=MoleculePosition(1,atm)-CenterOfMass(1)
        comvecs(2,nvecs)=MoleculePosition(2,atm)-CenterOfMass(2)
        comvecs(3,nvecs)=MoleculePosition(3,atm)-CenterOfMass(3)
        charges(nvecs)=Atom(atype)%Charge
      end if
    end do
    if(nvecs == 0)return

    call CalculateStructureFactor(sys,SimulationCell,StructureFactor)
  end subroutine CalculateMoleculeStructureFactor

  !=============================================================================================================================
  subroutine CalculateStructureFactor(sys,scell,StructureFactor)
    integer, intent(in)                       :: sys
    type(SimcellInfo), intent(in)             :: scell
    type(StructureFactorInfo), intent(inout)  :: StructureFactor

    integer :: kx, ky, kz, nk, a
    complex(PR) :: temp
    real(PR) :: kxrx,kyry,kzrz
    real(PR), dimension(3) :: vec1,vec2

    !** Calculate kx, ky, kz = 0 , -1 and 1 explicitly 
    do a = 1, nvecs
      eikx(a,0) = (1._PR, 0._PR)
      eiky(a,0) = (1._PR, 0._PR)
      eikz(a,0) = (1._PR, 0._PR)

      !** Transform vector to account for non-orthorhombic cells
      if(scell%NonOrthorhombic)then
        vec1=Matrix_x_Vector(scell%RLatVecs(:,:),atvecs(:,a))
        vec2=Matrix_x_Vector(scell%RLatVecs(:,:),comvecs(:,a))
      else
        vec1=atvecs(:,a)
        vec2=comvecs(:,a)
      end if
      kxrx=TWOPI/scell%BoxLength(1)*vec1(1)
      kyry=TWOPI/scell%BoxLength(2)*vec1(2)
      kzrz=TWOPI/scell%BoxLength(3)*vec1(3)

      kxrx2(a,1)=TWOPI/scell%BoxLength(1)*vec2(1)
      kyry2(a,1)=TWOPI/scell%BoxLength(2)*vec2(2)
      kzrz2(a,1)=TWOPI/scell%BoxLength(3)*vec2(3)

      eikx(a,1) = cmplx(cos(kxrx), sin(kxrx), PR)
      eiky(a,1) = cmplx(cos(kyry), sin(kyry), PR)
      eikz(a,1) = cmplx(cos(kzrz), sin(kzrz), PR)

      eiky(a,-1) = conjg(eiky(a,1))
      eikz(a,-1) = conjg(eikz(a,1))

      kyry2(a,-1)=-TWOPI/scell%BoxLength(2)*vec2(2)
      kzrz2(a,-1)=-TWOPI/scell%BoxLength(3)*vec2(3)
    end do

    !** Calculate remaining kx, ky and kz by recurrence 
    do kx = 2, kmax(1,sys)
      do a = 1, nvecs
        eikx(a,kx) = eikx(a,kx-1)*eikx(a,1)
        kxrx2(a,kx)=kxrx2(a,kx-1)+kxrx2(a,1)
      end do
    end do

    do ky = 2, kmax(2,sys)
      do a = 1, nvecs
        eiky(a,ky) = eiky(a,ky-1)*eiky(a,1)
        eiky(a,-ky) = conjg(eiky(a,ky))
        kyry2(a,ky)=kyry2(a,ky-1)+kyry2(a,1)
        kyry2(a,-ky)=-kyry2(a,ky)
      end do
    end do

    do kz = 2, kmax(3,sys)
      do a = 1, nvecs
        eikz(a,kz) = eikz(a,kz-1)*eikz(a,1)
        eikz(a,-kz) = conjg(eikz(a,kz))
        kzrz2(a,kz)=kzrz2(a,kz-1)+kzrz2(a,1)
        kzrz2(a,-kz)=-kzrz2(a,kz)
      end do
    end do

    !** Calculate the structure factor
    nk=0
    do kx=0,kmax(1,sys)!!kmax dimensions
      do ky=-kmax(2,sys),kmax(2,sys)
        do a=1,nvecs
          eikxky(a)=eikx(a,kx)*eiky(a,ky)
        end do
        do kz=-kmax(3,sys),kmax(3,sys)
          if (FourierMask(sys)%mask(kx,ky,kz))cycle
          nk=nk+1
          StructureFactor%RhoK(nk)=0._PR
          StructureFactor%RhokV(nk)=0._PR
          do a=1,nvecs
            eikr(a)=eikxky(a)*eikz(a,kz)
            temp=charges(a)*eikr(a)
            StructureFactor%RhoK(nk)=StructureFactor%RhoK(nk)+temp
            StructureFactor%RhokV(nk)=StructureFactor%RhokV(nk)+temp*(kxrx2(a,kx)+kyry2(a,ky)+kzrz2(a,kz))
          end do
        end do
      end do
    end do

    !** Calculate Dipole
    StructureFactor%Dipole=0._PR
    StructureFactor%DipoleV=0._PR
    do a=1,nvecs
      StructureFactor%Dipole=StructureFactor%Dipole+charges(a)*atvecs(1:3,a)
      StructureFactor%DipoleV=StructureFactor%DipoleV+charges(a)*comvecs(1:3,a)
    end do

  end subroutine CalculateStructureFactor

  !=============================================================================================================================
  subroutine EwaldFourierInteraction(sys,SimulationVolume,Kfactor,StructureFactor,Interactions)
    integer, intent(in)                    :: sys
    real(PR), intent(in)                   :: SimulationVolume
    type(KfactorInfo), intent(in)          :: Kfactor
    type(StructureFactorInfo), intent(in)  :: StructureFactor
    type(StorageInteractions), intent(out) :: Interactions

    integer  :: nk
    real(PR) :: temp

    !** Calculate Ewald Fourier Energy
    do nk=1,NumberOfKVectors(sys)
      temp=real(StructureFactor%RhoK(nk))*real(StructureFactor%RhoK(nk))+aimag(StructureFactor%RhoK(nk))*aimag(StructureFactor%RhoK(nk))
      Interactions%Energy%CoulFourier=Interactions%Energy%CoulFourier+Kfactor%Kfactor(nk)*temp
      Interactions%Virial%CoulFourier=Interactions%Virial%CoulFourier+Kfactor%KfactorV(nk)*temp
      temp=real(StructureFactor%RhokV(nk))*aimag(StructureFactor%RhoK(nk))-aimag(StructureFactor%RhokV(nk))*real(StructureFactor%RhoK(nk))
      Interactions%Virial%CoulFourier=Interactions%Virial%CoulFourier+2._PR*Kfactor%Kfactor(nk)*temp
    end do
    Interactions%Virial%CoulFourier=Interactions%Virial%CoulFourier/3._PR

    !** Surface part
    if(.not. TIN_FOIL_BC)then
      temp=StructureFactor%Dipole(1)**2+StructureFactor%Dipole(2)**2+StructureFactor%Dipole(3)**2
      temp=TWOPI/3._PR/SimulationVolume*temp
      Interactions%Energy%CoulFourier=Interactions%Energy%CoulFourier+temp
      Interactions%Virial%CoulFourier=Interactions%Virial%CoulFourier+temp
      temp=StructureFactor%DipoleV(1)*StructureFactor%Dipole(1)+StructureFactor%DipoleV(2)*StructureFactor%Dipole(2)+ &
        StructureFactor%DipoleV(3)*StructureFactor%Dipole(3)
      temp=2._PR*TWOPI/9._PR/SimulationVolume*temp
      Interactions%Virial%CoulFourier=Interactions%Virial%CoulFourier-temp
    end if

    !** Convert energy units
    if(.not. REDUCED_UNITS)then
      Interactions%Energy%CoulFourier=Interactions%Energy%CoulFourier*COULOMBIC_CONVERSION_FACTOR
      Interactions%Virial%CoulFourier=Interactions%Virial%CoulFourier*COULOMBIC_CONVERSION_FACTOR
    end if

  end subroutine EwaldFourierInteraction

  !=============================================================================================================================
  subroutine TotalIntraMolecularFourierInteraction(sys,Coordinates,energy)
    integer, intent(in)                                :: sys
    type(SpeciesCoordinates), dimension(:), intent(in) :: Coordinates
    real(PR), intent(Out)                              :: energy

    integer :: spc,mol
    real(PR) ::nrg

    energy=0._PR
    do spc=1,NumberOfSpecies
      if(Molecule(spc)%HasPartialCharges)then
        do mol=1,Coordinates(spc)%NumberOfMolecules
          call IntramolecularFourier(spc,Coordinates(spc)%Positions(1:3,mol,:),nrg)
          energy=energy+nrg
        end do
      end if
    end do
  end subroutine TotalIntraMolecularFourierInteraction

  !=============================================================================================================================
  subroutine IntramolecularFourier(spc,Positions,energy)
    integer, intent(In)                  :: spc
    real(PR), dimension(:,:), intent(in) :: Positions
    real(PR), intent(Out)                :: energy

    real(PR)  :: charge1,charge2
    real(PR)  :: sepvec(3),rijsq(1)
    real(PR)  :: rij,erfc_kr
    integer   :: atm1,atm2,type1,type2

    energy=0._PR

    if(COULOMB_SCREENING)then
      do atm1=1,Molecule(spc)%NumberOfAtoms
        type1=Molecule(spc)%AtomType(atm1)
        if(.not. PairCoulombInteraction(type1,type1))cycle
        charge1=Atom(type1)%Charge
        energy=energy+(Ewald_alpha_by_sqrt_pi*exp(-Ewald_kappasq_by_four_alphasq)-Ewald_kappa*ErrorFunctionComplement(Ewald_kappa_by_two_alpha)/2._PR)*charge1**2
        do atm2=atm1+1,Molecule(spc)%NumberOfAtoms
          type2=Molecule(spc)%AtomType(atm2)
          if(.not. PairCoulombInteraction(type1,type2))cycle
          charge2=Atom(type2)%Charge
          sepvec=Positions(:,atm1)-Positions(:,atm2)
          rijsq(1)=sepvec(1)**2+sepvec(2)**2+sepvec(3)**2
          rij=sqrt(rijsq(1))
          erfc_kr=(ErrorFunctionComplement(Ewald_alpha*rij+Ewald_kappa_by_two_alpha)*exp(Ewald_kappa*rij)+ &
            ErrorFunctionComplement(Ewald_alpha*rij-Ewald_kappa_by_two_alpha)*exp(-Ewald_kappa*rij))/2._PR
          energy=energy+charge1*charge2*(1._pr-erfc_kr)/rij
        end do
      end do
    else
      do atm1=1,Molecule(spc)%NumberOfAtoms
        type1=Molecule(spc)%AtomType(atm1)
        if(.not. PairCoulombInteraction(type1,type1))cycle
        charge1=Atom(type1)%Charge
        energy=energy+Ewald_alpha_by_sqrt_pi*charge1**2
        do atm2=atm1+1,Molecule(spc)%NumberOfAtoms
          type2=Molecule(spc)%AtomType(atm2)
          if(.not. PairCoulombInteraction(type1,type2))cycle
          charge2=Atom(type2)%Charge
          sepvec=Positions(:,atm1)-Positions(:,atm2)
          rijsq(1)=sepvec(1)**2+sepvec(2)**2+sepvec(3)**2
          rij=sqrt(rijsq(1))
          erfc_kr=ErrorFunctionComplement(Ewald_alpha*rij)
          energy=energy+charge1*charge2*(1._pr-erfc_kr)/rij
        end do
      end do
    end if

    !** Convert energy units
    if(.not. REDUCED_UNITS)energy=energy*COULOMBIC_CONVERSION_FACTOR

  end subroutine IntramolecularFourier

end module ewaldsum
