module inter
  use consts, only: PR, MAX_NO_OF_SPECIES
  use utils, only: ErrorMessage, Cross_Product
  use atoms_and_molecules, only: GetAtomType, Molecule, NumberOfSpecies
  use config, only: SpeciesCoordinates
  use simcell, only: SimcellInfo, GetMinimumImage
  use storage, only: StorageInteractions, InteractionAccumulator, ResetStorageInteractions, ResetInteractionAccumulator, Operator(+)
  use cell_list, only: CellListInfo, GetCellNumber
  use forcefield, only: PairPotential, TwoBodyNonBondedInteraction, EwaldRealSpaceInteraction, TailCorrections, &
    BondInteraction, AngleInteraction, TorsionInteraction, Molecule_intra, IntraPairInteraction, IntraCoulInteraction
  implicit none
  public
  save

  integer, parameter                    :: MAX_NO_OF_SEPARATION_VECTORS=70000
  integer                               :: NumberOfSeparationVectors
  real(PR), dimension(:,:), allocatable :: SeparationVector, ComSepVec
  real(PR), dimension(:), allocatable   :: SeparationDistanceSq, VirialFactor

contains
  subroutine AllocateMemoryToSeparationVectors(error)
    logical, intent(out) :: error
    integer :: ierror

    error=.false.
    allocate(SeparationVector(3,MAX_NO_OF_SEPARATION_VECTORS), ComSepVec(3,MAX_NO_OF_SEPARATION_VECTORS), &
      SeparationDistanceSq(MAX_NO_OF_SEPARATION_VECTORS), VirialFactor(MAX_NO_OF_SEPARATION_VECTORS), STAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
        'Failed to allocate memory to Arrays. STAT = ',ierror
      error=.true.
      return
    end if

  end subroutine AllocateMemoryToSeparationVectors
  !=============================================================================================
  subroutine TotalShortRangePairwiseInteraction(Coordinates,SimulationCell,CellList,Interactions,overlap,error)
    type(SpeciesCoordinates), dimension(:), intent(in) :: Coordinates
    type(SimcellInfo), intent(in)                      :: SimulationCell
    type(CellListInfo), intent(in)                     :: CellList
    type(StorageInteractions), intent(out)             :: Interactions
    logical, intent(out)                               :: overlap
    logical, intent(out)                               :: error

    integer :: n,neigh,nvecs
    integer :: spc1, spc2
    integer :: atm1, atm2
    integer :: type1, type2
    type(InteractionAccumulator) :: Accumulator

    error=.false.
    overlap=.false.
    !** Loop over atoms of species 1
    do spc1=1,NumberOfSpecies
      do atm1=1,Molecule(spc1)%NumberOfAtoms
        type1=Molecule(spc1)%AtomType(atm1)
        !** Loop over atoms of species 2
        do spc2=spc1,NumberOfSpecies
          do atm2=1,Molecule(spc2)%NumberOfAtoms
            type2=Molecule(spc2)%AtomType(atm2)

            if(CellList%active)then
              call AllSeparationVectorsWithCellList(CellList,Coordinates,SimulationCell,spc1,spc2,atm1,atm2,error)
            else
              call AllSeparationVectors(Coordinates,SimulationCell,spc1,spc2,atm1,atm2,error)
            end if
            if(error)return
            Accumulator=TwoBodyNonBondedInteraction(PairPotential(type1,type2),SeparationDistanceSq(1:NumberOfSeparationVectors),VirialFactor(1:NumberOfSeparationVectors))
            if(Accumulator%overlap)then
              overlap=Accumulator%overlap
              return
            end if
            Interactions%Energy%vdW=Interactions%Energy%vdW+Accumulator%Energy
            Interactions%Virial%vdW=Interactions%Virial%vdW+Accumulator%Virial
            Interactions%dUdKappa%vdW=Interactions%dUdKappa%vdW+Accumulator%dUdKappa

            Accumulator=EwaldRealSpaceInteraction(type1,type2,SeparationDistanceSq(1:NumberOfSeparationVectors),VirialFactor(1:NumberOfSeparationVectors))
            Interactions%Energy%CoulReal=Interactions%Energy%CoulReal+Accumulator%Energy
            Interactions%Virial%CoulReal=Interactions%Virial%CoulReal+Accumulator%Virial
          end do
        end do
      end do
    end do

  end subroutine TotalShortRangePairwiseInteraction

  !================================================================================
  subroutine TotalLongRangePairwiseInteraction(NumberOfMolecules,SimulationVolume,Interactions)
    integer, dimension(MAX_NO_OF_SPECIES), intent(in)  :: NumberOfMolecules
    real(PR), intent(in)                               :: SimulationVolume
    type(StorageInteractions), intent(out)             :: Interactions

    integer :: spc1, spc2
    integer :: atm1, atm2
    integer :: type1, type2
    integer :: nmoles1, nmoles2
    real(PR) :: vol
    type(InteractionAccumulator) :: Accumulator

    do spc1=1,NumberOfSpecies
      do spc2=spc1,NumberOfSpecies
        call ResetInteractionAccumulator(Accumulator)
        do atm1=1,Molecule(spc1)%NumberOfAtoms
          type1=Molecule(spc1)%AtomType(atm1)
          do atm2=1,Molecule(spc2)%NumberOfAtoms
            type2=Molecule(spc2)%AtomType(atm2)
            Accumulator=Accumulator+TailCorrections(PairPotential(type1,type2))
          end do
        end do

        nmoles1=NumberOfMolecules(spc1)
        nmoles2=NumberOfMolecules(spc2)
        vol=SimulationVolume
        if(spc1 == spc2)then
          Interactions%Energy%vdW=Interactions%Energy%vdW+real(nmoles1*nmoles2,PR)/vol*Accumulator%Energy
          Interactions%Virial%vdW=Interactions%Virial%vdW+real(nmoles1*nmoles2,PR)/vol*Accumulator%Virial
          Interactions%dUdKappa%vdW=Interactions%dUdKappa%vdW+real(nmoles1*nmoles2,PR)/vol*Accumulator%dUdKappa
        else
          Interactions%Energy%vdW=Interactions%Energy%vdW+2._PR*real(nmoles1*nmoles2,PR)/vol*Accumulator%Energy
          Interactions%Virial%vdW=Interactions%Virial%vdW+2._PR*real(nmoles1*nmoles2,PR)/vol*Accumulator%Virial
          Interactions%dUdKappa%vdW=Interactions%dUdKappa%vdW+2._PR*real(nmoles1*nmoles2,PR)/vol*Accumulator%dUdKappa
        end if
      end do
    end do
  end subroutine TotalLongRangePairwiseInteraction

  !================================================================================
  subroutine MoleculeSystemShortRangePairwiseInteraction(MoleculePosition,CenterOfMass,Coordinates,SimulationCell,CellList,spc1,mol1,Interactions,overlap)
    real(PR), dimension(:,:), intent(in)               :: MoleculePosition
    real(PR), dimension(3), intent(in)                 :: CenterOfMass
    type(SpeciesCoordinates), dimension(:), intent(in) :: Coordinates
    type(SimcellInfo), intent(in)                      :: SimulationCell
    type(CellListInfo), intent(in)                     :: CellList
    integer, intent(in)                                :: spc1,mol1
    type(StorageInteractions), intent(out)             :: Interactions
    logical, intent(out)                               :: overlap

    integer :: atm1
    real(PR), dimension(3)  :: AtomPosition
    type(StorageInteractions)    :: Interactions1

    overlap=.false.
    call ResetStorageInteractions(Interactions)
    do atm1=1,Molecule(spc1)%NumberOfAtoms
      AtomPosition=MoleculePosition(:,atm1)
      call AtomSystemShortRangePairwiseInteraction(AtomPosition,CenterOfMass,Coordinates,SimulationCell,CellList,spc1,mol1,atm1,Interactions1,overlap)
      if(overlap)return
      Interactions=Interactions+Interactions1
    end do

  end subroutine MoleculeSystemShortRangePairwiseInteraction

  !================================================================================
  subroutine AtomSystemShortRangePairwiseInteraction(AtomPosition,CenterOfMass,Coordinates,SimulationCell,CellList,spc1,mol1,atm1,Interactions,overlap)
    real(PR), dimension(3), intent(in)                 :: AtomPosition,CenterOfMass
    type(SpeciesCoordinates), dimension(:), intent(in) :: Coordinates
    type(SimcellInfo), intent(in)                      :: SimulationCell
    type(CellListInfo), intent(in)                     :: CellList
    integer, intent(in)                                :: spc1,mol1,atm1
    type(StorageInteractions), intent(out)             :: Interactions
    logical, intent(out)                               :: overlap

    integer :: spc2,mol2,atm2
    integer :: type1,type2
    type(InteractionAccumulator) :: Accumulator

    overlap=.false.
    type1=Molecule(spc1)%AtomType(atm1)
    do spc2=1,NumberOfSpecies
      do atm2=1,Molecule(spc2)%NumberOfAtoms
        type2=Molecule(spc2)%AtomType(atm2)

        if(CellList%active)then
          call AtomSystemSeparationVectorsWithCellList(CellList,AtomPosition,CenterOfMass,Coordinates,SimulationCell,spc1,spc2,atm2,mol1)
        else
          call AtomSystemSeparationVectors(AtomPosition,CenterOfMass,Coordinates,SimulationCell,spc1,spc2,atm2,mol1)
        end if

        Accumulator=TwoBodyNonBondedInteraction(PairPotential(type1,type2),SeparationDistanceSq(1:NumberOfSeparationVectors),VirialFactor(1:NumberOfSeparationVectors))
        if(Accumulator%overlap)then
          overlap=Accumulator%overlap
          return
        end if
        Interactions%Energy%vdW=Interactions%Energy%vdW+Accumulator%Energy
        Interactions%Virial%vdW=Interactions%Virial%vdW+Accumulator%Virial
        Interactions%dUdKappa%vdW=Interactions%dUdKappa%vdW+Accumulator%dUdKappa

        Accumulator=EwaldRealSpaceInteraction(type1,type2,SeparationDistanceSq(1:NumberOfSeparationVectors),VirialFactor(1:NumberOfSeparationVectors))
        Interactions%Energy%CoulReal=Interactions%Energy%CoulReal+Accumulator%Energy
        Interactions%Virial%CoulReal=Interactions%Virial%CoulReal+Accumulator%Virial
      end do
    end do

  end subroutine AtomSystemShortRangePairwiseInteraction

  !================================================================================
  subroutine MoleculeSystemLongRangePairwiseInteraction(spc1,nmoles1,NumberOfMolecules,SimulationVolume,Interactions)
    integer, intent(in)                                :: spc1,nmoles1
    integer, dimension(MAX_NO_OF_SPECIES), intent(in) :: NumberOfMolecules
    real(PR), intent(in)                               :: SimulationVolume
    type(StorageInteractions), intent(out)             :: Interactions

    integer :: spc2
    integer :: atm1, atm2
    integer :: type1, type2
    integer :: nmoles2
    real(PR) :: vol
    type(InteractionAccumulator) :: Accumulator, Accumulator1

    do spc2=1,NumberOfSpecies
      call ResetInteractionAccumulator(Accumulator)
      do atm1=1,Molecule(spc1)%NumberOfAtoms
        type1=Molecule(spc1)%AtomType(atm1)
        do atm2=1,Molecule(spc2)%NumberOfAtoms
          type2=Molecule(spc2)%AtomType(atm2)
          Accumulator=Accumulator+TailCorrections(PairPotential(type1,type2))
        end do
      end do

      nmoles2=NumberOfMolecules(spc2)
      vol=SimulationVolume
      if(spc1 == spc2)then
        Interactions%Energy%vdW=Interactions%Energy%vdW+real(2*nmoles1-1,PR)/vol*Accumulator%Energy
        Interactions%Virial%vdW=Interactions%Virial%vdW+real(2*nmoles1-1,PR)/vol*Accumulator%Virial
        Interactions%dUdKappa%vdW=Interactions%dUdKappa%vdW+real(2*nmoles1-1,PR)/vol*Accumulator%dUdKappa
      else
        Interactions%Energy%vdW=Interactions%Energy%vdW+2._PR*real(nmoles2,PR)/vol*Accumulator%Energy
        Interactions%Virial%vdW=Interactions%Virial%vdW+2._PR*real(nmoles2,PR)/vol*Accumulator%Virial
        Interactions%dUdKappa%vdW=Interactions%dUdKappa%vdW+2._PR*real(nmoles2,PR)/vol*Accumulator%dUdKappa
      end if
    end do
  end subroutine MoleculeSystemLongRangePairwiseInteraction

  !================================================================================
  !** Intra-molecular interaction
  !================================================================================
  subroutine TotalIntraMolecularInteractions(Coordinates,Interactions,overlap)
    type(SpeciesCoordinates), dimension(:), intent(in) :: Coordinates
    type(StorageInteractions), intent(out)             :: Interactions
    logical, intent(out)                               :: overlap

    integer :: spc, mol

    do spc=1,NumberOfSpecies
      do mol=1,Coordinates(spc)%NumberOfMolecules
        Interactions=Interactions+IntraMolecularInteraction(Coordinates(spc)%Positions(:,mol,:),spc,mol,overlap)
      end do
    end do

  end subroutine TotalIntraMolecularInteractions

  !================================================================================
  function IntraMolecularInteraction(MoleculePosition,spc,mol,overlap) result(Interactions)
    real(PR), dimension(:,:), intent(in)               :: MoleculePosition
    integer, intent(in)                                :: spc,mol
    logical, intent(out)                               :: overlap
    type(StorageInteractions)                          :: Interactions

    integer :: i, atm1, atm2, atm3, atm4
    real(PR) :: rij, rijsq, theta, phi
    real(PR), dimension(3) :: vec1, vec2, vecX, vecY, vecZ
    real(PR) :: energy

    overlap=.false.
    !** Bond interaction
    do i=1,Molecule_intra(spc)%NumberOfBonds
      atm1=Molecule_intra(spc)%BondAtom(i,1)
      atm2=Molecule_intra(spc)%BondAtom(i,2)
      rij=norm2(MoleculePosition(:,atm1)-MoleculePosition(:,atm2))
      energy=BondInteraction(Molecule_intra(spc)%Bond(i),rij)
      Interactions%Energy%Bond=Interactions%Energy%Bond+energy
    end do

    !** Angle interaction
    do i=1,Molecule_intra(spc)%NumberOfAngles
      atm1=Molecule_intra(spc)%AngleAtom(i,1)
      atm2=Molecule_intra(spc)%AngleAtom(i,2)
      atm3=Molecule_intra(spc)%AngleAtom(i,3)
      vec1=MoleculePosition(:,atm1)-MoleculePosition(:,atm2)
      vec2=MoleculePosition(:,atm3)-MoleculePosition(:,atm2)
      theta=dot_product(vec1,vec2)/(norm2(vec1)*norm2(vec2))
      energy=AngleInteraction(Molecule_intra(spc)%Angle(i),theta)
      Interactions%Energy%Angle=Interactions%Energy%Angle+energy
    end do

    !** Torsion interaction
    do i=1,Molecule_intra(spc)%NumberOfTorsions
      atm1=Molecule_intra(spc)%TorsionAtom(i,1)
      atm2=Molecule_intra(spc)%TorsionAtom(i,2)
      atm3=Molecule_intra(spc)%TorsionAtom(i,3)
      atm4=Molecule_intra(spc)%TorsionAtom(i,4)
      vec1=MoleculePosition(:,1)-MoleculePosition(:,2)
      vecZ=MoleculePosition(:,2)-MoleculePosition(:,3)
      vec2=MoleculePosition(:,3)-MoleculePosition(:,4)
      vecZ=vecZ/norm2(vecZ)
      vecY=Cross_Product(vec2,vecZ)
      vecY=vecZ/norm2(vecY)
      vecX=Cross_Product(vecY,vecZ)
      phi=atan2(dot_product(vec1,vecY),dot_product(vec1,vecX))
      energy=TorsionInteraction(Molecule_intra(spc)%Torsion(i),phi)
      Interactions%Energy%Torsion=Interactions%Energy%Torsion+energy
    end do

    !** Intra Pair interaction
    do i=1,Molecule_intra(spc)%NumberOfIntraPairs
      atm1=Molecule_intra(spc)%IntraPairAtom(i,1)
      atm2=Molecule_intra(spc)%IntraPairAtom(i,2)
      rijsq=sum((MoleculePosition(:,atm1)-MoleculePosition(:,atm2))**2)
      energy=IntraPairInteraction(Molecule_intra(spc)%IntraPair(i),rijsq,overlap)
      if(overlap)return
      Interactions%Energy%IntraPair=Interactions%Energy%IntraPair+energy
    end do

    !** Intra Coul interaction
    do i=1,Molecule_intra(spc)%NumberOfIntraCouls
      atm1=Molecule_intra(spc)%IntraCoulAtom(i,1)
      atm2=Molecule_intra(spc)%IntraCoulAtom(i,2)
      rijsq=sum((MoleculePosition(:,atm1)-MoleculePosition(:,atm2))**2)
      energy=IntraCoulInteraction(Molecule_intra(spc)%IntraCoul(i),rijsq)
      Interactions%Energy%IntraCoul=Interactions%Energy%IntraCoul+energy
    end do

  end function IntraMolecularInteraction

  !================================================================================
  !      Separation Vectors
  !================================================================================
  subroutine AllSeparationVectors(Coordinates,SimulationCell,spc1,spc2,atm1,atm2,error)
    type(SpeciesCoordinates), dimension(:), intent(in) :: Coordinates
    type(SimcellInfo), intent(in)                      :: SimulationCell
    integer, intent(in)    :: spc1, spc2, atm1, atm2
    logical, intent(out)   :: error

    integer :: mol1, mol2
    integer :: n,low
    integer :: required_size, ierror

    error=.false.
    required_size=Coordinates(spc1)%NumberOfMolecules*Coordinates(spc2)%NumberOfMolecules
    if(required_size > size(SeparationVector,2))then
      deallocate(SeparationVector,SeparationDistanceSq,ComSepVec,VirialFactor,STAT=ierror)
      if(ierror /=0)then
        write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
          'Deallocation Failure. ierror = ',ierror
        error=.true.
        return
      end if
      allocate(SeparationVector(3,required_size), &
        SeparationDistanceSq(required_size), &
        ComSepVec(3,required_size), &
        VirialFactor(required_size), STAT=ierror)
      if(ierror /=0)then
        write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
          'Allocation Failure. ierror = ',ierror
        error=.true.
        return
      end if
    end if

    n=0
    if(Molecule(spc1)%NumberOfAtoms == 1 .and. Molecule(spc2)%NumberOfAtoms == 1)then
      do mol1=1,Coordinates(spc1)%NumberOfMolecules
        if(spc1 == spc2)then
          low=mol1+1
        else
          low=1
        end if
        do mol2=low,Coordinates(spc2)%NumberOfMolecules
          n=n+1
          SeparationVector(1,n)=Coordinates(spc1)%Positions(1,mol1,atm1)-Coordinates(spc2)%Positions(1,mol2,atm2)
          SeparationVector(2,n)=Coordinates(spc1)%Positions(2,mol1,atm1)-Coordinates(spc2)%Positions(2,mol2,atm2)
          SeparationVector(3,n)=Coordinates(spc1)%Positions(3,mol1,atm1)-Coordinates(spc2)%Positions(3,mol2,atm2)
          VirialFactor(n)=1._PR
        end do
      end do
      NumberOfSeparationVectors=n
      call GetMinimumImage(SimulationCell,NumberOfSeparationVectors,SeparationVector,SeparationDistanceSq)
    else
      do mol1=1,Coordinates(spc1)%NumberOfMolecules
        if(spc1 == spc2)then
          low=mol1+1
        else
          low=1
        end if
        do mol2=low,Coordinates(spc2)%NumberOfMolecules
          n=n+1
          SeparationVector(1,n)=Coordinates(spc1)%Positions(1,mol1,atm1)-Coordinates(spc2)%Positions(1,mol2,atm2)
          SeparationVector(2,n)=Coordinates(spc1)%Positions(2,mol1,atm1)-Coordinates(spc2)%Positions(2,mol2,atm2)
          SeparationVector(3,n)=Coordinates(spc1)%Positions(3,mol1,atm1)-Coordinates(spc2)%Positions(3,mol2,atm2)
          ComSepVec(1,n)=Coordinates(spc1)%CenterOfMass(1,mol1)-Coordinates(spc2)%CenterOfMass(1,mol2)
          ComSepVec(2,n)=Coordinates(spc1)%CenterOfMass(2,mol1)-Coordinates(spc2)%CenterOfMass(2,mol2)
          ComSepVec(3,n)=Coordinates(spc1)%CenterOfMass(3,mol1)-Coordinates(spc2)%CenterOfMass(3,mol2)
        end do
      end do
      NumberOfSeparationVectors=n
      call GetMinimumImage(SimulationCell,NumberOfSeparationVectors,SeparationVector,SeparationDistanceSq,ComSepVec,VirialFactor)
    end if

  end subroutine AllSeparationVectors

  !================================================================================
  subroutine AtomSystemSeparationVectors(AtomPosition,CenterOfMass,Coordinates,SimulationCell,spc1,spc2,atm2,mol1)
    real(PR), dimension(3), intent(in) :: AtomPosition,CenterOfMass
    type(SpeciesCoordinates), dimension(:), intent(in) :: Coordinates
    type(SimcellInfo), intent(in)                      :: SimulationCell
    integer, intent(in) :: spc1,spc2,atm2,mol1

    integer :: mol2
    integer :: n,loop,nloops,low(2),high(2)

    n=0
    if(spc1 == spc2)then
      nloops=2
      low(1)=1
      high(1)=mol1-1
      low(2)=mol1+1
      high(2)=Coordinates(spc2)%NumberOfMolecules
    else
      nloops=1
      low(1)=1
      high(1)=Coordinates(spc2)%NumberOfMolecules
    end if

    if(Molecule(spc1)%NumberOfAtoms == 1 .and. Molecule(spc2)%NumberOfAtoms == 1)then
      do loop=1,nloops
        do mol2=low(loop),high(loop)
          n=n+1
          SeparationVector(1,n)=AtomPosition(1)-Coordinates(spc2)%Positions(1,mol2,atm2)
          SeparationVector(2,n)=AtomPosition(2)-Coordinates(spc2)%Positions(2,mol2,atm2)
          SeparationVector(3,n)=AtomPosition(3)-Coordinates(spc2)%Positions(3,mol2,atm2)
          VirialFactor(n)=1._PR
        end do
      end do
      NumberOfSeparationVectors=n
      call GetMinimumImage(SimulationCell,NumberOfSeparationVectors,SeparationVector,SeparationDistanceSq)
    else
      do loop=1,nloops
        do mol2=low(loop),high(loop)
          n=n+1
          SeparationVector(1,n)=AtomPosition(1)-Coordinates(spc2)%Positions(1,mol2,atm2)
          SeparationVector(2,n)=AtomPosition(2)-Coordinates(spc2)%Positions(2,mol2,atm2)
          SeparationVector(3,n)=AtomPosition(3)-Coordinates(spc2)%Positions(3,mol2,atm2)
          ComSepVec(1,n)=CenterOfMass(1)-Coordinates(spc2)%CenterOfMass(1,mol2)
          ComSepVec(2,n)=CenterOfMass(2)-Coordinates(spc2)%CenterOfMass(2,mol2)
          ComSepVec(3,n)=CenterOfMass(3)-Coordinates(spc2)%CenterOfMass(3,mol2)
        end do
      end do
      NumberOfSeparationVectors=n
      call GetMinimumImage(SimulationCell,NumberOfSeparationVectors,SeparationVector,SeparationDistanceSq,ComSepVec,VirialFactor)
    end if

  end subroutine AtomSystemSeparationVectors

  !=============================================================================================
  !          Calculation of Separation Vectors between two species
  !=============================================================================================
  subroutine AllSeparationVectorsWithCellList(CellList,Coordinates,SimulationCell,spc1,spc2,atm1,atm2,error)
    type(CellListInfo), intent(in)                     :: CellList
    type(SpeciesCoordinates), dimension(:), intent(in) :: Coordinates
    type(SimcellInfo), intent(in)                      :: SimulationCell
    integer, intent(in)    :: spc1, spc2, atm1, atm2
    logical, intent(out)   :: error

    integer :: mol1, mol2, n
    integer :: icel,jcel,ncel
    integer :: required_size, ierror

    error=.false.
    required_size=Coordinates(spc1)%NumberOfMolecules*Coordinates(spc2)%NumberOfMolecules
    if(required_size > size(SeparationVector,2))then
      deallocate(SeparationVector,SeparationDistanceSq,ComSepVec,VirialFactor,STAT=ierror)
      if(ierror /=0)then
        write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
          'Deallocation Failure. ierror = ',ierror
        error=.true.
        return
      end if
      allocate(SeparationVector(3,required_size), &
        SeparationDistanceSq(required_size), &
        ComSepVec(3,required_size), &
        VirialFactor(required_size), STAT=ierror)
      if(ierror /=0)then
        write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
          'Allocation Failure. ierror = ',ierror
        error=.true.
        return
      end if
    end if

    n=0
    if(Molecule(spc1)%NumberOfAtoms == 1 .and. Molecule(spc2)%NumberOfAtoms == 1)then
      do mol1=1,Coordinates(spc1)%NumberOfMolecules
        icel=GetCellNumber(Coordinates(spc1)%Positions(:,mol1,atm1),SimulationCell,CellList)
        do ncel=1,27
          jcel=CellList%CellNeighbors(ncel,icel)
          mol2=CellList%SpeciesCellList(spc2)%HeadOfCell(jcel,atm2)
          do while(mol2 /= 0)
            if(spc2 /= spc1 .or. mol2 > mol1)then
              n=n+1
              SeparationVector(1,n)=Coordinates(spc1)%Positions(1,mol1,atm1)-Coordinates(spc2)%Positions(1,mol2,atm2)
              SeparationVector(2,n)=Coordinates(spc1)%Positions(2,mol1,atm1)-Coordinates(spc2)%Positions(2,mol2,atm2)
              SeparationVector(3,n)=Coordinates(spc1)%Positions(3,mol1,atm1)-Coordinates(spc2)%Positions(3,mol2,atm2)
              VirialFactor(n)=1._PR
            end if
            mol2=CellList%SpeciesCellList(spc2)%LinkedList(mol2,atm2)
          end do
        end do
      end do
      NumberOfSeparationVectors=n
      call GetMinimumImage(SimulationCell,NumberOfSeparationVectors,SeparationVector,SeparationDistanceSq)
    else
      do mol1=1,Coordinates(spc1)%NumberOfMolecules
        icel=GetCellNumber(Coordinates(spc1)%Positions(:,mol1,atm1),SimulationCell,CellList)
        do ncel=1,27
          jcel=CellList%CellNeighbors(ncel,icel)
          mol2=CellList%SpeciesCellList(spc2)%HeadOfCell(jcel,atm2)
          do while(mol2 /= 0)
            if(spc2 /= spc1 .or. mol2 > mol1)then
              n=n+1
              SeparationVector(1,n)=Coordinates(spc1)%Positions(1,mol1,atm1)-Coordinates(spc2)%Positions(1,mol2,atm2)
              SeparationVector(2,n)=Coordinates(spc1)%Positions(2,mol1,atm1)-Coordinates(spc2)%Positions(2,mol2,atm2)
              SeparationVector(3,n)=Coordinates(spc1)%Positions(3,mol1,atm1)-Coordinates(spc2)%Positions(3,mol2,atm2)
              ComSepVec(1,n)=Coordinates(spc1)%CenterOfMass(1,mol1)-Coordinates(spc2)%CenterOfMass(1,mol2)
              ComSepVec(2,n)=Coordinates(spc1)%CenterOfMass(2,mol1)-Coordinates(spc2)%CenterOfMass(2,mol2)
              ComSepVec(3,n)=Coordinates(spc1)%CenterOfMass(3,mol1)-Coordinates(spc2)%CenterOfMass(3,mol2)
            end if
            mol2=CellList%SpeciesCellList(spc2)%LinkedList(mol2,atm2)
          end do
        end do
      end do
      NumberOfSeparationVectors=n
      call GetMinimumImage(SimulationCell,NumberOfSeparationVectors,SeparationVector,SeparationDistanceSq,ComSepVec,VirialFactor)
    end if

  end subroutine AllSeparationVectorsWithCellList

  !=============================================================================================
  !          Calculation of Separation Vectors for one molecule
  !================================================================================
  subroutine AtomSystemSeparationVectorsWithCellList(CellList,AtomPosition,CenterOfMass,Coordinates,SimulationCell,spc1,spc2,atm2,mol1)
    type(CellListInfo), intent(in)                     :: CellList
    real(PR), dimension(3), intent(in)                 :: AtomPosition,CenterOfMass
    type(SpeciesCoordinates), dimension(:), intent(in) :: Coordinates
    type(SimcellInfo), intent(in)                      :: SimulationCell
    integer, intent(in) :: spc1,spc2,atm2,mol1

    integer :: mol2, n
    integer :: icel,jcel,ncel

    n=0
    icel=GetCellNumber(AtomPosition,SimulationCell,CellList)
    if(Molecule(spc1)%NumberOfAtoms == 1 .and. Molecule(spc2)%NumberOfAtoms == 1)then
      do ncel=1,27
        jcel=CellList%CellNeighbors(ncel,icel)
        mol2=CellList%SpeciesCellList(spc2)%HeadOfCell(jcel,atm2)
        do while(mol2 /= 0)
          if(spc1 /= spc2 .or. mol1 /= mol2)then
            n=n+1
            SeparationVector(1,n)=AtomPosition(1)-Coordinates(spc2)%Positions(1,mol2,atm2)
            SeparationVector(2,n)=AtomPosition(2)-Coordinates(spc2)%Positions(2,mol2,atm2)
            SeparationVector(3,n)=AtomPosition(3)-Coordinates(spc2)%Positions(3,mol2,atm2)
            VirialFactor(n)=1._PR
          end if
          mol2=CellList%SpeciesCellList(spc2)%LinkedList(mol2,atm2)
        end do
      end do
      NumberOfSeparationVectors=n
      call GetMinimumImage(SimulationCell,NumberOfSeparationVectors,SeparationVector,SeparationDistanceSq)
    else
      do ncel=1,27
        jcel=CellList%CellNeighbors(ncel,icel)
        mol2=CellList%SpeciesCellList(spc2)%HeadOfCell(jcel,atm2)
        do while(mol2 /= 0)
          if(spc1 /= spc2 .or. mol1 /= mol2)then
            n=n+1
            SeparationVector(1,n)=AtomPosition(1)-Coordinates(spc2)%Positions(1,mol2,atm2)
            SeparationVector(2,n)=AtomPosition(2)-Coordinates(spc2)%Positions(2,mol2,atm2)
            SeparationVector(3,n)=AtomPosition(3)-Coordinates(spc2)%Positions(3,mol2,atm2)
            ComSepVec(1,n)=CenterOfMass(1)-Coordinates(spc2)%CenterOfMass(1,mol2)
            ComSepVec(2,n)=CenterOfMass(2)-Coordinates(spc2)%CenterOfMass(2,mol2)
            ComSepVec(3,n)=CenterOfMass(3)-Coordinates(spc2)%CenterOfMass(3,mol2)
          end if
          mol2=CellList%SpeciesCellList(spc2)%LinkedList(mol2,atm2)
        end do
      end do
      NumberOfSeparationVectors=n
      call GetMinimumImage(SimulationCell,NumberOfSeparationVectors,SeparationVector,SeparationDistanceSq,ComSepVec,VirialFactor)
    end if

  end subroutine AtomSystemSeparationVectorsWithCellList

end module inter

