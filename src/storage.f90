module storage
  use consts, only: PR, NO_OF_SUBSETS, UNIT_ENERGY, AVOGADRO_NUMBER
  implicit none

  !=============================================================================
  !** Storage Data Structure
  type :: TopStorage
    real(PR)  :: Total=0._PR
    real(PR)  :: vdW=0._PR
    real(PR)  :: CoulReal=0._PR
    real(PR)  :: CoulFourier=0._PR
    real(PR)  :: IntraFourier=0._PR
    real(PR)  :: Einstein=0._PR
    real(PR)  :: Bond=0._PR
    real(PR)  :: Angle=0._PR
    real(PR)  :: Torsion=0._PR
    real(PR)  :: IntraPair=0._PR
    real(PR)  :: IntraCoul=0._PR
  end type TopStorage

  !=============================================================================
  type :: StorageInteractions
    type(TopStorage) :: Energy           !** Potential Energy of System
    type(TopStorage) :: Virial           !** Virial: used to calculate Pressure
    type(TopStorage) :: dUdKappa         !** Derivative of Energy w.r.t. kappa in Yukawa forcefield
    logical          :: overlap=.false.  !** Overlap of atoms
  end type StorageInteractions

  !=============================================================================
  type :: InteractionAccumulator
    real(PR) :: Energy=0._PR, Virial=0._PR, dUdKappa=0._pr
    logical :: overlap=.false.
  end type InteractionAccumulator

  !=============================================================================
  interface operator (+)
    module procedure storageinteractions_p_storageinteractions
    module procedure accumulator_p_accumulator
  end interface

  interface operator (-)
    module procedure storageinteractions_m_storageinteractions
  end interface

  public :: StorageInteractions
  public :: ResetStorageInteractions, ResetInteractionAccumulator, UpdateStorageInteractions
  public :: operator(+), operator(-)
  public :: DisplayStorage

contains
  !=============================================================================
  function accumulator_p_accumulator(acc1,acc2) result(acc)
    type(InteractionAccumulator), intent(in) :: acc1,acc2
    type(InteractionAccumulator)             :: acc

    acc%Energy=acc1%Energy+acc2%Energy
    acc%Virial=acc1%Virial+acc2%Virial
    acc%dUdKappa=acc1%dUdKappa+acc2%dUdKappa
    acc%overlap=(acc1%overlap .or. acc2%overlap)
  end function accumulator_p_accumulator
  !=============================================================================
  elemental function storage_p_storage(inter1,inter2) result(inter)
    type(TopStorage), intent(In) :: inter1,inter2
    type(TopStorage)             :: inter

    inter%vdW=inter1%vdW+inter2%vdW
    inter%CoulReal=inter1%CoulReal+inter2%CoulReal
    inter%CoulFourier=inter1%CoulFourier+inter2%CoulFourier
    inter%IntraFourier=inter1%IntraFourier+inter2%IntraFourier
    inter%Einstein=inter1%Einstein+inter2%Einstein
    inter%Bond=inter1%Bond+inter2%Bond
    inter%Angle=inter1%Angle+inter2%Angle
    inter%Torsion=inter1%Torsion+inter2%Torsion
    inter%IntraPair=inter1%IntraPair+inter2%IntraPair
    inter%IntraCoul=inter1%IntraCoul+inter2%IntraCoul

  end function storage_p_storage

  !=============================================================================
  elemental function storage_m_storage(inter1,inter2) result(inter)
    type(TopStorage), intent(In) :: inter1,inter2
    type(TopStorage)             :: inter

    inter%vdW=inter1%vdW-inter2%vdW
    inter%CoulReal=inter1%CoulReal-inter2%CoulReal
    inter%CoulFourier=inter1%CoulFourier-inter2%CoulFourier
    inter%IntraFourier=inter1%IntraFourier-inter2%IntraFourier
    inter%Einstein=inter1%Einstein-inter2%Einstein
    inter%Bond=inter1%Bond-inter2%Bond
    inter%Angle=inter1%Angle-inter2%Angle
    inter%Torsion=inter1%Torsion-inter2%Torsion
    inter%IntraPair=inter1%IntraPair-inter2%IntraPair
    inter%IntraCoul=inter1%IntraCoul-inter2%IntraCoul

  end function storage_m_storage

  !=============================================================================
  subroutine UpdateStorage(inter)
    type(TopStorage), intent(InOut) :: inter

    inter%Total=inter%vdW+inter%CoulReal+inter%CoulFourier-inter%IntraFourier + &
      inter%Bond + inter%Angle + inter%Torsion + inter%IntraPair + inter%IntraCoul
  end subroutine UpdateStorage

  !=============================================================================
  elemental subroutine ResetStorageInteractions(interactions)
    type(StorageInteractions), Intent(Out) :: interactions
    return
  end subroutine ResetStorageInteractions

  !=============================================================================
  function storageinteractions_p_storageinteractions(interactions1,interactions2) result(interactions)
    type(StorageInteractions), intent(In) :: interactions1,interactions2
    type(StorageInteractions)             :: interactions

    interactions%Energy=storage_p_storage(interactions1%Energy,interactions2%Energy)
    interactions%Virial=storage_p_storage(interactions1%Virial,interactions2%Virial)
    interactions%dUdKappa=storage_p_storage(interactions1%dUdKappa,interactions2%dUdKappa)
    interactions%overlap=(interactions1%overlap .or. interactions2%overlap)
  end function storageinteractions_p_storageinteractions

  !=============================================================================
  function storageinteractions_m_storageinteractions(interactions1,interactions2) result(interactions)
    type(StorageInteractions), intent(In) :: interactions1,interactions2
    type(StorageInteractions)             :: interactions

    interactions%Energy=storage_m_storage(interactions1%Energy,interactions2%Energy)
    interactions%Virial=storage_m_storage(interactions1%Virial,interactions2%Virial)
    interactions%dUdKappa=storage_m_storage(interactions1%dUdKappa,interactions2%dUdKappa)
    interactions%overlap=(interactions1%overlap .or. interactions2%overlap)
  end function storageinteractions_m_storageinteractions

  !=============================================================================
  subroutine UpdateStorageInteractions(interactions)
    type(StorageInteractions), Intent(InOut) :: interactions

    call UpdateStorage(interactions%Energy)
    call UpdateStorage(interactions%Virial)
    call UpdateStorage(interactions%dUdKappa)
  end subroutine UpdateStorageInteractions

  !=============================================================================
  elemental subroutine ResetInteractionAccumulator(Accumulator)
    type(InteractionAccumulator), Intent(Out) :: Accumulator
    return
  end subroutine ResetInteractionAccumulator

  !=============================================================================================================================
  subroutine DisplayStorage(str,strname,unitno,REDUCED_UNITS)
    type(TopStorage), intent(In) :: str
    character(*), intent(In) :: strname
    integer, intent(IN) :: unitno
    logical, intent(in) :: REDUCED_UNITS

    integer :: spc1,spc2

    if(REDUCED_UNITS)then
      write(unitno,*)
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,'(a)')trim(strname)
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,101)'System','Total ',':',str%Total
      write(unitno,100)'van der Waals ',':',str%vdW
      write(unitno,100)'Real Coulombic ',':',str%CoulReal
      write(unitno,100)'Fourier Coulombic ',':',str%CoulFourier-str%IntraFourier
      write(unitno,100)'Einstein ',':',str%Einstein
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,*)
    else
      write(unitno,*)
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,'(a)')trim(strname)
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,101)'System','Total(J/mol) ',':',str%Total*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,100)'van der Waals(J/mol) ',':',str%vdW*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,100)'Real Coulombic(J/mol) ',':',str%CoulReal*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,100)'Fourier Coulombic(J/mol) ',':',(str%CoulFourier-str%IntraFourier)*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,100)'Einstein(J/mol) ',':',str%Einstein*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,*)
    end if

100 format(t20,a,t50,a,e20.6)
101 format(a,t20,a,t50,a,e20.6)
  end subroutine DisplayStorage

  !=============================================================================
end module storage
