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
  public :: DisplayStorage,CompareAndDisplayStorage

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
      write(unitno,'(t30,a)')trim(strname)
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,100)'Total ',':',str%Total
      write(unitno,100)'van der Waals ',':',str%vdW
      write(unitno,100)'Real Coulombic ',':',str%CoulReal
      write(unitno,100)'Fourier Coulombic ',':',str%CoulFourier-str%IntraFourier
      write(unitno,100)'Einstein ',':',str%Einstein
      write(unitno,100)'Bond (J/mol) ',':',str%Bond
      write(unitno,100)'Angle (J/mol) ',':',str%Angle
      write(unitno,100)'Torsion (J/mol) ',':',str%Torsion
      write(unitno,100)'IntraPair (J/mol) ',':',str%IntraPair
      write(unitno,100)'IntraCoul (J/mol) ',':',str%IntraCoul
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,*)
    else
      write(unitno,*)
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,'(t30,a)')trim(strname)
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,100)'Total(J/mol) ',':',str%Total*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,100)'van der Waals(J/mol) ',':',str%vdW*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,100)'Real Coulombic(J/mol) ',':',str%CoulReal*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,100)'Fourier Coulombic(J/mol) ',':',(str%CoulFourier-str%IntraFourier)*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,100)'Einstein(J/mol) ',':',str%Einstein*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,100)'Bond (J/mol) ',':',str%Bond*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,100)'Angle (J/mol) ',':',str%Angle*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,100)'Torsion (J/mol) ',':',str%Torsion*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,100)'IntraPair (J/mol) ',':',str%IntraPair*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,100)'IntraCoul (J/mol) ',':',str%IntraCoul*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,*)
    end if
 
100 format(a,t30,a,e20.6)
  end subroutine DisplayStorage

  !=============================================================================================================================
  subroutine CompareAndDisplayStorage(str1,strname1,str2,strname2,unitno,REDUCED_UNITS)
    type(TopStorage), intent(In) :: str1,str2
    character(*), intent(In) :: strname1,strname2
    integer, intent(IN) :: unitno
    logical, intent(in) :: REDUCED_UNITS

    integer :: spc1,spc2

    if(REDUCED_UNITS)then
      write(unitno,*)
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,'(t30,a,t60,a)')trim(strname1),trim(strname2)
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,200)'Total ',':',str1%Total,str2%Total
      write(unitno,200)'van der Waals ',':',str1%vdW,str2%vdW
      write(unitno,200)'Real Coulombic ',':',str1%CoulReal,str2%CoulReal
      write(unitno,200)'Fourier Coulombic ',':',str1%CoulFourier-str1%IntraFourier,str2%CoulFourier-str2%IntraFourier
      write(unitno,200)'Einstein ',':',str1%Einstein,str2%Einstein
      write(unitno,200)'Bond (J/mol) ',':',str1%Bond,str2%Bond
      write(unitno,200)'Angle (J/mol) ',':',str1%Angle,str2%Angle
      write(unitno,200)'Torsion (J/mol) ',':',str1%Torsion,str2%Torsion
      write(unitno,200)'IntraPair (J/mol) ',':',str1%IntraPair,str2%IntraPair
      write(unitno,200)'IntraCoul (J/mol) ',':',str1%IntraCoul,str2%IntraCoul
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,*)
    else
      write(unitno,*)
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,'(t30,a,t60,a)')trim(strname1),trim(strname2)
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,200)'Total(J/mol) ',':',str1%Total*UNIT_ENERGY*AVOGADRO_NUMBER,str2%Total*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,200)'van der Waals(J/mol) ',':',str1%vdW*UNIT_ENERGY*AVOGADRO_NUMBER,str2%vdW*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,200)'Real Coulombic(J/mol) ',':',str1%CoulReal*UNIT_ENERGY*AVOGADRO_NUMBER,str2%CoulReal*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,200)'Fourier Coulombic(J/mol) ',':',(str1%CoulFourier-str1%IntraFourier)*UNIT_ENERGY*AVOGADRO_NUMBER,(str2%CoulFourier-str2%IntraFourier)*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,200)'Einstein(J/mol) ',':',str1%Einstein*UNIT_ENERGY*AVOGADRO_NUMBER,str2%Einstein*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,200)'Bond (J/mol) ',':',str1%Bond*UNIT_ENERGY*AVOGADRO_NUMBER,str2%Bond*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,200)'Angle (J/mol) ',':',str1%Angle*UNIT_ENERGY*AVOGADRO_NUMBER,str2%Angle*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,200)'Torsion (J/mol) ',':',str1%Torsion*UNIT_ENERGY*AVOGADRO_NUMBER,str2%Torsion*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,200)'IntraPair (J/mol) ',':',str1%IntraPair*UNIT_ENERGY*AVOGADRO_NUMBER,str2%IntraPair*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,200)'IntraCoul (J/mol) ',':',str1%IntraCoul*UNIT_ENERGY*AVOGADRO_NUMBER,str2%IntraCoul*UNIT_ENERGY*AVOGADRO_NUMBER
      write(unitno,'(a)')'-----------------------------------------------------------------'
      write(unitno,*)
    end if
 
200 format(a,t30,a,e20.6,t60,e20.6)
  end subroutine CompareAndDisplayStorage

  !=============================================================================
end module storage
