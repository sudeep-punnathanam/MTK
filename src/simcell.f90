module simcell
  use consts, only: PR, strlen
  use utils, only: Invert3x3Matrix, CROSS_PRODUCT, Matrix_x_Vector
  use random, only: RandomNumber
  implicit none
  public
  save
  
  !=============================================================================================================================
  !** Simulation Cell Definition
  type :: SimcellInfo
    real(PR), dimension(3)   :: BoxLength,BoxWidth,BoxAngle
    real(PR), dimension(3)   :: InverseBoxLength,HalfBoxLength
    real(PR), dimension(3,3) :: LatVecs,RLatVecs

    real(PR) :: Volume
    logical  :: NonOrthorhombic=.false.
  end type SimcellInfo
  !=============================================================================================================================

  interface ApplyBoundaryCondition
    module procedure :: ApplyBoundaryCondition_Atom
    module procedure :: ApplyBoundaryCondition_Molecule
  end interface ApplyBoundaryCondition

  interface GetMinimumImage
    module procedure :: GetMinimumImage_Atom
    module procedure :: GetMinimumImage_Molecule
  end interface GetMinimumImage


contains
  !=============================================================================================================================
  subroutine SetSimulationCellBoxLengthsAndBoxAngles(scell,boxlength,boxangle)
    type(SimcellInfo), intent(out) :: scell
    real(PR), dimension(3), intent(in) :: boxlength, boxangle

    scell%BoxLength=boxlength
    scell%BoxAngle=boxangle

    !** 1st-lattice vector
    scell%LatVecs(1,1)=1._PR
    scell%LatVecs(2,1)=0._PR
    scell%LatVecs(3,1)=0._PR

    !** 2nd-lattice vector
    scell%LatVecs(1,2)=cos(scell%BoxAngle(3))
    scell%LatVecs(2,2)=sin(scell%BoxAngle(3))
    scell%LatVecs(3,2)=0._PR

    !** 3rd-lattice vector
    scell%LatVecs(1,3)=cos(scell%BoxAngle(2))
    scell%LatVecs(2,3)=(cos(scell%BoxAngle(1))-cos(scell%BoxAngle(2))*cos(scell%BoxAngle(3))/sin(scell%BoxAngle(3)))
    scell%LatVecs(3,3)=sqrt(1._PR-scell%LatVecs(1,3)**2-scell%LatVecs(2,3)**2)

    !** Reciprocal lattice vectors
    scell%RLatVecs=Invert3x3Matrix(scell%LatVecs)

    !** Box Widths
    scell%BoxWidth=scell%BoxLength/sqrt(scell%RLatVecs(:,1)**2+scell%RLatVecs(:,2)**2+scell%RLatVecs(:,3)**2)

    !** Other variables
    scell%InverseBoxLength=1._PR/scell%BoxLength
    scell%HalfBoxLength=scell%BoxLength/2._PR

    !** Simulation Cell Volume
    scell%Volume=scell%BoxLength(1)*scell%BoxLength(2)*scell%BoxLength(3)
    if(scell%NonOrthorhombic)then
      scell%Volume = scell%Volume * dot_product(scell%LatVecs(:,1),CROSS_PRODUCT(scell%LatVecs(:,2),scell%LatVecs(:,3)))
    end if
  end subroutine SetSimulationCellBoxLengthsAndBoxAngles

  !=============================================================================================================================
  subroutine SetSimulationCellBoxLengths(scell,boxlength)
    type(SimcellInfo), intent(inout) :: scell
    real(PR), dimension(3), intent(in) :: boxlength

    scell%BoxLength=boxlength

    !** Box Widths
    scell%BoxWidth=scell%BoxLength/sqrt(scell%RLatVecs(:,1)**2+scell%RLatVecs(:,2)**2+scell%RLatVecs(:,3)**2)

    !** Other variables
    scell%InverseBoxLength=1._PR/scell%BoxLength
    scell%HalfBoxLength=scell%BoxLength/2._PR

    !** Simulation Cell Volume
    scell%Volume=scell%BoxLength(1)*scell%BoxLength(2)*scell%BoxLength(3)
    if(scell%NonOrthorhombic)then
      scell%Volume = scell%Volume * dot_product(scell%LatVecs(:,1),CROSS_PRODUCT(scell%LatVecs(:,2),scell%LatVecs(:,3)))
    end if
  end subroutine SetSimulationCellBoxLengths

  !=============================================================================================================================
  subroutine SetSimulationCellBoxAngles(scell,boxangle)
    type(SimcellInfo), intent(inout) :: scell
    real(PR), dimension(3), intent(in) :: boxangle

    scell%BoxAngle=boxangle

    !** 2nd-lattice vector
    scell%LatVecs(1,2)=cos(scell%BoxAngle(3))
    scell%LatVecs(2,2)=sin(scell%BoxAngle(3))
    scell%LatVecs(3,2)=0._PR

    !** 3rd-lattice vector
    scell%LatVecs(1,3)=cos(scell%BoxAngle(2))
    scell%LatVecs(2,3)=(cos(scell%BoxAngle(1))-cos(scell%BoxAngle(2))*cos(scell%BoxAngle(3))/sin(scell%BoxAngle(3)))
    scell%LatVecs(3,3)=sqrt(1._PR-scell%LatVecs(1,3)**2-scell%LatVecs(2,3)**2)

    !** Reciprocal lattice vectors
    scell%RLatVecs=Invert3x3Matrix(scell%LatVecs)

    !** Box Widths
    scell%BoxWidth=scell%BoxLength/sqrt(scell%RLatVecs(:,1)**2+scell%RLatVecs(:,2)**2+scell%RLatVecs(:,3)**2)

    !** Simulation Cell Volume
    !** Change of Box angles assumes that the simulation cell is non-orthorhombic
    scell%Volume=scell%BoxLength(1)*scell%BoxLength(2)*scell%BoxLength(3)
    scell%Volume = scell%Volume * dot_product(scell%LatVecs(:,1),CROSS_PRODUCT(scell%LatVecs(:,2),scell%LatVecs(:,3)))

  end subroutine SetSimulationCellBoxAngles

  !=============================================================================================================================
  subroutine ApplyBoundaryCondition_Atom(scell,atvec)
    type(SimcellInfo), intent(in)         :: scell
    real(PR), dimension(3), intent(inout) :: atvec

    integer :: i
    real(PR) :: dX,dY,dZ

    if(scell%NonOrthorhombic)atvec=Matrix_x_Vector(scell%RLatVecs,atvec)

    atvec(1)=atvec(1)-scell%BoxLength(1)*real(floor(atvec(1)*scell%InverseBoxLength(1)+0.5_PR),PR)
    atvec(2)=atvec(2)-scell%BoxLength(2)*real(floor(atvec(2)*scell%InverseBoxLength(2)+0.5_PR),PR)
    atvec(3)=atvec(3)-scell%BoxLength(3)*real(floor(atvec(3)*scell%InverseBoxLength(3)+0.5_PR),PR)

    if(scell%NonOrthorhombic)atvec=Matrix_x_Vector(scell%LatVecs,atvec)

  end subroutine ApplyBoundaryCondition_Atom

  !=============================================================================================================================
  subroutine ApplyBoundaryCondition_Molecule(scell,natoms,atvec,com)
    type(SimcellInfo), intent(in)                :: scell
    integer, intent(in)                          :: natoms
    real(PR), dimension(3,natoms), intent(inout) :: atvec
    real(PR), dimension(3), intent(inout)        :: com

    integer :: i
    real(PR) :: dX,dY,dZ

    if(scell%NonOrthorhombic)then
      do i=1,natoms
        atvec(:,i)=Matrix_x_Vector(scell%RLatVecs,atvec(:,i))
      end do
      com=Matrix_x_Vector(scell%RLatVecs,com)
    end if

    dX=scell%BoxLength(1)*real(floor(com(1)*scell%InverseBoxLength(1)+0.5_PR),PR)
    dY=scell%BoxLength(2)*real(floor(com(2)*scell%InverseBoxLength(2)+0.5_PR),PR)
    dZ=scell%BoxLength(3)*real(floor(com(3)*scell%InverseBoxLength(3)+0.5_PR),PR)

    do i=1,natoms
      atvec(1,i)=atvec(1,i)-dX
      atvec(2,i)=atvec(2,i)-dY
      atvec(3,i)=atvec(3,i)-dZ
    end do
    com(1)=com(1)-dX
    com(2)=com(2)-dY
    com(3)=com(3)-dZ

    if(scell%NonOrthorhombic)then
      do i=1,natoms
        atvec(:,i)=Matrix_x_Vector(scell%LatVecs,atvec(:,i))
      end do
      com=Matrix_x_Vector(scell%LatVecs,com)
    end if

  end subroutine ApplyBoundaryCondition_Molecule

  !=============================================================================================================================
  subroutine GetMinimumImage_Atom(scell,n,sepvec,distsq)
    type(SimcellInfo), intent(in)           :: scell
    integer, intent(in)                     :: n
    real(PR), dimension(3,n), intent(inout) :: sepvec
    real(PR), dimension(n), intent(out)     :: distsq

    integer :: i
    real(PR) :: dX,dY,dZ

    if(scell%NonOrthorhombic)then
      do i=1,n
        sepvec(1:3,i)=Matrix_x_Vector(scell%RLatVecs,sepvec(1:3,i))
      end do
    end if


    do i=1,n
      sepvec(1,i)=sepvec(1,i)-scell%BoxLength(1)*real(floor(sepvec(1,i)*scell%InverseBoxLength(1)+0.5_PR),PR)
      sepvec(2,i)=sepvec(2,i)-scell%BoxLength(2)*real(floor(sepvec(2,i)*scell%InverseBoxLength(2)+0.5_PR),PR)
      sepvec(3,i)=sepvec(3,i)-scell%BoxLength(3)*real(floor(sepvec(3,i)*scell%InverseBoxLength(3)+0.5_PR),PR)
    end do

    if(scell%NonOrthorhombic)then
      do i=1,n
        sepvec(1:3,i)=Matrix_x_Vector(scell%LatVecs,sepvec(1:3,i))
      end do
    end if

    do i=1,n
      distsq(i)=sepvec(1,i)**2+sepvec(2,i)**2+sepvec(3,i)**2
    end do

  end subroutine  GetMinimumImage_Atom

  !=============================================================================================================================
  subroutine GetMinimumImage_Molecule(scell,n,sepvec,distsq,comsepvec,vfactor)
    type(SimcellInfo), intent(in)           :: scell
    integer, intent(in)                     :: n
    real(PR), dimension(3,n), intent(inout) :: sepvec,comsepvec
    real(PR), dimension(n), intent(out)     :: distsq,vfactor

    integer :: i
    real(PR) :: dX,dY,dZ

    if(scell%NonOrthorhombic)then
      do i=1,n
        sepvec(1:3,i)=Matrix_x_Vector(scell%RLatVecs,sepvec(1:3,i))
        comsepvec(1:3,i)=Matrix_x_Vector(scell%RLatVecs,comsepvec(1:3,i))
      end do
    end if


    do i=1,n
      dX=scell%BoxLength(1)*real(floor(sepvec(1,i)*scell%InverseBoxLength(1)+0.5_PR),PR)
      dY=scell%BoxLength(2)*real(floor(sepvec(2,i)*scell%InverseBoxLength(2)+0.5_PR),PR)
      dZ=scell%BoxLength(3)*real(floor(sepvec(3,i)*scell%InverseBoxLength(3)+0.5_PR),PR)

      sepvec(1,i)=sepvec(1,i)-dX
      sepvec(2,i)=sepvec(2,i)-dY
      sepvec(3,i)=sepvec(3,i)-dZ
      comsepvec(1,i)=comsepvec(1,i)-dX
      comsepvec(2,i)=comsepvec(2,i)-dY
      comsepvec(3,i)=comsepvec(3,i)-dZ
    end do

    if(scell%NonOrthorhombic)then
      do i=1,n
        sepvec(1:3,i)=Matrix_x_Vector(scell%LatVecs,sepvec(1:3,i))
        comsepvec(1:3,i)=Matrix_x_Vector(scell%LatVecs,comsepvec(1:3,i))
      end do
    end if

    do i=1,n
      distsq(i)=sepvec(1,i)**2+sepvec(2,i)**2+sepvec(3,i)**2
      vfactor(i)=sepvec(1,i)*comsepvec(1,i)+sepvec(2,i)*comsepvec(2,i)+sepvec(3,i)*comsepvec(3,i)
      vfactor(i)=vfactor(i)/distsq(i)
    end do

  end subroutine  GetMinimumImage_Molecule

  !=============================================================================================================================
  subroutine RandomPointInCell(scell,vec)
    type(SimcellInfo), intent(in)       :: scell
    real(PR), dimension(3), intent(OUT) :: vec

    vec(1)=(RandomNumber()-0.5_PR)*scell%BoxLength(1)
    vec(2)=(RandomNumber()-0.5_PR)*scell%BoxLength(2)
    vec(3)=(RandomNumber()-0.5_PR)*scell%BoxLength(3)

    if(scell%NonOrthorhombic)vec=matmul(scell%LatVecs,vec)

  end subroutine RandomPointInCell

  !=============================================================================================================================
end module simcell
