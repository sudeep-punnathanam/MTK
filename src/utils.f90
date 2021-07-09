module utils
  use consts, only: PR, PI, TWOPI, strlen, lstrlen, xlstrlen
  implicit none
  public

  interface swap
    module procedure swap_i, swap_r
  end interface

  character(len=xlstrlen) :: ErrorMessage=''

contains
  subroutine ReadStringFromFile(String,Tag,InputFile,filename,lineno,EndTag,error)
    character(len=lstrlen), intent(out) :: String
    character(len=*), intent(in) :: Tag
    character(len=lstrlen), dimension(:), intent(in) :: InputFile
    character(len=*), intent(in) :: filename
    integer, intent(out) :: lineno
    character(len=*), intent(in), optional :: EndTag
    logical, intent(out), optional :: error

    integer :: iend

    if(present(error))error=.false.

    iend=size(InputFile)
    do lineno=1,iend
      if(ReadString(InputFile(lineno),trim(Tag),String)==1)return
      if(present(EndTag))then
        if(ReadString(InputFile(lineno),trim(EndTag),String)==1)exit
      end if
    end do
    lineno=0

    if(present(error))then
      write(ErrorMessage,'(2a,i5,4x,4a)')__FILE__,':',__LINE__, &
        'Tag ',trim(adjustl(Tag)),' missing in ',trim(filename)
      error=.true.
    end if
      
  end subroutine ReadStringFromFile

  function ReadString(line,tag,string) result(pos)
    character(LEN=*), intent(IN) :: line
    character(LEN=*), intent(IN) :: tag
    character(LEN=*), intent(OUT) :: string
    integer :: pos

    integer :: n
    character(LEN=xlstrlen) :: lline
    character(LEN=strlen) :: ltag
    character(LEN=strlen) :: fmtstr

    lline=' '
    ltag=' '
    lline=lowercase(line)
    ltag=lowercase(tag)
    ltag=adjustl(ltag)
    pos=index(lline,trim(ltag))
    if(pos == 1)then
      n=len_trim(adjustl(ltag))
      write(fmtstr,'(a,i2,a)')"(",n,"x,a)"
      read(line,fmtstr)string
      if(index(string,' ') == 1)then
        string=trim(adjustl(string))
      else
        pos=0
      end if
    else
      string=' '
    end if 
  end function ReadString

  function split(line,fields,opt_delim) result(nfields)
    character(len=*), intent(In) :: line
    character(len=*), dimension(:), intent(Out) :: fields
    character, intent(In), optional :: opt_delim
    integer                         :: nfields

    character(len=xlstrlen)         :: tr_line
    character                       :: delim
    integer                         :: length,i,j

    if(present(opt_delim))then
      delim=opt_delim
    else
      delim=' '
    end if

    tr_line=trim(adjustl(line))
    nfields=0
    length=len_trim(tr_line)
    if(length == 0)return
    if(length == len(tr_line))then
      write(0,'(2a,i4,a)')__FILE__,':',__LINE__,'Increase length of line'
      stop
    end if

    i=1
    nfields=0
    do
      j=index(tr_line(i:),delim)
      if(j /= 1)then
        nfields=nfields+1
        if(nfields > size(fields))then
          write(0,'(2a,i4,2a)')__FILE__,':',__LINE__,'Too many fields in line :',trim(tr_line)
          stop
        end if
        fields(nfields)=tr_line(i:i+j-2)
      end if
      i=i+j
      if(i > length)exit
    end do
  end function split

  subroutine FindStringInFile(tag,unitno,lineno,string,error)
    character(len=*), intent(In)     :: tag
    integer, intent(In)              :: unitno
    integer, intent(Out)             :: lineno
    character(lstrlen), intent(Out)  :: string
    logical, intent(Out)             :: error

    character(len=lstrlen)       :: line
    character(len=xlstrlen)      :: filename
    integer                      :: indx, ierror
    logical                      :: lopen

    error=.false.
    inquire(unit=unitno,name=filename,opened=lopen,iostat=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,a,i4)')__FILE__,':',__LINE__, &
        'No file associated with unit number ',unitno
      error=.true.
      return
    end if
    if(.not. lopen)then
      write(ErrorMessage,'(2a,i5,4x,3a)')__FILE__,':',__LINE__, &
        'File ',trim(filename),' not open.'
      error=.true.
      return
    end if

    rewind(unitno)
    lineno = 0    
    do 
      lineno = lineno + 1
      read(unitno,'(a)', IOSTAT=ierror)line
      line=adjustl(line)

      if (ierror == -1) then
        !** End-of-file
        lineno = 0
        line = ' '
        exit
      else if(ierror /= 0)then
        write(ErrorMessage,'(2a,i5,4x,4a)')__FILE__,':',__LINE__, &
          'Read failure while searching for string ',trim(tag),' in file ',trim(filename)
        error=.true.
        return
      end if
      
      indx=ReadString(line,tag,string)
      if(indx == 1)return
      
    end do
  end subroutine FindStringInFile

  pure function lowercase(string1)  result(string2)
    character(len=*), intent(in)  :: string1
    character(len=:), allocatable :: string2

    integer :: i,length

    length=len(string1)
    allocate(character(len=length) :: string2)
    do i=1,length
      if(lge(string1(i:i),'A') .and. lle(string1(i:i),'Z'))then
        string2(i:i)=achar(iachar(string1(i:i)) + 32 )
      else
        string2(i:i)=string1(i:i)
      end if
    end do
  end function lowercase

  pure function uppercase(string1)  result(string2)
    character(len=*), intent(In)  :: string1
    character(len=:), allocatable :: string2

    integer :: i,length

    length=len(string1)
    allocate(character(len=length) :: string2)
    do i=1,length
      if(lge(string1(i:i),'a') .and. lle(string1(i:i),'z'))then
        string2(i:i)=achar(iachar(string1(i:i)) - 32 )
      else
        string2(i:i)=string1(i:i)
      end if
    end do
  end function uppercase

  !==============================================================================================================
  !**   Rotate Vectors
  !==============================================================================================================
  subroutine RotateVectors(Positions,RotationMatrix,RotatePoint)
    real(PR), dimension(:,:), intent(inout) :: Positions
    real(PR), dimension(3,3), intent(IN) :: RotationMatrix
    real(PR), dimension(3), intent(IN) :: RotatePoint

    real(PR), dimension(3) :: vec
    integer :: a,n
    
    n=size(Positions,2)

    !** Rotate atoms of the molecule
    do a=1,n
      vec=Positions(:,a)-RotatePoint
      vec=Matrix_x_Vector(RotationMatrix,vec)
      Positions(:,a)=RotatePoint+vec
    end do

  end subroutine RotateVectors

  !==========================================================================
  !** Gets a unit number for a file.
  !==========================================================================
  integer function GetFileUnit(filename,lopen,error)
    character(len=*), intent(In) :: filename
    logical, intent(Out)         :: lopen
    logical, intent(Out)         :: error

    integer :: unitno,ierror

    error=.false.

    inquire(FILE=trim(filename),OPENED=lopen,NUMBER=unitno,IOSTAT=ierror)
    if(ierror /= 0)then
      write(ErrorMessage,'(2a,i5,4x,3a,i3)')__FILE__,':',__LINE__,&
        'INQUIRE statement failure for file ',trim(filename),'. IOSTAT =',ierror
      error=.true.
      return
    end if
    if(lopen)then
      GetFileUnit=unitno
      return
    end if
    unitno=10
    do
      inquire(UNIT=unitno,OPENED=lopen,IOSTAT=ierror)
      if(ierror /= 0)then
        write(ErrorMessage,'(2a,i5,4x,a,i3,a,i3)')__FILE__,':',__LINE__,&
          'INQUIRE statement failure for unit number ',unitno,'. IOSTAT =',ierror
        error=.true.
        return
      end if
      if(.not. lopen)exit
      unitno=unitno+1
    end do
    GetFileUnit=unitno
  end function GetFileUnit

  !============================================================================
  !** Subset Number
  !============================================================================
  integer function GetSubsetNumber(spc1,spc2)
    integer, intent(In) :: spc1,spc2

    integer :: n,m

    n=max(spc1,spc2)
    m=min(spc1,spc2)
    GetSubsetNumber=(n*(n+1))/2-m+1
  end function GetSubsetNumber 

  !============================================================================
  ! Determines the format specification for diplaying a real number with n 
  ! significant digits (default is 6)
  !============================================================================
  function NumberFormat(a,n) result(fmtstr)
    real(PR), intent(IN) :: a
    integer, intent(in), optional :: n
    character(LEN=strlen)   :: fmtstr

    real(PR) :: b
    integer :: i,m

    if(present(n))then
      m=n
    else
      m=6
    end if
    b=log10(abs(a))
    if(b < real(m,pr) .and. b > -1._PR)then
      i=m-ceiling(b)
      if(i<10)then
        write(fmtstr,'(a,i1,a)')'(f15.',i,')'
      else
        write(fmtstr,'(a,i2,a)')'(f15.',i,')'
      end if
    else
      fmtstr='(es15.5e3)'
    end if
  end function NumberFormat

  !============================================================================
  ! Convert a real number to a string
  !============================================================================
  function RealToString(a,n)
    real(PR), intent(in) :: a
    integer, intent(in), optional :: n
    character(len=strlen) :: RealToString

    character(len=strlen) :: fmtstr

    if(present(n))then
      fmtstr=NumberFormat(a,n)
    else
      fmtstr=NumberFormat(a)
    end if

    write(RealToString,fmtstr)a

  end function RealToString

  !============================================================================
  ! Inverts a 3x3 matrix numerically
  !============================================================================
  function Invert3x3Matrix(a) result(b)
    real(PR), dimension(3,3), intent(In) :: a
    real(PR), dimension(3,3)             :: b

    real(PR) :: d

    !**Calculate the adjoint matrix
    b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
    b(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
    b(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1)
    b(1,2)=a(3,2)*a(1,3)-a(3,3)*a(1,2)
    b(2,2)=a(3,3)*a(1,1)-a(3,1)*a(1,3)
    b(3,2)=a(3,1)*a(1,2)-a(3,2)*a(1,1)
    b(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
    b(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
    b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)

    !** Calculate the determinant
    d=a(1,1)*b(1,1) + a(1,2)*b(2,1) + a(1,3)*b(3,1)

    !** Obtain the inverse
    b=b/d
    
  end function Invert3x3Matrix

  !============================================================================
  ! Multiply 3x3 matrix with vector
  !============================================================================
  function Matrix_x_Vector(a,b) result(c)
    real(PR), dimension(3,3), intent(In) :: a
    real(PR), dimension(3), intent(In)   :: b
    real(PR), dimension(3)               :: c

    c(1)=a(1,1)*b(1)+a(1,2)*b(2)+a(1,3)*b(3)
    c(2)=a(2,1)*b(1)+a(2,2)*b(2)+a(2,3)*b(3)
    c(3)=a(3,1)*b(1)+a(3,2)*b(2)+a(3,3)*b(3)

  end function Matrix_x_Vector

  !============================================================================
  ! Multiply vector with 3x3 matrix
  !============================================================================
  function Vector_x_Matrix(b,a) result(c)
    real(PR), dimension(3), intent(In)   :: b
    real(PR), dimension(3,3), intent(In) :: a
    real(PR), dimension(3)               :: c

    c(1)=b(1)*a(1,1)+b(2)*a(2,1)+b(3)*a(3,1)
    c(2)=b(1)*a(1,2)+b(2)*a(2,2)+b(3)*a(3,2)
    c(3)=b(1)*a(1,3)+b(2)*a(2,3)+b(3)*a(3,3)

  end function Vector_x_Matrix

  !============================================================================
  ! Multiply 3x3 matrix with 3x3 matrix
  !============================================================================
  function Matrix_x_Matrix(a,b) result(c)
    real(PR), dimension(3,3), intent(In) :: a,b
    real(PR), dimension(3,3)             :: c

    c(1,1)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
    c(2,1)=a(2,1)*b(1,1)+a(2,2)*b(2,1)+a(2,3)*b(3,1)
    c(3,1)=a(3,1)*b(1,1)+a(3,2)*b(2,1)+a(3,3)*b(3,1)
    c(1,2)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
    c(2,2)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
    c(3,2)=a(3,1)*b(1,2)+a(3,2)*b(2,2)+a(3,3)*b(3,2)
    c(1,3)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
    c(2,3)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
    c(3,3)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)

  end function Matrix_x_Matrix

  !=========================================================================================
  ! Calculates the cross product of two vectors
  !=========================================================================================
  pure function CROSS_PRODUCT(array1,array2) result(array)
    real(Kind=PR), dimension(3), intent(In) :: array1,array2
    real(Kind=PR), dimension(3)             :: array

    array(1)=array1(2)*array2(3)- array1(3)*array2(2)
    array(2)=array1(3)*array2(1)- array1(1)*array2(3)
    array(3)=array1(1)*array2(2)- array1(2)*array2(1)
  end function CROSS_PRODUCT

  !=========================================================================================
  ! Determines the coefficients and points for Gauss-Legendre Quadrature
  !=========================================================================================
  subroutine gauleg(x1,x2,x,w)
    implicit none
    real(PR), parameter :: pi=3.141592653589793238462643383279502884197_PR
    real(PR), intent(in) :: x1,x2
    real(PR), dimension(:), intent(out) :: x,w
    real(PR), parameter :: EPS=3.0e-14_PR

    !** Given the lower and upper limits of integration x1 and x2, this routine returns arrays x and w
    !** of length N containing the abscissas and weights of the Gauss-Legendre N-point quadrature
    !** formula. The parameter EPS is the relative precision. Note that internal computations are
    !** done in double precision.

    integer :: its,j,m,n
    integer, parameter :: MAXIT=10
    real(PR) :: xl,xm
    real(PR), dimension((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
    logical, dimension((size(x)+1)/2) :: unfinished

    n=size(x)
    m=(n+1)/2

    xm=0.5_PR*(x2+x1)
    xl=0.5_PR*(x2-x1)
    z=cos(pi*(arth(1,1,m)-0.25_PR)/(n+0.5_PR))
    unfinished=.true.
    do its=1,MAXIT
      where(unfinished)
        p1=1._PR
        p2=0._PR
      end where
      do j=1,n
        where(unfinished)
          p3=p2
          p2=p1
          p1=((2._PR*j-1._PR)*z*p2-(j-1._PR)*p3)/j
        end where
      end do
      !** p1 now contains the desired Legendre polynomials. We next compute pp, the derivatives,
      !** by a standard relation involving also p2, the polynomials of one lower order.
      where(unfinished)
        pp=n*(z*p1-p2)/(z*z-1._PR)
        z1=z
        z=z1-p1/pp
        unfinished=(abs(z-z1) > EPS)
      end where
      if(.not. any(unfinished))exit
    end do

    if(its == MAXIT+1)then
      write(*,*)'too many iterations in gauleg'
      stop
    endif
    x(1:m)=xm-xl*z
    x(n:n-m+1:-1)=xm+xl*z
    w(1:m)=2._PR*xl/((1._pr-z**2)*pp**2)
    w(n:n-m+1:-1)=w(1:m)
    return
  contains
    !** Array function returning an arithmetic progression
    function arth(first,increment,n)
      integer, intent(in) :: first,increment,n
      integer, dimension(n) :: arth
      integer :: k,k2,temp
      integer, parameter :: NPAR_ARTH=16,NPAR2_ARTH=8

      if(n > 0)arth(1)=first
      if(n <= NPAR_ARTH)then
        do k=2,n
          arth(k)=arth(k-1)+increment
        enddo
      else
        do k=2,NPAR2_ARTH
          arth(k)=arth(k-1)+increment
        enddo
        temp=increment*NPAR2_ARTH
        k=NPAR2_ARTH
        do
          if(k >= n)exit
          k2=k+k
          arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
          temp=temp+temp
          k=k2
        enddo
      endif
    end function arth
  end subroutine gauleg

  pure real(PR) function ErrorFunctionComplement(x) result(y)
    implicit none
    real(PR), intent(in) :: x

    REAL(PR) :: t,u
    REAL(PR), parameter ::  PA=3.97886080735226000_PR,P0=2.75374741597376782e-1_PR,P1=4.90165080585318424e-1_PR
    REAL(PR), parameter ::  P2=7.74368199119538609e-1_PR,P3=1.07925515155856677_PR,P4=1.31314653831023098_PR
    REAL(PR), parameter ::  P5=1.37040217682338167_PR,P6=1.18902982909273333_PR,P7=8.05276408752910567e-1_PR
    REAL(PR), parameter ::  P8=3.57524274449531043e-1_PR,P9=1.66207924969367356e-2_PR,P10=-1.19463959964325415e-1_PR
    REAL(PR), parameter ::  P11=-8.38864557023001992e-2_PR,P12=2.49367200053503304e-3_PR,P13=3.90976845588484035e-2_PR
    REAL(PR), parameter ::  P14=1.61315329733252248e-2_PR,P15=-1.33823644533460069e-2_PR,P16=-1.27223813782122755e-2_PR
    REAL(PR), parameter ::  P17=3.83335126264887303e-3_PR,P18=7.73672528313526668e-3_PR,P19=-8.70779635317295828e-4_PR
    REAL(PR), parameter ::  P20=-3.96385097360513500e-3_PR,P21=1.19314022838340944e-4_PR,P22=1.27109764952614092e-3_PR

    t=PA/(PA+abs(x))
    u=t-0.5_PR
    y=(((((((((P22*u+P21)*u+P20)*u+P19)*u+P18)*u+P17)*u+P16)*u+P15)*u+P14)*u+P13)*u+P12
    y=((((((((((((y*u+P11)*u+P10)*u+P9)*u+P8)*u+P7)*u+P6)*u+P5)*u+P4)*u+P3)*u+P2)*u+P1)*u+P0)*t*exp(-x*x)
    if(x<0._PR) y=2._PR-y

  end function ErrorFunctionComplement


  
  !** Swaps integers
  subroutine swap_i(a,b)
    integer, intent(inout) :: a,b
    integer :: temp
    temp=a
    a=b
    b=temp
  end subroutine swap_i

  !** Swaps real numbers
  subroutine swap_r(a,b)
    real(PR), intent(inout) :: a,b
    real(PR) :: temp
    temp=a
    a=b
    b=temp
  end subroutine swap_r
  
end module utils
