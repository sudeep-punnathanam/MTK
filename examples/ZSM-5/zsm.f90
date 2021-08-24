program zsm
  implicit none
  integer, parameter :: pr=kind(1.D0)
  real(pr), parameter :: sx=20.07_pr,sy=19.92_pr,sz=13.42_pr
  real(pr), parameter :: lx=40.14_pr,ly=39.84_pr,lz=40.26_pr
  real(pr), dimension(3456) :: rx,ry,rz
  integer :: i,j,k,l
  integer :: nx,ny,nz

  open(unit=10,file='sili.mol',action='read')
  read(10,'(3/)')
  do i=1,288
    read(10,*)k,rx(i),ry(i),rz(i)
  end do

  k=0
  do nx=0,1
    do ny=0,1
      do nz=0,2
        do i=1,288
          rx(k+i)=rx(i)+real(nx,pr)*sx
          ry(k+i)=ry(i)+real(ny,pr)*sy
          rz(k+i)=rz(i)+real(nz,pr)*sz
        end do
        k=k+288
      end do
    end do
  end do
  close(unit=10)
  rx=rx-lx/2._pr
  ry=ry-ly/2._pr
  rz=rz-lz/2._pr

  open(unit=20,file='zsm5.xyz',action='write')
  open(unit=30,file='zsm5.config',action='write')
  open(unit=40,file='zsm5.mol',action='write')
  write(20,*)3456
  write(20,*)
  write(30,'(a10,i10)')'ZSM5',2304
  write(40,'(a10,i10)')'Group',2304
  write(40,*)'Rigid'
  l=0
  do j=1,12
    do i=1,288
      k=i+(j-1)*288
!!$      rx(k)=rx(k)-lx*anint(rx(k)/lx)
!!$      ry(k)=ry(k)-ly*anint(ry(k)/ly)
!!$      rz(k)=rz(k)-lz*anint(rz(k)/lz)
      if(i <= 96)then
        write(20,'(a4,3(4x,f10.4))')'Si',rx(k),ry(k),rz(k)
      else
        write(20,'(a4,3(4x,f10.4))')'O',rx(k),ry(k),rz(k)
        l=l+1
        write(30,100)1,l,rx(k),ry(k),rz(k)
        write(40,'(i10,a10,3f15.4)')l,'O',rx(k),ry(k),rz(k)
      end if
    end do
  end do
  close(unit=20)
  close(unit=30)
  close(unit=40)
100 format(i5,2x,i5,3(4x,f15.4))
end program zsm
