! Haoyu Shang Peking University
! 2023.10.12
! Email shy@stu.pku.edu.cn

module para
! this module contains the parameters used in this code   
    implicit none

    real(kind=8) ::C=1000.0d0
    real(kind=8) :: C2p=0.005d0 !C2 is used in radius calculation, C2=C2p*C
    real(kind=8) ::pi=3.141592653589793d0
    integer :: gaussnum=735
    real(kind=8) :: masshat=938.91852d0
    real(kind=8) :: massbar=938.91897d0
    real(kind=8) :: hbarc=197.327053d0
end module


module gauss_legr
      ! this module set the guass legr points
          private
          public gauleg,zerotoinf
          contains
            
       
            SUBROUTINE gauleg(x1,x2,x,w,n)  
            INTEGER n 
            integer, parameter :: rk = kind ( 1.0D+00 )
            real ( kind = rk ) x1,x2,x(n),w(n)  
            DOUBLE PRECISION EPS  
            PARAMETER (EPS=3.d-14)  
            INTEGER i,j,m  
            DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1  
            m=(n+1)/2  
            xm=0.5d0*(x2+x1)  
            xl=0.5d0*(x2-x1)  
            do  i=1,m  
              z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))  
      1       continue  
                p1=1.d0  
                p2=0.d0  
                do j=1,n  
                  p3=p2  
                  p2=p1  
                  p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j  
               end do  
                pp=n*(z*p1-p2)/(z*z-1.d0)  
                z1=z  
                z=z1-p1/pp  
              if(abs(z-z1).gt.EPS)goto 1  
              x(i)=xm-xl*z  
              x(n+1-i)=xm+xl*z  
              w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)  
              w(n+1-i)=w(i)  
            end do  
            return  
            end subroutine
      
            subroutine zerotoinf(n,k,s)
                  use para,only:C,pi
                  implicit none
                  real(kind=8) :: k(n),s(n)
                  ! k is x point ,s is weight
                  integer :: n
                  ! n is the number we set the guass Legendre
                  real(kind=8) :: x(n),w(n)
                  ! x and w are points and weights in [-1,1]
                  integer ::i
                  x=0.000000000000000d0
                  w=0.000000000000000d0
                  k=0.000000000000000d0
                  s=0.000000000000000d0
                  call gauleg(-1.0d0,1.0d0,x,w,n)
                  do i=1,n
                        k(i)=C*tan(pi*(x(i)+1.0d0)/4.0d0)
                        s(i)=C*pi/4.0d0*w(i)/(cos(pi*(x(i)+1.0d0)/4.0d0)**2)
                  end do
                  end
      end module

module sphericalbessel
      implicit none  
      public sphericalbessel0,sphericalbessel2  
      contains
      function sphericalbessel0(x)
            implicit none
            ! input
            real(kind=8) :: x
            ! output
            real(kind=8) ::sphericalbessel0

            sphericalbessel0=sin(x)/x

            return
      end function
      function sphericalbessel2(x)
            implicit none
            ! input
            real(kind=8) :: x
            ! output
            real(kind=8) ::sphericalbessel2

            sphericalbessel2=(3.0d0/x**2-1.0d0)*sin(x)/x-3.0d0*cos(x)/x**2

            return
      end function

      
end module 
      !Sort of bubble method
subroutine index_sort(array, n, indices)
      implicit none
      integer, intent(in) :: n
      real(kind=8), intent(in) :: array(n)
      integer, intent(out) :: indices(n)
      integer :: i, j, temp_index
      
      indices = [(i, i=1, n)]  ! initialize indices
      
      do i = 1, n-1
            do j = 1, n-i
                  if (array(indices(j)) > array(indices(j+1))) then
                  temp_index = indices(j)
                  indices(j) = indices(j+1)
                  indices(j+1) = temp_index
                  end if
            end do
      end do
      end subroutine index_sort

module wavefunction
      implicit none
      
      contains
      function ufunction(r,k,s,psi,n)
            use para,only:pi,hbarc
            use sphericalbessel
            implicit none 
            !input
            real(kind=8) :: r
            ! array k(n),s(n),gauss-legendre points
            integer ::n
            real(kind=8),dimension(n) :: k,s
            real(kind=8),dimension(n) :: psi !the momentum phase 
            !output
            real(kind=8) ::ufunction

            integer ::i

            ufunction=0.0d0
            do i=1,n
                  ufunction=ufunction+k(i)**2*sphericalbessel0((k(i)*r)/hbarc)*psi(i)*hbarc**(-1.5d0)*s(i)*r*dsqrt(2.0d0/pi)
            end do
            !if (r .gt. 20.0d0) ufunction=0.0d0
            return
      end function
      function wfunction(r,k,s,psi,n)
            use para,only:pi,hbarc
            use sphericalbessel
            implicit none 
            !input
            real(kind=8) :: r
            ! array k(n),s(n),gauss-legendre points
            integer ::n
            real(kind=8),dimension(n) :: k,s
            real(kind=8),dimension(n) :: psi !the momentum phase 
            !output
            real(kind=8) ::wfunction

            integer ::i
            
            wfunction=0.0d0
            do i=1,n
                  wfunction=wfunction+k(i)**2*sphericalbessel2((k(i)*r)/hbarc)*psi(i)*hbarc**(-1.5d0)*s(i)*r*dsqrt(2.0d0/pi)
            end do
           !if (r .gt. 20.0d0) wfunction=0.0d0

            return
      end function
      
end module wavefunction
program deuteron
      use gauss_legr
      use para
      use sphericalbessel
      use wavefunction
      implicit none

! these are varables used in dsyev
      real(kind=8),allocatable ::A(:,:),work(:),evector(:,:)
      real(kind=8),allocatable :: evr(:),evi(:) !the eigenvalues
      integer :: info,lwork,dim

!these are varables used in infinite integral
      real(kind=8),allocatable :: k(:),s(:)

!these are varables used to get potential elements
      external n3lo500new
      real(kind=8) :: v(6),xmev,ymev
      integer :: j,inn
      logical heform,sing,trip,coup,endep
      character*4 label
      common /cpot/   v,xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep,label
      common /cnn/ inn
      real(kind=8) :: momentumi,momentumj
      
!file inout varables
      integer ::kread,kwrite,kwrite1,kwrite2,kwrite3,ios     
      character(len=20) :: format_str

!result varables
      real(kind=8) :: bindingenergy
      real(kind=8),allocatable :: eigenvector(:)
      !3S1 u(r)
      real(kind=8) :: r 
      real(kind=8) ::u_r(100) 
      real(kind=8) ::asymptotic_u_r(100) 
      !3D1 w(r)
      real(kind=8) ::w_r(100)
      real(kind=8) ::asymptotic_w_r1(100) 
      real(kind=8) ::asymptotic_w_r2(100) 
      !radius
      real(kind=8) ::radius,radius_temp
      real(kind=8),allocatable ::sr(:),kr(:)
      !quadrupole moment
      real(kind=8) ::quadrupolem
      !asymptotic As
      real(kind=8) :: asymptotic_As
      !asymptotic Ad
      real(kind=8) :: asymptotic_Ad
      ! gamma
      real(kind=8) :: gamma
      !ete
      real(kind=8) :: eta
      !pd
      real(kind=8) :: pd,ps
!local varables      
      integer :: i,row,column
      real(kind=8) :: delta
      integer,allocatable :: indices(:) !sort by indices
      real(kind=8) ::norm
      real(kind=8) ::avesum
      integer ::count
      real(kind=8) :: uf,wf

!open file 
      kwrite=10
      open(unit=10,file='deuteron_output')     
      kwrite3=13
      open(unit=13,file='deuteron_output_detail') 

      
! set the gauss-Legendre points
      allocate(k(gaussnum))
      allocate(s(gaussnum))
      call zerotoinf(gaussnum,k,s)
      !write the ki,kj points
      write(kwrite3,*) 'ki,kj points:'
      do i=1,gaussnum
            write(kwrite3,"(F20.10)") k(i)
      end do
      
! set the potential values
      heform=.false.
      sing=.true.
      trip=.true.
      coup=.true.
      inn=2  !np interaction
      j=1 ! 3S1 and 3D1

!set the matrix elements     
      dim=2*gaussnum 
      !create format_str
      write(format_str, "('(',I0,'F20.10)')") dim
      allocate(A(dim,dim))
      allocate(evector(dim,dim))
      allocate(evr(dim))
      allocate(evi(dim))
      A=0.0d0
      evi=0.0d0
      evr=0.0d0
      do row=1,gaussnum
            do column=1,gaussnum
                  momentumi=k(row)
                  momentumj=k(column)
                  xmev=momentumi
                  ymev=momentumj
                  call n3lo500new
                  !v(4)=v00(i,j) v(6)=v02(i,j) v(5)=v20(i,j) v(3)=22(i,j)
                  !delta_ij
                  if (row .eq. column)then
                        delta=1.0d0
                  else
                        delta=0.0d0
                  end if

                  !top-left block
                  A(row,column)=A(row,column)+v(4)*momentumj**2*s(column)+delta*momentumi**2/masshat
                  !top-right block
                  A(row,gaussnum+column)=A(row,gaussnum+column)+v(6)*momentumj**2*s(column)
                  !lower-left block
                  A(gaussnum+row,column)=A(gaussnum+row,column)+v(5)*momentumj**2*s(column)
                  !lower-right block
                  A(gaussnum+row,gaussnum+column)=A(gaussnum+row,gaussnum+column)+v(3)*momentumj**2*s(column)+delta*momentumi**2/masshat

            end do
      end do
!write the matrix elements
!      write(kwrite,*) 'matrix elements:'
!      do row=1,dim
!            write(kwrite,format_str) (A(row,i),i=1,dim)       
!      end do

!solve the matrix's eigenvalues

      !initial dsyev
      lwork = -1
      allocate(work(1))
      call dgeev('N', 'V', dim, A, dim, evr, evi, A, dim, evector, dim, work, lwork, info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))

      !compute the eigenvalues and eigenvectors 
      call dgeev('N', 'V', dim, A, dim, evr, evi, A, dim, evector, dim, work, lwork, info)
      !write the results
!      write(kwrite,*) 'eigenvalues-realpart:'
!      write(kwrite,format_str) evr
      
      !Order from small to large
      allocate(indices(dim))
      call index_sort(evr, dim, indices)
      write(kwrite3,*) 'eigenvalues-realpart-ordered:'
      write(kwrite3,format_str) (evr(indices(i)),i=1,dim)
!      write(kwrite,*) 'eigenvalues-imaginarypart:'
!      write(kwrite,format_str) evi
!      write(kwrite,*) 'eigenvectors:'
!      do row=1,dim
!            write(kwrite,format_str) (evector(row,i),i=1,dim)       
!      end do  
      
      !change to binding energy
      bindingenergy=2.0d0*massbar-dsqrt((4.0d0*massbar)**2+16.0d0*massbar*evr(indices(1)))/2.0d0
      write(kwrite,*) 'bindingenergy:'
      write(kwrite,*) bindingenergy

      !the wave function in momentum phase
      allocate(eigenvector(dim))
      eigenvector=evector(:,indices(1))
      ! write(kwrite3,*) 'eigenvector in momentum phase:'
      ! write(kwrite3,format_str) eigenvector
      !normalization
      norm=0.0d0
      do i=1,gaussnum
            ! psi0
            norm=norm+k(i)**2*eigenvector(i)**2*s(i)
            ! psi2
            norm=norm+k(i)**2*eigenvector(i+gaussnum)**2*s(i)                  
      end do
      norm=dsqrt(norm)
      eigenvector=eigenvector/norm
      write(kwrite3,*) 'normalization N:'
      write(kwrite3,*) norm
      write(kwrite3,*) 'eigenvector in momentum phase:'
      write(kwrite3,format_str) eigenvector

      !the asymptotic_As
      gamma=dsqrt(-evr(indices(1))*masshat)/hbarc
      write(kwrite3,*) 'gamma:'
      write(kwrite3,*) gamma
      count=0
      avesum=0.0d0
      do i=1,50
            r=7.0d0+0.2d0*i
            u_r(i)=ufunction(r,k,s,eigenvector,gaussnum)
            asymptotic_u_r(i)=exp(-gamma*r)
            asymptotic_As=u_r(i)/asymptotic_u_r(i)
            if((r.gt.8.0d0) .and. (r.lt.14.0d0) )then
                  avesum=avesum+asymptotic_As
                  count=count+1
            else
                  continue
            end if
      end do
      asymptotic_As=avesum/count
      write(kwrite,*) 'As:'
      write(kwrite,*) asymptotic_As

      !the asymptotic_Ad
      count=0
      avesum=0.0d0
      do i=1,30
            r=10.0d0+0.2d0*i
            w_r(i)=wfunction(r,k,s,eigenvector(gaussnum+1),gaussnum)
            asymptotic_w_r2(i)=exp(-gamma*r)*(1.0d0+3.0d0/(gamma*r)+3.0d0/(gamma*r)**2.0d0)
            asymptotic_Ad=w_r(i)/asymptotic_w_r2(i)
            if((r.gt.11.0d0) .and. (r.lt.15.0d0) )then
                  avesum=avesum+asymptotic_Ad
                  count=count+1
            else
                  continue
            end if
      end do
      asymptotic_Ad=avesum/count
      write(kwrite,*) 'Ad:'
      write(kwrite,*) asymptotic_Ad
      !the ratio eta
      write(kwrite,*) 'eta:'
      eta=asymptotic_Ad/asymptotic_As
      write(kwrite,*) eta

      !the wave function in space phase
      kwrite1=11
      kwrite2=12
      kread=14
      open(unit=kwrite1,file='3S1_wave_function')
      open(unit=kwrite2,file='3D1_wave_function')
      write(kwrite1,*) '3S1 wave function:'
      write(kwrite2,*) '3D1 wave function:'
      open(unit=kread,file='wavefunction_r')
      read(kread,*,iostat=ios) r
      do while(ios .eq. 0)
            uf=ufunction(r,k,s,eigenvector,gaussnum)
            wf=wfunction(r,k,s,eigenvector(gaussnum+1),gaussnum)
            write(kwrite1,'(F6.2,F12.7)') r,uf
            write(kwrite2,'(F6.2,F12.7)') r,wf
            read(kread,*,iostat=ios) r
      end do
      
      ! the radius calculation
      !sr=s*C2p
      !kr=k*C2p
      allocate(kr(gaussnum))
      allocate(sr(gaussnum))

      radius_temp=0.0d0
      do i=1,gaussnum
            r=k(i)
            if(r.lt.15.0d0)then
                  uf=ufunction(r,k,s,eigenvector,gaussnum)
                  wf=wfunction(r,k,s,eigenvector(gaussnum+1),gaussnum)
            else
                  uf=asymptotic_As*exp(-gamma*r)
                  wf=asymptotic_Ad*exp(-gamma*r)*(1.0d0+3.0d0/(gamma*r)+3.0d0/(gamma*r)**2.0d0)
            end if
            radius_temp=radius_temp+s(i)*r**2*(uf**2+wf**2)
      end do
      radius_temp=dsqrt(radius_temp)
      radius=0.5d0*radius_temp
      write(kwrite,*) 'radius:'
      write(kwrite,*)  radius
      
      !the quadrupole moment
      quadrupolem=0.0d0
      do i=1,gaussnum
            r=k(i)
            if(r.lt.15.0d0)then
                  uf=ufunction(r,k,s,eigenvector,gaussnum)
                  wf=wfunction(r,k,s,eigenvector(gaussnum+1),gaussnum)
            else
                  uf=asymptotic_As*exp(-gamma*r)
                  wf=asymptotic_Ad*exp(-gamma*r)*(1.0d0+3.0d0/(gamma*r)+3.0d0/(gamma*r)**2.0d0)
            end if
            quadrupolem=quadrupolem+s(i)*r**2*wf*(dsqrt(8.0d0)*uf-wf)
      end do 
      quadrupolem=quadrupolem/20.0d0
      write(kwrite,*) 'quadrupole moment:'
      write(kwrite,*)  quadrupolem

      !the pd
      pd=0.0d0
      ps=0.0d0
      do i=1,gaussnum
            r=k(i)
            if(r.lt.15.0d0)then
                  uf=ufunction(r,k,s,eigenvector,gaussnum)
                  wf=wfunction(r,k,s,eigenvector(gaussnum+1),gaussnum)
            else
                  uf=asymptotic_As*exp(-gamma*r)
                  wf=asymptotic_Ad*exp(-gamma*r)*(1.0d0+3.0d0/(gamma*r)+3.0d0/(gamma*r)**2.0d0)
            end if
            pd=pd+s(i)*wf**2
            ps=ps+s(i)*uf**2
      end do
      write(kwrite,*) 'pd:'
      write(kwrite,*)  pd
      write(kwrite,*) 'norm test:'
      write(kwrite,*)  pd+ps

      deallocate(kr)
      deallocate(sr)
      deallocate(work)
      deallocate(A)
      deallocate(evr)
      deallocate(evi)
      deallocate(evector)
      deallocate(indices)
      deallocate(eigenvector)
!   k, s reset,they are not right

      inn=2  !np interaction
      j=0 ! 1S0 and 3D1

!set the matrix elements     
      dim=gaussnum 
      !create format_str
      write(format_str, "('(',I0,'F20.10)')") dim
      allocate(A(dim,dim))
      allocate(evector(dim,dim))
      allocate(evr(dim))
      allocate(evi(dim))
      A=0.0d0
      evi=0.0d0
      evr=0.0d0
      do row=1,gaussnum
            do column=1,gaussnum
                  momentumi=k(row)
                  momentumj=k(column)
                  xmev=momentumi
                  ymev=momentumj
                  call n3lo500new
                  !v(4)=v00(i,j) v(6)=v02(i,j) v(5)=v20(i,j) v(3)=22(i,j)
                  !delta_ij
                  if (row .eq. column)then
                        delta=1.0d0
                  else
                        delta=0.0d0
                  end if

                  !top-left block
                  A(row,column)=A(row,column)+v(1)*momentumj**2*s(column)+delta*momentumi**2/masshat

            end do
      end do
!write the matrix elements
      write(kwrite,*) '============================================================================='
      write(kwrite,*) '1S0'
!      write(kwrite,*) 'matrix elements:'
!      do row=1,dim
!            write(kwrite,format_str) (A(row,i),i=1,dim)       
!      end do

!solve the matrix's eigenvalues

      !initial dsyev
      lwork = -1
      allocate(work(1))
      call dgeev('N', 'V', dim, A, dim, evr, evi, A, dim, evector, dim, work, lwork, info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))

      !compute the eigenvalues and eigenvectors 
      call dgeev('N', 'V', dim, A, dim, evr, evi, A, dim, evector, dim, work, lwork, info)
      !write the results
      ! write(kwrite,*) 'eigenvalues-realpart:'
      ! write(kwrite,format_str) evr
      
      !Order from small to large
      allocate(indices(dim))
      call index_sort(evr, dim, indices)
      write(kwrite,*) 'eigenvalues-realpart-ordered:'
      write(kwrite,*) evr(indices(1))
      ! write(kwrite,*) 'eigenvalues-imaginarypart:'
      ! write(kwrite,format_str) evi
!      write(kwrite,*) 'eigenvectors:'
!      do row=1,dim
!            write(kwrite,format_str) (evector(row,i),i=1,dim)       
!      end do  
      deallocate(k)
      deallocate(s)
      deallocate(work)
      deallocate(A)
      deallocate(evr)
      deallocate(evi)
      close(kread)
      close(kwrite)
      close(kwrite1)
      close(kwrite2)
      close(kwrite3)


      end program deuteron


