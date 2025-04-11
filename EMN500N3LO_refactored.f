c***********************************************************************
c    version one :March 2022
c    version two :June 2022 (change the representation of theparameters
c    in lsj representation)
c    version three: 2023.3.29. LO
c    version four: 2023.4.6 NLO 
c    version five:2023.4.7 N2LO
c    version six:2023.4.10 N2LO rebulid codes
c    version seven:2023.4.14 N3LO except two loop terms
c    version eight:2023.4.18 N3LO pion terms finished
c    version nine:2023.4.20 N3LO all
c    version ten:2023.5.28 charge dependent
c    
c***********************************************************************  
c
c    author: shang haoyu 
c            physics school
c            Peking university
c    email:  shy@stu.pku.edu.cn
c
c
c***********************************************************************   
c    module kqxy contains the variable we use mostly in this code

c    function initialize evaluate the variables we use in this code

c    function normk, normq means the length of vector k and q
c    thier variable z means cos(theta)
      module potential_global
c     this module contains the global variable of the subroutine potential

c     variables:xmev,ymev,conta(24),lambda,j,v(6)
c     v(6) means
c     in the following order:
c     0v(singlet), 1v(uncoupled triplet), v++, v--, v+-, v-+ (coupled)
c     input
common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
c        arguments and values of this subroutine
c
      implicit real*8 (a-h,o-z)
      private

c     input
      public xmev,ymev,j,inn

c     output
      public v 
      
      logical heform,sing,trip,coup,endep
      character*4 label
      real*8 xmev,ymev
      integer j
      real*8 v(6)     
      common /cpot/   v,xmev,ymev
      common /cstate/ j,heform,sing,trip,coup,endep,label
      common /cnn/ inn
      integer :: inn=2
             
      end module


         module const

c   variable:pi,ga,mpi,mpi0,mpipm,mass,fpi,tidelambda,c1,c2,c3,c4
c   subroutine :ini_const
         public
         real*8 ::alpha 
         real*8 ::pi 
         real*8 ::ga 
         real*8 ::mpi 
         real*8 ::mpi0 
         real*8 ::mpipm 
         real*8 ::mnucleon(3)
         real*8 ::mass 
         real*8 ::fpi 
         real*8 ::tidelambda 
         real*8 ::lambda
         real*8 ::c1 
         real*8 ::c2 
         real*8 ::c3 
         real*8 ::c4 
         real*8 ::d15m14 
         real*8 ::d1p2 
         real*8 ::d3 
         real*8 ::d5 
         real*8 ::c1s0(3)
         real*8 ::c(24)
         contains
         subroutine ini_const
            use potential_global,only:inn
c      this subroutine should been used in the main programm            
            alpha=137.03599976d0
            pi=3.141592653589793d0
            ga=1.29d0
            mpi=138.0390d0
            mpi0=134.9766d0
            mpipm=139.5702d0
c       p mass           
            mnucleon(1)=938.2720d0
c       average mass
            mnucleon(2)=938.9182d0          
c       n mass
            mnucleon(3)=939.5654d0
            mass=mnucleon(inn)
            fpi=92.4d0
            tidelambda=650.0d0
            lambda=500.0d0
            c1=-1.07d0
            c2=3.20d0
            c3=-5.32d0
            c4=3.56d0
            d15m14=1.90d0
            d1p2=1.04d0 
            d3=-0.48d0 
            d5=0.14d0
c        p-p 1s0            
            c1s0(1)=-0.13753500d0
c        n-p 1s0            
            c1s0(2)=-0.13956293d0
c        n-n 1s0            
            c1s0(3)=-0.1385887d0
c         conta in p phase      
c         the sequence accoding to machleidt 2011
           c(1)=c1s0(inn)
         c(2)=2.417201616d0
         c(3)=-2.332142565d0
         c(4)=-16.73748278d0
         c(5)=1.18499088d0
         c(6)=4.989911864d0
         c(7)=0.187002997d0
         c(8)=9.976039933d0
         c(9)=-0.819122188d0
         c(10)=4.812595281d0
         c(11)=-0.159635365d0
         c(12)=0.823265695d0
         c(13)=-4.319199809d0
         c(14)=-19.17105287d0
         c(15)=-5.59034819d0
         c(16)=0.502604767d0
         c(17)=1.162495137d0
         c(18)=1.759338786d0
         c(19)=-1.946384037d0
         c(20)=-3.200942439d0
         c(21)=-0.757935632d0
         c(22)=6.02385494d0
         c(23)=0.010519022d0
         c(24)=-1.3336602d0
            mpi=mpi/mass
            mpi0=mpi0/mass
            mpipm=mpipm/mass
            fpi=fpi/mass
            tidelambda=tidelambda/mass
            c1=c1*1.0d-3*mass
            c2=c2*1.0d-3*mass
            c3=c3*1.0d-3*mass
            c4=c4*1.0d-3*mass
            d15m14=d15m14*1.0d-6*mass*mass
            d1p2=d1p2*1.0d-6*mass*mass
            d3=d3*1.0d-6*mass*mass
            d5=d5*1.0d-6*mass*mass
         end subroutine
      end module
     

      module paravari

c variable:  pi,ga,mpi,mpi0,mpipm,mass,fpi,tidelambda,c1,c2,c3,c4(const)
c             x,y,dwn,wnq,wn3,dwnq,x2,y2,c(24)
c subroutine: ini_paravari      
       
        use const
        real*8,save :: x,y,dwn,wnq,wn3,dwnq,x2,y2

        contains

          subroutine ini_paravari
            use potential_global,only:xmev,ymev
            real*8 t(24)
            real*8 matrix1(9,9),matrix2(15,15)
            logical :: parlsj=.true.          
            data((matrix1(j,i),i=1,9),j=1,9)/
     1 0.0198943d0,       0.0d0,       0.0d0,       0.0d0,        0.0d0,
     * 0.0596831d0,       0.0d0,       0.0d0,       0.0d0,
     2       0.0d0, 0.0099471d0,-0.0049735d0,-0.0149207d0, -0.0149207d0,
     *       0.0d0, 0.0298415d0,       0.0d0, -0.024867d0,
     3       0.0d0, 0.0397887d0, 0.0198943d0, 0.0596831d0,  0.0596831d0,
     *       0.0d0, 0.1193662d0,       0.0d0, 0.0994718d0,
     4 -0.019894d0,       0.0d0,       0.0d0,       0.0d0,        0.0d0,
     * 0.0198943d0,       0.0d0,       0.0d0,       0.0d0,
     5       0.0d0,-0.0099471d0,-0.0049735d0, 0.0149207d0,        0.0d0,
     *       0.0d0, 0.0099471d0, 0.0140674d0,-0.0099471d0,
     6       0.0d0, -0.039788d0, 0.0198943d0,-0.0596831d0,        0.0d0,
     *       0.0d0, 0.0397887d0, 0.0562697d0, 0.0397887d0,
     7       0.0d0,       0.0d0,-0.0397887d0,       0.0d0, -0.0596831d0,
     *       0.0d0,       0.0d0,       0.0d0, 0.0994718d0,
     8       0.0d0,       0.0d0, 0.0099471d0,       0.0d0, -0.0149207d0,
     *       0.0d0,       0.0d0,-0.0422023d0, 0.0049735d0,
     9       0.0d0,       0.0d0,-0.0397887d0,       0.0d0,  0.0596831d0,
     *       0.0d0,       0.0d0,-0.1688093d0,-0.0198943d0/
c   here we evaluate the number we use
            data((matrix2(j,i),i=1,15),j=1,15)/
     1  0.002486796d0,	0.001243398d0,	-0.002486796d0,	-0.007460388d0,
     *  -0.007460388d0,	0.007460388d0,	0.003730194d0,	0.003730194d0,	
     *  0.0d0,	0.0d0,	0.00621699d0,	0.00621699d0,	-0.01243398d0,	
     *  0.0d0,	0.008703786d0,	
     2  0.039788736d0,	0.019894368d0,	0.039788736d0,	0.119366207d0,	
     *  0.119366207d0,	0.119366207d0,	0.059683104d0,	0.059683104d0,	
     *  0.0d0,	0.0d0,	0.099471839d0,	0.099471839d0,	0.198943679d0,	
     *  0.0d0,	0.139260575d0,	
     3  0.059683104d0,	-0.009947184d0,	0.0d0,	0.0d0,	0.0d0,	
     *  0.179049311d0,	-0.029841552d0,	-0.029841552d0,	0.0d0,
     *  0.0d0,	-0.04973592d0,	-0.04973592d0,	0.0d0,	0.0d0,	
     * -0.069630288d0,	
     4 -0.039788736d0,	0.019894368d0,	0.0d0,	0.0d0,	0.0d0,	
     * -0.119366207d0,	0.059683104d0,	-0.029841552d0,	0.0d0,	
     * 0.0d0,	-0.04973592d0,	-0.04973592d0,	0.0d0,	0.0d0,
     *	-0.069630288d0,	
     5 -0.002486796d0,	-0.001243398d0,	-0.002486796d0,
     * 0.007460388d0,	0.0d0,	0.002486796d0,	0.001243398d0,
     * 0.002486796d0,	0.003516861d0,	0.003516861d0,	-0.00621699d0,
     * 0.0d0,	-0.004973592d0,	-0.006091381d0,	0.003730194d0,	
     6 -0.039788736d0,	-0.019894368d0,	0.039788736d0,	
     * -0.119366207d0,	0.0d0,	0.039788736d0,	0.019894368d0,
     * 0.039788736d0,	0.05626977d0,	0.05626977d0,	-0.099471839d0,
     *	0.0d0,	0.079577472d0,	0.0974621d0,	0.059683104d0,	
     7 -0.059683104d0,	0.009947184d0,	0.0d0,	0.0d0,	0.0d0,	
     * 0.059683104d0,	-0.009947184d0,	-0.019894368d0,
     * 0.084404655d0,	-0.028134885d0,	0.04973592d0,	0.0d0,
     * 0.0d0,	0.0d0,	-0.029841552d0,	
     8 0.039788736d0,	-0.019894368d0,	0.0d0,	0.0d0,	
     * 0.0d0,	-0.039788736d0,	0.019894368d0,	0.009947184d0,	
     * 0.028134885d0,	-0.028134885d0,	0.04973592d0,	-0.04973592d0,
     * 0.0d0,	0.0d0,	-0.009947184d0,	
     9 0.0d0,	0.0d0,	-0.019894368d0,	0.0d0,	-0.029841552d0,	
     * 0.0d0,	0.0d0,	0.044762328d0,	0.0d0,	0.0d0,	0.0d0,
     * 0.02486796d0,	0.04973592d0,	0.0d0,	-0.069630288d0,	
     * 0.0d0,	0.0d0,	-0.079577472d0,	0.0d0,	-0.119366207d0,
     * 0.0d0,	0.0d0,	-0.179049311d0,	0.0d0,	0.0d0,	0.0d0,
     *	-0.099471839d0,	0.198943679d0,	0.0d0,	0.27852115d0,	
     1 0.0d0,	0.0d0,	0.004973592d0,	0.0d0,	-0.007460388d0,	
     * 0.0d0,	0.0d0,	-0.003730194d0,	-0.010550582d0,	
     * -0.010550582d0,	0.0d0,	0.00621699d0,	0.002486796d0,
     * 0.018274144d0,	-0.002486796d0,	
     2 0.0d0,	0.0d0,	0.019894368d0,	0.0d0,	-0.029841552d0,	
     * 0.0d0,	0.0d0,	0.014920776d0,	-0.126606982d0,	0.042202327d0,
     * 0.0d0,	-0.02486796d0,	0.009947184d0,	-0.170558675d0,
     * 0.009947184d0,	
     3 0.0d0,	0.0d0,	-0.019894368d0,	0.0d0,	0.029841552d0,
     * 0.0d0,	0.0d0,	0.014920776d0,	-0.126606982d0,	0.042202327d0,
     * 0.0d0,	-0.02486796d0,	-0.009947184d0,	0.170558675d0,
     * 0.009947184d0,	
     4 0.0d0,	0.0d0,	-0.079577472d0,	0.0d0,	0.119366207d0,	
     * 0.0d0, 0.0d0,	-0.059683104d0,	-0.168809309d0,	
     * -0.168809309d0,0.0d0,	0.099471839d0,	-0.039788736d0,
     * -0.2923863d0,	-0.039788736d0,	
     5 0.0d0,	0.0d0,	0.0d0,	0.0d0,	0.0d0,	0.0d0,	0.0d0,
     *	-0.059683104d0,	-0.084404655d0,	0.084404655d0,	0.0d0,	
     * 0.099471839d0,	0.0d0,	0.0d0,	-0.039788736d0/
c this subroutine should be used in potential subroutine
cas it contains the variables often been used            
            dwn=1.0d0/mass
            wnq=mass*mass
            wn3=mass*mass*mass
            dwnq=dwn*dwn

            x=xmev*dwn
            y=ymev*dwn
            x2=x*x
            y2=y*y


         if(parlsj)then
            t(1)=c(1)
            t(2)=c(2)
            t(10)=c(3)
            t(11)=c(4)
            t(3)=c(5)
            t(12)=c(6)
            t(4)=c(7)
            t(13)=c(8)
            t(5)=c(9)
            t(14)=c(10)
            t(6)=c(11)
            t(7)=c(12)
            t(15)=c(13)
            t(16)=c(14)
            t(17)=c(15)
            t(8)=c(16)
            t(18)=c(17)
            t(19)=c(18)
            t(20)=c(19)
            t(21)=c(20)
            t(9)=c(21)
            t(22)=c(22)
            t(23)=c(23)
            t(24)=c(24)
            c=0.0d0
            do i=1,9
            do j=1,9
               c(i)=c(i)+matrix1(i,j)*t(j)
            end do 
            end do
            do i=10,24
            do j=10,24
               c(i)=c(i)+matrix2(i-9,j-9)*t(j)
            end do
            end do
         else
            t(1)=c(2)
            t(2)=c(3)
            t(3)=c(4)
            c(2)=t(2)
            c(3)=t(3)
            c(4)=t(1)            
         end if            
           c(1)=c(1)*0.01d0*wnq
           c(4)=c(4)*0.01d0*wnq
           c(2)=c(2)*wnq*wnq*1.0d-8
           c(3)=c(3)*wnq*wnq*1.0d-8
           c(5:9)=c(5:9)*wnq*wnq*1.0d-8
           c(10:24)=c(10:24)*wn3*wn3*1.0d-14
          end subroutine
         end module
         
         module genfunc

c variable: none
c function:normk,normq,qdotk,kcrossq2,wfunc,lfunc,afunc
            
            use paravari
            private pi,ga,mpi,mpi0,mpipm,mass,fpi,tidelambda,
     +       c1,c2,c3,c4,x,y,dwn,wnq,wn3,dwnq,x2,y2,c
c this module contains general functions will be used in
c potential subroutine
         contains
          real*8 function normk(z)
            implicit real*8 (a-h,o-z)
            real*8 z
            normk=dsqrt(x*x+y*y+2.0d0*x*y*z)/2.0d0
            return
          end function
        
          real*8 function normq(z)
            implicit real*8 (a-h,o-z)
            real*8 z    
            normq=dsqrt(x*x+y*y-2.0d0*x*y*z)
            return
          end function

          real*8 function qdotk()
            implicit real*8 (a-h,o-z)
            qdotk=(x*x-y*y)/2.0d0
            return
          end function

          real*8 function kcrossq2(z)
            implicit real*8 (a-h,o-z)
            real*8 z
            kcrossq2=normq(z)**2*normk(z)**2
     1      -(qdotk())**2
            return
          end function

c   function w(q) in (20) PRC 66,014002(2002)  
         
          real*8 function wfunc(z)
             implicit real*8 (a-h,o-z)
             real*8 z    
             wfunc=dsqrt(4.0d0*(mpi)**2+normq(z)**2)
             return
          end function
c   function L(q) in (19) PRC 66,014002(2002)         
          real*8 function lfunc(z)
            implicit real*8 (a-h,o-z)
            real*8 z  
            logical :: SFC=.true.
c   SFC means spectrum function cutoff            
            if (SFC)then
            lfunc=wfunc(z)/(2.0d0*normq(z))*dlog(((tidelambda)**2
     1       *(2.0d0*(mpi)**2+normq(z)**2)-2.0d0*(mpi)**2
     2      *normq(z)**2+tidelambda*dsqrt((tidelambda)**2
     3      -4.0d0*(mpi)**2)*normq(z)*wfunc(z))
     4      /(2.0d0*(mpi)**2*((tidelambda)**2+normq(z)**2)))
            else  
            lfunc=wfunc(z)/normq(z)*dlog((wfunc(z)+normq(z))
     1       /(2.0d0*mpi))
            end if 
            return
          end function
         
          real*8 function afunc(z)
          real*8 z
           logical :: SFC=.true.
           if (SFC)then 
            afunc=datan((normq(z)*(tidelambda-2.0d0*mpi))/
     1       (normq(z)**2+2.0d0*tidelambda*mpi))
     2      /(2.0d0*normq(z))
           else
c   tidelambda = infinity 
            afunc=datan(normq(z)/(2.0d0*mpi))/(2.0d0*normq(z))
           end if 
           return
           end function
        end module
  
      module velementm
         private
         public velement
         type velement
         real*8 ::vc(6)
         real*8 ::wc(6)
         real*8 ::vss(6)
         real*8 ::wss(6)
         real*8 ::vt(6)
         real*8 ::wt(6)
         real*8 ::vls(6)
         real*8 ::wls(6)
         real*8 ::vsk(6)
         real*8 ::vslsl(6)
         real*8 ::sum(6)
         contains 
         procedure :: init => ini_velement
         procedure :: add => add_velement
         end type
         contains
         subroutine ini_velement(this)
            class(velement) ::this
            this%vc=0.0d0
            this%vls=0.0d0
            this%vss=0.0d0
            this%vt=0.0d0
            this%wc=0.0d0
            this%wls=0.0d0
            this%wss=0.0d0
            this%wt=0.0d0
            this%vsk=0.0d0
            this%vslsl=0.0d0
            this%sum=0.0d0
            end subroutine
         subroutine add_velement(this)
            class(velement) ::this
            this%sum=this%vc+this%vls+this%vss+this%vt
     +    +this%wc+this%wls+this%wss+this%wt+this%vsk+this%vslsl
         end subroutine

      end module

      module lopot 
         use genfunc
         use paravari,only:c
         use velementm
         use const
         
         implicit none 
         private

c    the interface
c    variable
        public lo_ct,lo_onepi 
c    function
        public  onepii0,onepii1,onepi1
        public lo_vc_ct,lo_vs_ct      

        type(velement) ::lo_ct
        type(velement) ::lo_onepi
        contains
c    one-pi exchange term         
         real*8 function onepi(z,masspi)
          real*8 z,masspi
          onepi=-ga**2/(4.0d0*(fpi)**2
     +     *(normq(z)**2+(masspi)**2))
          return
          end function
c    the nn and pp onepi interaction
          real*8 function onepi1(z)
          real*8 z
          onepi1=onepi(z,mpi0)
          return
         end function
c    the np interaction         
c     charge dependent one-pi terms
c         I=0
          real*8 function onepii0(z)
          real*8 z
          onepii0=-onepi(z,mpi0)-2.0d0*onepi(z,mpipm)
          return
          end function
c         I=1
          real*8 function onepii1(z)
          real*8 z
          onepii1=-onepi(z,mpi0)+2.0d0*onepi(z,mpipm)
          return
          end function

c       contact terms
          real*8 function lo_vc_ct(z)
          real*8 z
          lo_vc_ct=c(1)
          return
          end function

          real*8 function lo_vs_ct(z)
          real*8 z
          lo_vs_ct=c(4)
          return
          end function          
      
      end module 

      module nlopot 
         use genfunc
         use paravari,only:c
         use velementm
         use const
         implicit none 
         private
         
c     the interface 
c     variable
        public nlo_ct,nlo_tp
c     function
        public nlowc,nlovss,nlovt
        public nlo_vc_ct,nlo_vs_ct,nlo_vls_ct,
     +  nlo_vt_ct,nlo_vsk_ct

        type(velement) :: nlo_ct
        type(velement) ::nlo_tp !two pi
        contains
c     two pion exchange terms        
          real*8 function nlowc(z)
          real*8 z
          nlowc=-lfunc(z)/(384.0d0*pi**2*(fpi)**4)
     +    *(4.0d0*(mpi)**2*(5.0d0*ga**4-4.0d0*ga**2-1.0d0)
     +    +normq(z)**2*(23.0d0*ga**4-10.0d0*ga**2-1.0d0)
     +    +48.0d0*ga**4*(mpi)**4/wfunc(z)**2)
          return
         end function

         real*8 function nlovss(z)
         real*8 z
         nlovss=-normq(z)**2*nlovt(z)
         return
        end function

        real*8 function nlovt(z)
        real*8 z
        nlovt=-3.0d0*ga**4/(64.0d0*pi**2*(fpi)**4)*lfunc(z)
        return
       end function

c     contact terms
          real*8 function nlo_vc_ct(z)
          real*8 z
          nlo_vc_ct=c(2)*normq(z)**2+c(3)*normk(z)**2
          return
          end function

          real*8 function nlo_vs_ct(z)
          real*8 z
          nlo_vs_ct=c(5)*normq(z)**2+c(6)*normk(z)**2 
          return
         end function

          real*8 function nlo_vls_ct(z)
          real*8 z
          nlo_vls_ct=c(7)
          return
          end function
         
          real*8 function nlo_vt_ct(z)
          real*8 z
          nlo_vt_ct=c(8)
          return
         end function

          real*8 function nlo_vsk_ct(z)
          real*8 z
          nlo_vsk_ct=c(9)
          return
          end function

      end module         
      
      module n2lopot
         use genfunc
         use paravari,only:c
         use velementm
         use const
         implicit none 
         private

c     the interface 
c     variable
         public n2lo_ol
c     function         
         public n2lo_vc_ol,n2lo_ws_ol,n2lo_wt_ol
         type(velement) ::n2lo_ol !ol means two pion one loop   
        contains
c     one loop two pion exchange         
        real*8 function n2lo_vc_ol(z)
          real*8 z
          n2lo_vc_ol=3.0d0*ga**2/(16.0d0*pi*(fpi)**4)
     1     *(2.0d0*(mpi)**2*(c3-2.0d0*c1)
     2     +c3*normq(z)**2)*(2.0d0*(mpi)**2
     3     +normq(z)**2)*afunc(z)
          return
        end function

          real*8 function n2lo_ws_ol(z)
          real*8 z
          n2lo_ws_ol=-normq(z)**2*n2lo_wt_ol(z)
          return
          end function

          real*8 function n2lo_wt_ol(z)
          real*8 z
          n2lo_wt_ol=-ga**2/(32.0d0*pi*(fpi)**4)
     1   *c4*wfunc(z)**2*afunc(z)
          return
          end function          


                      
      end module
      module n3lopot
c to improve the speed of calculation,you should use 
c n3lo_ini before using subroutine n3lo500new(in the main programm)                 
         use velementm
         use const
         use genfunc
         use paravari,only:c

         implicit none 
         private

c     the interface
c     variable
         public n3lo_fd,n3lo_rc,n3lo_tl,n3lo_cM,n3lo_ct
c     function
         public n3lo_vc_fd,n3lo_wt_fd,n3lo_ws_fd 
         public n3lo_vc_rc,n3lo_wc_rc,n3lo_vt_rc,n3lo_wt_rc,
     +    n3lo_vs_rc,n3lo_ws_rc,n3lo_vls_rc,n3lo_wls_rc
         public n3lo_vc_tl,n3lo_ws_tl,n3lo_wt_tl,n3lo_vs_tl,n3lo_vt_tl,
     +    n3lo_wc_tl
         public n3lo_vc_cM,n3lo_wc_cM,n3lo_wt_cM,n3lo_ws_cM,
     +    n3lo_vls_cM,n3lo_wls_cM
         public n3lo_vc_ct,n3lo_vt_ct,n3lo_vls_ct,n3lo_vs_ct,
     +    n3lo_vsk_ct,n3lo_vslsl_ct
         public  pi_gamma  
c     subroutine         
         public n3lo_ini

         type(velement) ::n3lo_fd
         type(velement) ::n3lo_rc
         type(velement) ::n3lo_tl
         type(velement) ::n3lo_cM
         type(velement) ::n3lo_ct
         real*8 :: imvs(96),imwc(96),imvc(96),
     +   imws(96),imwt(96),imvt(96)    

         contains
         subroutine n3lo_ini
            implicit none
            real*8 wt(96),ct(96)
            real*8 xlb
            integer i
            xlb=2.0d0*mpi
            call fset(ct,wt,xlb,tidelambda,96)
            do i=1,96
            imvs(i)=n3lo_imvs(ct(i))
            imwc(i)=n3lo_imwc(ct(i))
            imvc(i)=n3lo_imvc(ct(i))
            imws(i)=n3lo_imws(ct(i))
            imwt(i)=n3lo_imwt(ct(i))
            imvt(i)=n3lo_imvt(ct(i))            
            end do
         end subroutine
c     football diagram(the 'fd')
         real*8 function n3lo_vc_fd(z)
         real*8 z
      n3lo_vc_fd=3.0d0/(16.0d0*pi**2*fpi**4)*((c2/6.0d0*wfunc(z)**2
     +   +c3*(2.0d0*mpi**2+normq(z)**2)-4.0d0*c1*mpi**2)**2+c2**2
     +   /45.0d0*wfunc(z)**4)*lfunc(z)
        return
        end function

        real*8 function n3lo_wt_fd(z)
        real*8 z
        n3lo_wt_fd=c4**2/(96.0d0*pi**2*fpi**4)*wfunc(z)**2*lfunc(z)
        return
      end function

        real*8 function n3lo_ws_fd(z)
        real*8 z
        n3lo_ws_fd=-normq(z)**2*n3lo_wt_fd(z)
        return
        end function
c     relativistic corrections(the 'rc')
        
        real*8 function n3lo_vc_rc(z)
        real*8 z
        n3lo_vc_rc=3.0d0*ga**4/(128.0d0*pi*fpi**4)
     +  *(mpi**5/(2.0d0*wfunc(z)**2)+(2.0d0*mpi**2
     +  +normq(z)**2)*(normq(z)**2-mpi**2)*afunc(z))
        return
      end function
      
        real*8 function n3lo_wc_rc(z)
        real*8 z
        n3lo_wc_rc=ga**2/(64.0d0*pi*fpi**4)*(3*ga**2
     +  *mpi**5/(2.0d0*wfunc(z)**2)+(ga**2*(3.0d0*mpi**2
     + +2.0d0*normq(z)**2)-2.0d0*mpi**2-normq(z)**2)
     +  *(2.0d0*mpi**2+normq(z)**2)*afunc(z))
        return
       end function
      
        real*8 function n3lo_vt_rc(z)
        real*8 z
        n3lo_vt_rc=3.0d0*ga**4/(256.0d0*pi*fpi**4)
     +   *(5.0d0*mpi**2+2.0d0*normq(z)**2)*afunc(z)
        return
      end function

        real*8 function n3lo_wt_rc(z)
        real*8 z
        n3lo_wt_rc=ga**2/(128.0d0*pi*fpi**4)
     +  *(ga**2*(3.0d0*mpi**2+normq(z)**2)-wfunc(z)**2)
     +  *afunc(z)
        return 
        end function

        real*8 function n3lo_vs_rc(z)
        real*8 z
        n3lo_vs_rc=-normq(z)**2*n3lo_vt_rc(z)
        return
      end function

        real*8 function n3lo_ws_rc(z)
        real*8 z
        n3lo_ws_rc=-normq(z)**2*n3lo_wt_rc(z)
        return
        end function

        real*8 function n3lo_vls_rc(z)
        real*8 z
        n3lo_vls_rc=3.0d0*ga**4/(32.0d0*pi*fpi**4)
     +  *(2.0d0*mpi**2+normq(z)**2)*afunc(z)
        return
      end function

        real*8 function n3lo_wls_rc(z)
        real*8 z
        n3lo_wls_rc=ga**2*(1.0d0-ga**2)/(32.0d0*pi*
     +   fpi**4)*wfunc(z)**2*afunc(z)
        return
       end function

c     spectral functions,name as n3lo_imv(mu)
      
        real*8 function n3lo_imvc(mu)
        real*8 mu
        n3lo_imvc=3*ga**4*(2.0d0*mpi**2-mu**2)/(pi*mu
     +  *(4.0d0*fpi)**6)*((mpi**2-2.0d0*mu**2)*(2.0d0*mpi
     +  +(2.0d0*mpi**2-mu**2)/(2.0d0*mu)*dlog((mu+2.0d0*mpi)
     + /(mu-2.0d0*mpi)))+4.0d0*ga**2*mpi*(2.0d0*mpi**2-mu**2))
        return
      end function

        real*8 function n3lo_imws(mu)
        real*8 mu
        n3lo_imws=ga**4*(4*mpi**2-mu**2)/(pi*(4.0d0*fpi)**6)
     +   *((mpi**2-mu**2/4.0d0)*dlog((mu+2.0d0*mpi)/(mu-2.0d0*mpi))
     +  +(1.0d0+2.0d0*ga**2)*mu*mpi)
        return
        end function

        real*8 function n3lo_imwt(mu)
        real*8 mu
        n3lo_imwt=n3lo_imws(mu)/mu**2
        return
      end function

        real*8 function n3lo_imvs(mu)
        implicit none
        real*8 mu
        real*8 rk  !rootk
        real*8 integral,xr  !x regulation
        real*8 wt(96),ct(96)
        integer i
        call fset(ct,wt,0.0d0,1.0d0,12)
        integral=0.0d0
        rk=dsqrt(mu**2/4.0d0-mpi**2)
        do i=1,12
        xr=ct(i)
        integral=integral+wt(i)*((1.0d0-xr**2)*(1.0d0/6.0d0
     +  -mpi**2/(rk**2*xr**2)+(1.0d0+mpi**2/(rk**2*xr**2))**(1.5d0)
     +   *dlog((rk*xr+dsqrt(mpi**2+rk**2*xr**2))/mpi)))  
        end do
        n3lo_imvs=ga**2*mu*rk**3/(8.0d0*pi*fpi**4)*d15m14
     +  +2.0d0*ga**6*mu*rk**3/(8.0d0*pi*fpi**2)**3*integral
        return
        end function

        real*8 function n3lo_imvt(mu)
        real*8 mu
        n3lo_imvt=n3lo_imvs(mu)/mu**2
        return
      end function
        
        real*8 function n3lo_imwc(mu)
        implicit none
        real*8 mu
        real*8 rk,wt(96),ct(96),xr
        real*8 term1,term21,term22,term23,term24,term
        real*8 integral
        integer ::i
        integral=0.0d0
        call fset(ct,wt,0.0d0,1.0d0,12)
        rk=dsqrt(mu**2/4.0d0-mpi**2)
        do i=1,12
        xr=ct(i)
        term1=ga**2*(mu**2-2.0d0*mpi**2)+2.0d0*
     +   (1.0d0-ga**2)*rk**2*xr**2
        term21=96.0d0*pi**2*fpi**2*((2.0d0*mpi**2-mu**2)
     +   *d1p2-2.0d0*rk**2*xr**2*d3+4.0d0*mpi**2*d5)
        term22=(4.0d0*mpi**2*(1.0d0+2.0d0*ga**2)-mu**2
     +   *(1.0d0+5.0d0*ga**2))*rk/mu*dlog((mu+2.0d0*rk)
     +  /(2.0d0*mpi))+mu**2/12.0d0*(5.0d0+13.0d0*ga**2)
     +  -2.0d0*mpi**2*(1.0d0+2.0*ga**2)
        term23=-3.0d0*rk**2*xr**2+6.0d0*rk*xr*dsqrt(mpi**2
     +   +rk**2*xr**2)*dlog((rk*xr+dsqrt(mpi**2+rk**2*xr**2))/mpi)
        term24=ga**4*(mu**2-2.0d0*rk**2*xr**2-2.0d0*mpi**2)
     +   *(5.0d0/6.0d0+mpi**2/(rk**2*xr**2)-(1.0d0+mpi**2
     +   /(rk**2*xr**2))**(1.5d0)*dlog((rk*xr+dsqrt(mpi**2
     +   +rk**2*xr**2))/mpi))
        term=term1*(term21+term22+term23+term24)
        integral=integral+wt(i)*term
        end do
        n3lo_imwc=(2.0d0*rk)/(3.0d0*mu*(8.0d0*pi*fpi**2)**3)*integral
        return
        end function
         
c     two loop contributions(tl)

c     v_{C,S}=-2q^2/pi*integrate_{2mpi,Lambda}(imvcs(i mu)/(mu^5(mu^2+q^2)))
        real*8 function vcstl(func,z)
c      input
        real*8 z
        real*8 func(96)
c     local 
        real*8 wt(96),ct(96)
        real*8 xlb
        real*8 integral
        integer i
        xlb=2.0d0*mpi
        call fset(ct,wt,xlb,tidelambda,96)
c     the integral
        integral=0.0d0
        do i=1,96
         integral=integral+wt(i)*(func(i)/
     +    (ct(i)**5*(ct(i)**2+normq(z)**2)))
        end do
        vcstl=-2.0d0*normq(z)**6/pi*integral
        return
      end function

c     v_{T,LS}
        real*8 function vtlstl(func,z)
c      input
        real*8 z
        real*8 func(96)
c     local 
        real*8 wt(96),ct(96)
        real*8 xlb
        real*8 integral
        integer i
        xlb=2.0d0*mpi
        call fset(ct,wt,xlb,tidelambda,96)
c     the integral
        integral=0.0d0
        do i=1,96
         integral=integral+wt(i)*(func(i)/
     +    (ct(i)**3*(ct(i)**2+normq(z)**2)))
        end do
        vtlstl=2.0d0*normq(z)**4/pi*integral
        return
      end function

      real*8 function n3lo_vc_tl(z)
      real*8 z
      n3lo_vc_tl=vcstl(imvc,z)
      return
      end function
      
      real*8 function n3lo_ws_tl(z)
      real*8 z
      n3lo_ws_tl=vcstl(imws,z)
      return
      end function 

      real*8 function n3lo_wt_tl(z)
      real*8 z
      n3lo_wt_tl=vtlstl(imwt,z)
      return
      end function

      real*8 function n3lo_vs_tl(z)
      real*8 z
      n3lo_vs_tl=vcstl(imvs,z)
      return
      end function 

      real*8 function n3lo_vt_tl(z)
      real*8 z
      n3lo_vt_tl=vtlstl(imvt,z)
      return
      end function

      real*8 function n3lo_wc_tl(z)
      real*8 z
      n3lo_wc_tl=vcstl(imwc,z)
      return
      end function

c    ci/M contributions('cM')
      real*8 function n3lo_vc_cM(z)
      real*8 z
      n3lo_vc_cM=-ga**2*lfunc(z)/(32.0d0*pi**2*fpi**4)
     +  *((c2-6.0d0*c3)*normq(z)**4+4.0d0*(6.0d0*c1+c2-3.0d0*c3)
     +  *normq(z)**2*mpi**2+6.0d0*(c2-2.0d0*c3)*mpi**4
     +  +24.0d0*(2.0d0*c1+c3)*mpi**6/wfunc(z)**2)
      return
      end function

      real*8 function n3lo_wc_cM(z)
      real*8 z
      n3lo_wc_cM=-c4*normq(z)**2*lfunc(z)/(192.0d0*pi**2*fpi**4)
     + *(ga**2*(8.0d0*mpi**2+5.0d0*normq(z)**2)+wfunc(z)**2)
      return
      end function
      
      real*8 function n3lo_wt_cM(z)
      real*8 z
      n3lo_wt_cM=-c4*lfunc(z)/(192.0d0*pi**2*fpi**4)*
     + (ga**2*(16.0d0*mpi**2+7.0d0*normq(z)**2)-wfunc(z)**2)
      return
      end function

      real*8 function n3lo_ws_cM(z)
      real*8 z
      n3lo_ws_cM=-n3lo_wt_cM(z)*normq(z)**2
      return
      end function

      real*8 function n3lo_vls_cM(z)
      real*8 z
      n3lo_vls_cM=c2*ga**2/(8.0d0*pi**2*fpi**4)
     + *wfunc(z)**2*lfunc(z)
      return
      end function

      real*8 function n3lo_wls_cM(z)
      real*8 z
      n3lo_wls_cM=-c4*lfunc(z)/(48.0d0*pi**2*fpi**4)*
     + (ga**2*(8.0d0*mpi**2+5.0d0*normq(z)**2)+wfunc(z)**2)
      return
      end function

c   contact terms(ct)
      real*8 function n3lo_vc_ct(z)
      real*8 z
       n3lo_vc_ct=c(10)*normq(z)**4+c(11)*normk(z)**4
     + +c(12)*normq(z)**2*normk(z)**2+c(13)*kcrossq2(z)
       return
      end function

      real*8 function n3lo_vs_ct(z)
      real*8 z
      n3lo_vs_ct=c(14)*normq(z)**4+c(15)*normk(z)**4
     + +c(16)*normq(z)**2*normk(z)**2+c(17)*kcrossq2(z)
       return
       end function

       real*8 function n3lo_vls_ct(z)
       real*8 z
       n3lo_vls_ct=c(18)*normq(z)**2+c(19)*normk(z)**2 
       return
      end function

       real*8 function n3lo_vt_ct(z)
       real*8 z
       n3lo_vt_ct=c(20)*normq(z)**2+c(21)*normk(z)**2
       return
       end function

       real*8 function n3lo_vsk_ct(z)
       real*8 z
       n3lo_vsk_ct=c(22)*normq(z)**2+c(23)*normk(z)**2
       return
      end function

       real*8 function n3lo_vslsl_ct(z)
       real*8 z
       n3lo_vslsl_ct=c(24)
       return
       end function


c    pi-gamma (charge dependent) (tensor)
        real*8 function pi_gamma(z)
        real*8 z
        real*8 beta
        beta=normq(z)/mpi
        pi_gamma=-ga**2/(2.0d0*fpi**2*mpi**2*pi*alpha)
     +  *(-(1.0d0-beta**2)**2/(2.0d0*beta**4
     +  *(1.0d0+beta**2))*dlog(1.0d0+beta**2)+1.0d0/
     +  (2.0d0*beta**2))
        return
      end function
      end module

      


c   we write another legendre 
        function legendre(x,j)
c
c
c   x is the independent variable
c
        real*8  x,legendrem1,a,b
        real*8  legendre
        integer j
c
c
c
c        compute legendre polynom for j equals zero
c
c
        if (j.gt.0) go to 1
        legendre=1.d0
        legendrem1=0.d0
        if (j.lt.0) legendre=0.d0
        return
c
c
c
c        compute legendre polynoms for j equals one
c
c
c
    1 legendre=x
        legendrem1=1.d0
        if (j.eq.1) return
c
c
c
c        compute legendre polynom for j greater or equal two
c
c
c
        do i=2,j
        a=x*legendre
        b=a-legendrem1
        legendrem1=legendre
        legendre=-b/dfloat(i)+b+a
        end do
c

        return
        end

        function alj(func,l,j)
c   function alj calculate the integral a^j(l), l is obit-angular
c   momentum and j is total angular momentum,func is the function 
c   being integrated

c   ALJ^J(q',q)=pi*INTEGRATE[fuc(q',q,z) z^l P_J(z)，-1,1]
        implicit none
c    input        
        real*8,external :: func
        real*8,external :: legendre
        integer :: l,j
c    output        
        real*8 :: alj
c    local
        real*8 :: pi   
        real*8 ::ct(96),wt(96)  
        integer ::i  

        alj=0.0d0
        pi=3.141592653589793d0
        call fset(ct,wt,-1.0d0,1.0d0,24)
        do i=1,24
        alj=alj+wt(i)*(legendre(ct(i),j)*func(ct(i))*ct(i)**l*pi)
        end do
        return
        end

c   this function calculate the cutoff
c   lambda is cutoff energy ,n adjust the sharp degree 

        function cutoff(lambda,n)
         use potential_global,only:xmev,ymev
        implicit real*8 (a-h,o-z)
        real*8 lambda,t,expo
        real*8 cutoff
        integer n
        t=dfloat(n)
        expo=(xmev/lambda)**(2.0d0*t)+(ymev/lambda)**(2.0d0*t)
        cutoff=dexp(-expo)  
        end
        module decompose
c   this module calculate lsj decomposition
c   we use pot(6) to record them
c   0v(singlet), 1v(uncoupled triplet), v++, v--, v+-, v-+ (coupled) 
         implicit none
c   all the subroutines are public                 

         contains
      subroutine lsjvcentral(pot,vcentral,j)

c        input         
         real*8,external :: vcentral
         integer j
c        output
         real*8 pot(6)
c        local
         real*8,external :: alj
         
         pot=0.0d0
         if (j .eq. 0) then
         pot(1)=2.0d0*alj(vcentral,0,j)
         pot(3)=2.0d0*alj(vcentral,0,j+1)
         else 
         pot(1)=2.0d0*alj(vcentral,0,j)
         pot(2)=2.0d0*alj(vcentral,0,j)
         pot(3)=2.0d0*alj(vcentral,0,j+1)
         pot(4)=2.0d0*alj(vcentral,0,j-1)
         pot(5)=0.0d0
         pot(6)=0.0d0
         end if
        end   

        subroutine lsjvspinspin(pot,vspinspin,j)

c        input         
         real*8,external :: vspinspin
         integer j
c        output
         real*8 pot(6)
c        local
         real*8,external :: alj

         pot=0.0d0        
         if (j .eq. 0)then
            pot(1)=-6.0d0*alj(vspinspin,0,j)
            pot(3)=2.0d0*alj(vspinspin,0,j+1)
            else 
            pot(1)=-6.0d0*alj(vspinspin,0,j)
            pot(2)=2.0d0*alj(vspinspin,0,j)
            pot(3)=2.0d0*alj(vspinspin,0,j+1)
            pot(4)=2.0d0*alj(vspinspin,0,j-1)
            pot(5)=0.0d0
            pot(6)=0.0d0
            end if   
        end   
       
      subroutine lsjvtensor(pot,vtensor,j)

c        global(x,y)
         use paravari,only:x,y,x2,y2
                
c        input         
         real*8,external :: vtensor
         integer j
c        output
         real*8 pot(6)
c        local
         real*8,external :: alj
         real*8 jd,jdp1,jd2p1
         
         pot=0.0d0
         jd=dfloat(j)
         jdp1=jd+1.0d0
         jd2p1=2.0d0*jd+1.0d0
         if (j .eq. 0) then
        pot(1)=2.0d0*(-(x2+y2)*alj(vtensor,0,j)
     1  +2.0d0*x*y*alj(vtensor,1,j))
        pot(3)=2.0d0/jd2p1*(-(x2+y2)*alj(vtensor,0,j+1)
     1  +2.0d0*x*y*alj(vtensor,0,j))
        else
        pot(1)=2.0d0*(-(x2+y2)*alj(vtensor,0,j)
     1  +2.0d0*x*y*alj(vtensor,1,j))
        pot(2)=2.0d0*((x2+y2)*alj(vtensor,0,j)
     1  -2.0d0*x*y/jd2p1*(jd*alj(vtensor,0,j+1)
     2  +jdp1*alj(vtensor,0,j-1)))
        pot(3)=2.0d0/jd2p1*(-(x2+y2)*alj(vtensor,0,j+1)
     1  +2.0d0*x*y*alj(vtensor,0,j))
        pot(4)=2.0d0/jd2p1*((x2+y2)*alj(vtensor,0,j-1)
     1  -2.0d0*x*y*alj(vtensor,0,j))
        pot(5)=-4.0d0*dsqrt(jd*jdp1)/jd2p1*(y2*alj(vtensor,0,j+1)
     1  +x2*alj(vtensor,0,j-1)-2.0d0*x*y*alj(vtensor,0,j) )   
        pot(6)=-4.0d0*dsqrt(jd*jdp1)/jd2p1*(y2*alj(vtensor,0,j-1)
     1  +x2*alj(vtensor,0,j+1)-2.0d0*x*y*alj(vtensor,0,j) )
        end if 
        end 

        subroutine lsjvspinobit(pot,vspinobit,j)
c        global(x,y)
         use paravari,only:x,y
                
c        input         
         real*8,external :: vspinobit
         integer j
c        output
         real*8 pot(6)
c        local
         real*8,external :: alj
         real*8 jd,jdp1,jd2p1
         
         pot=0.0d0
         jd=dfloat(j)
         jdp1=jd+1.0d0
         jd2p1=2.0d0*jd+1.0d0
        if (j .eq. 0) then
        pot(1)=pot(1)+0.0d0
        pot(3)=pot(3)+2.0d0*x*y*(jd+2.0d0)/(2.0d0*jd+3.0d0)
     1  *(alj(vspinobit,0,j+2)-alj(vspinobit,0,j))
        else
        pot(1)=pot(1)+0.0d0
        pot(2)=pot(2)+2.0d0*x*y/jd2p1*(alj(vspinobit,0,j+1)
     1  -alj(vspinobit,0,j-1))
        pot(3)=pot(3)+2.0d0*x*y*(jd+2.0d0)/(2.0d0*jd+3.0d0)
     1  *(alj(vspinobit,0,j+2)-alj(vspinobit,0,j))
        pot(4)=pot(4)+2.0d0*x*y*(jd-1.0d0)/(2.0d0*jd-1.0d0)
     1  *(alj(vspinobit,0,j-2)-alj(vspinobit,0,j))
        pot(5)=pot(5)+0.0d0
        pot(6)=pot(6)+0.0d0
        end if
       end subroutine

       subroutine lsjvsigmal(pot,vsigmaL,j)
c        global(x,y,x2,y2)
         use paravari,only:x2,y2
                
c        input         
         real*8,external :: vsigmaL
         integer j
c        output
         real*8 pot(6)
c        local
         real*8,external :: alj
         real*8 jd,jdp1,jd2p1
         
         pot=0.0d0
         jd=dfloat(j)
         jdp1=jd+1.0d0
         jd2p1=2.0d0*jd+1.0d0
         if (j .eq. 0) then
         pot(1)=pot(1)+2.0d0*x2*y2*(alj(vsigmaL,2,j)-alj(vsigmaL,0,j))
         pot(3)=pot(3)+2.0d0*x2*y2*((2.0d0*jd+3.0d0)/jd2p1
     1  *alj(vsigmaL,0,j+1)-2.0d0/jd2p1*alj(vsigmaL,1,j)
     2  -alj(vsigmaL,2,j+1))
         else
         pot(1)=pot(1)+2.0d0*x2*y2*(alj(vsigmaL,2,j)-alj(vsigmaL,0,j))
         pot(2)=pot(2)+2.0d0*x2*y2*(-alj(vsigmaL,0,j)+(jd-1.0d0)/jd2p1
     1  *alj(vsigmaL,1,j+1)+(jd+2.0d0)/jd2p1*alj(vsigmaL,1,j-1))
         pot(3)=pot(3)+2.0d0*x2*y2*((2.0d0*jd+3.0d0)/jd2p1
     1  *alj(vsigmaL,0,j+1)-2.0d0/jd2p1*alj(vsigmaL,1,j)
     2  -alj(vsigmaL,2,j+1))
         pot(4)=pot(4)+2.0d0*x2*y2*((2*jd-1.0d0)/jd2p1
     1  *alj(vsigmaL,0,j-1)+2.0d0/jd2p1*alj(vsigmaL,1,j)
     2  -alj(vsigmaL,2,j-1))
         pot(5)=pot(5)-4.0d0*x2*y2*dsqrt(jd*jdp1)/(jd2p1)**2*
     1  (alj(vsigmaL,0,j+1)-alj(vsigmaL,0,j-1))
         pot(6)=pot(5)
         end if 
         return
         end subroutine

         subroutine lsjvsigmak(pot,vsigmak,j)
c        global(x,y,x2,y2)
         use paravari,only:x,y,x2,y2
                
c        input         
         real*8,external :: vsigmak
         integer j
c        output
         real*8 pot(6)
c        local
         real*8,external :: alj
         real*8 jd,jdp1,jd2p1
         
         pot=0.0d0
         jd=dfloat(j)
         jdp1=jd+1.0d0
         jd2p1=2.0d0*jd+1.0d0
         if (j .eq. 0) then
        pot(1)=pot(1)+0.5d0*(-(x2+y2)*alj(vsigmak,0,j)
     1  -2.0d0*x*y*alj(vsigmak,1,j))
        pot(3)=pot(3)+0.5d0/jd2p1*(-(x2+y2)*alj(vsigmak,0,j+1)
     1  -2.0d0*x*y*alj(vsigmak,0,j))
        else
        pot(1)=pot(1)+0.5d0*(-(x2+y2)*alj(vsigmak,0,j)
     1  -2.0d0*x*y*alj(vsigmak,1,j))
        pot(2)=pot(2)+0.5d0*((x2+y2)*alj(vsigmak,0,j)
     1  +2.0d0*x*y/jd2p1*(jd*alj(vsigmak,0,j+1)
     2  +jdp1*alj(vsigmak,0,j-1)))
        pot(3)=pot(3)+0.5d0/jd2p1*(-(x2+y2)*alj(vsigmak,0,j+1)
     1  -2.0d0*x*y*alj(vsigmak,0,j))
        pot(4)=pot(4)+0.5d0/jd2p1*((x2+y2)*alj(vsigmak,0,j-1)
     1  +2.0d0*x*y*alj(vsigmak,0,j))
        pot(5)=pot(5)-1.0d0*dsqrt(jd*jdp1)/jd2p1*(y2*alj(vsigmak,0,j+1)
     1  +x2*alj(vsigmak,0,j-1)+2.0d0*x*y*alj(vsigmak,0,j) )   
        pot(6)=pot(6)-1.0d0*dsqrt(jd*jdp1)/jd2p1*(y2*alj(vsigmak,0,j-1)
     1  +x2*alj(vsigmak,0,j+1)+2.0d0*x*y*alj(vsigmak,0,j) )
        end if
        return
       end subroutine
    
        subroutine isospindependent(pot1,j,pot2)
c       this subroutine contains the tau1 dot tau2 terms factor

c       input 
         real*8 pot1(6)
         integer j

c       output
         real*8 pot2(6)
         pot2=0.0d0
         if (mod(j,2).eq.0)then
            pot2(1)=pot1(1)
            pot2(2)=-3.0d0*pot1(2)
            pot2(3:6)=pot1(3:6)
         else
            pot2(1)=-3.0d0*pot1(1)
            pot2(2)=pot1(2)
            pot2(3:6)=-3.0d0*pot1(3:6)
         end if                
         end
      end module
        subroutine lsjdecomposition(pot,j)
        use paravari
        use lopot
        use nlopot
        use n2lopot
        use n3lopot
        use decompose
        use potential_global,only:inn
        implicit none
        real*8 pot(6),jd,jdp1,jd2p1,temp1(6),temp2(6),pigamma(6)
        real*8,external :: cutoff
        real*8,external :: alj
        real*8 ::ex,ey,ree,expexp2,expexp3,expexp4
        integer ::j,i
        call ini_const
        call ini_paravari
        call n3lo_ini
        pot=0.0d0
        jd=dfloat(j)
        jdp1=jd+1.0d0
        jd2p1=2.0d0*jd+1.0d0

c   central force part
c   when j=0,there is only two terms do not equal to zero
            
        ex=dsqrt(1.0d0+x*x)
        ey=dsqrt(1.0d0+y*y)
        ree=dsqrt(ex*ey)
        expexp2=cutoff(lambda,2)
        expexp3=cutoff(lambda,3)
        expexp4=cutoff(lambda,4)
c       lo 

c       contact terms
        call lo_ct%init
        call lo_onepi%init
        call lsjvcentral(lo_ct%vc,lo_vc_ct,j)
        call lsjvspinspin(lo_ct%vss,lo_vs_ct,j)
        call lo_ct%add
        lo_ct%sum=lo_ct%sum*expexp3

c       one-pion terms
        select case(inn)
        case(1)
c       the pp onepi exchange
        call lsjvtensor(lo_onepi%vt,onepi1,j)        
        case(2)
c       the np onepi exchange         
        call lsjvtensor(temp1,onepii0,j)
        call lsjvtensor(temp2,onepii1,j)
        if(mod(j,2) .eq. 1)then 
         lo_onepi%vt(1)=temp1(1)
         lo_onepi%vt(2)=temp2(2)
         lo_onepi%vt(3:6)=temp1(3:6)
         else
         lo_onepi%vt(1)=temp2(1)
         lo_onepi%vt(2)=temp1(2)
         lo_onepi%vt(3:6)=temp2(3:6)
         end if
         case(3)
c        the nn onepi exchange
         call lsjvtensor(lo_onepi%vt,onepi1,j)
         end select            
         call lo_onepi%add
         lo_onepi%sum=lo_onepi%sum*expexp4
c        nlo

c        contact terms
         call nlo_ct%init
         call lsjvcentral(nlo_ct%vc,nlo_vc_ct,j)
         call lsjvspinspin(nlo_ct%vss,nlo_vs_ct,j)
         call lsjvspinobit(nlo_ct%vls,nlo_vls_ct,j)
         call lsjvtensor(nlo_ct%vt,nlo_vt_ct,j)
         call lsjvsigmak(nlo_ct%vsk,nlo_vsk_ct,j)
         call nlo_ct%add
c        the C_{3P1} term n=3         
         if(j .eq. 1)then
            nlo_ct%sum(2)=nlo_ct%sum(2)*expexp3
            nlo_ct%sum(1)=nlo_ct%sum(1)*expexp2
            nlo_ct%sum(3:6)=nlo_ct%sum(3:6)*expexp2
         else
            nlo_ct%sum=nlo_ct%sum*expexp2
         end if

c       2-pi terms
        call nlo_tp%init 
        call lsjvcentral(temp1,nlowc,j)
        call isospindependent(temp1,j,nlo_tp%wc)
        call lsjvspinspin(nlo_tp%vss,nlovss,j)
        call lsjvtensor(nlo_tp%vt,nlovt,j)
        call nlo_tp%add
        nlo_tp%sum=nlo_tp%sum*expexp2

c       n2lo terms 
c       one loop terms  
        call n2lo_ol%init      
        call lsjvcentral(n2lo_ol%vc,n2lo_vc_ol,j)
        call lsjvspinspin(temp1,n2lo_ws_ol,j)
        call isospindependent(temp1,j,n2lo_ol%wss)
        call lsjvtensor(temp1,n2lo_wt_ol,j)
        call isospindependent(temp1,j,n2lo_ol%wt)
        call n2lo_ol%add
        n2lo_ol%sum=n2lo_ol%sum*expexp2
c   sigmak force part
        call lsjvtensor(temp1,pi_gamma,j)
        call isospindependent(temp1,j,temp2)
c        pot=pot+temp2
        pot=pot*expexp2
c   n3lo
        
c      contact terms
        call n3lo_ct%init
        call lsjvcentral(n3lo_ct%vc,n3lo_vc_ct,j)
        call lsjvspinspin(n3lo_ct%vss,n3lo_vs_ct,j)
        call lsjvspinobit(n3lo_ct%vls,n3lo_vls_ct,j)
        call lsjvtensor(n3lo_ct%vt,n3lo_vt_ct,j)
        call lsjvsigmak(n3lo_ct%vsk,n3lo_vsk_ct,j)
        call lsjvsigmal(n3lo_ct%vslsl,n3lo_vslsl_ct,j)
        call n3lo_ct%add
        if(j .eq.0)then
c       D_{3P0}         
         n3lo_ct%sum(3)=n3lo_ct%sum(3)*expexp3
         n3lo_ct%sum(1:2)=n3lo_ct%sum(1:2)*expexp2
         n3lo_ct%sum(4:6)=n3lo_ct%sum(4:6)*expexp2
        else if(j.eq.1)then
c       D_{3P1}         
         n3lo_ct%sum(2)=n3lo_ct%sum(2)*expexp3
         n3lo_ct%sum(1)=n3lo_ct%sum(1)*expexp2
         n3lo_ct%sum(3:6)=n3lo_ct%sum(3:6)*expexp2
        else if(j.eq.2)then
c       D_{1D2}         
         n3lo_ct%sum(1)=n3lo_ct%sum(1)*expexp3
c       D_{3D2}         
         n3lo_ct%sum(2)=n3lo_ct%sum(2)*expexp3
c       D_{3PF2}         
         n3lo_ct%sum(5:6)=n3lo_ct%sum(5:6)*expexp4
         n3lo_ct%sum(3:4)=n3lo_ct%sum(3:4)*expexp2
        else
         n3lo_ct%sum=n3lo_ct%sum*expexp2
        end if



c      pi exchange terms
        call n3lo_rc%init
        call n3lo_fd%init
        call n3lo_tl%init
        call n3lo_cM%init
        call lsjvcentral(n3lo_rc%vc,n3lo_vc_rc,j)
        call lsjvtensor(n3lo_rc%vt,n3lo_vt_rc,j)
        call lsjvspinspin(n3lo_rc%vss,n3lo_vs_rc,j)
        call lsjvspinobit(n3lo_rc%vls,n3lo_vls_rc,j)
        call lsjvcentral(temp1,n3lo_wc_rc,j)
        call isospindependent(temp1,j,n3lo_rc%wc)
        call lsjvtensor(temp1,n3lo_wt_rc,j)
        call isospindependent(temp1,j,n3lo_rc%wt)
        call lsjvspinspin(temp1,n3lo_ws_rc,j)
        call isospindependent(temp1,j,n3lo_rc%wss)
        call lsjvspinobit(temp1,n3lo_wls_rc,j)
        call isospindependent(temp1,j,n3lo_rc%wls)
        call n3lo_rc%add
        n3lo_rc%sum=n3lo_rc%sum*expexp2
        call lsjvcentral(n3lo_fd%vc,n3lo_vc_fd,j)
        call lsjvtensor(temp1,n3lo_wt_fd,j)
        call isospindependent(temp1,j,n3lo_fd%wt)
        call lsjvspinspin(temp1,n3lo_ws_fd,j)
        call isospindependent(temp1,j,n3lo_fd%wss)
        call n3lo_fd%add
        n3lo_fd%sum=n3lo_fd%sum*expexp2
        call lsjvcentral(n3lo_tl%vc,n3lo_vc_tl,j)
        call lsjvspinspin(temp1,n3lo_ws_tl,j)
        call isospindependent(temp1,j,n3lo_tl%wss)
        call lsjvtensor(temp1,n3lo_wt_tl,j)
        call isospindependent(temp1,j,n3lo_tl%wt)
        call lsjvspinspin(n3lo_tl%vss,n3lo_vs_tl,j)
        call lsjvtensor(n3lo_tl%vt,n3lo_vt_tl,j)
        call lsjvcentral(temp1,n3lo_wc_tl,j)
        call isospindependent(temp1,j,n3lo_tl%wc)
        call n3lo_tl%add 
        n3lo_tl%sum=n3lo_tl%sum*expexp2
        call lsjvcentral(n3lo_cM%vc,n3lo_vc_cM,j)
        call lsjvcentral(temp1,n3lo_wc_cM,j)
        call isospindependent(temp1,j,n3lo_cM%wc)
        call lsjvtensor(temp1,n3lo_wt_cM,j)
        call isospindependent(temp1,j,n3lo_cM%wt)
        call lsjvspinspin(temp1,n3lo_ws_cM,j)
        call isospindependent(temp1,j,n3lo_cM%wss)
        call lsjvspinobit(n3lo_cM%vls,n3lo_vls_cM,j)
        call lsjvspinobit(temp1,n3lo_wls_cM,j)
        call isospindependent(temp1,j,n3lo_cM%wls)
        call n3lo_cM%add
        n3lo_cM%sum=n3lo_cM%sum*expexp2
        call lsjvtensor(pigamma,pi_gamma,j)
        if(inn .ne. 2)then
         pigamma=0.0d0
        end if
        pigamma=pigamma*expexp2
         
c     sum all
        pot=lo_ct%sum+lo_onepi%sum+nlo_ct%sum+nlo_tp%sum+n2lo_ol%sum
     + +n3lo_tl%sum+n3lo_fd%sum+n3lo_cM%sum+n3lo_rc%sum+n3lo_ct%sum
     + +pigamma 
        do i=1,6
        pot(i)=pot(i)/(2.0d0*pi)**3*dwnq/ree
        end do
c     set the irelavant terms in nn and np interaction to zero        
        select case(inn)
        case(1)
         if(mod(j,2).eq.0)then
            pot(2)=0.0d0
         else
         pot(1)=0.0d0
         pot(3:6)=0.0d0
         end if
        case(2)
         continue
        case(3)
         if(mod(j,2).eq.0)then
            pot(2)=0.0d0
         else
         pot(1)=0.0d0
         pot(3:6)=0.0d0
         end if
         end select
        return
        end
        subroutine n3lo500new
c
c
c
         use potential_global,only:v,j
         call lsjdecomposition(v,j)
         return
         end subroutine

      subroutine fset(ct,wt,xlb,xub,n)
c     the intergral is sun( wt(i)*f(ct(i))
c       the integrate low bound and up bound        
        real*8 ::xlb,xub
c     the gauss point number        
        integer ::n
c       output  
        real*8 ::ct(96),wt(96)        
c     local
        real*8 ::x(273),a(273)
        integer :: i
c     n=8
      data x(16)/0.960289856497536 d0/, a(16)/0.101228536290376 d0/
      data x(17)/0.796666477413627 d0/, a(17)/0.222381034453374 d0/
      data x(18)/0.525532409916329 d0/, a(18)/0.313706645877887 d0/
      data x(19)/0.183434642495650 d0/, a(19)/0.362683783378362 d0/
c     n=12
      data x(36)/0.981560634246719 d0/, a(36)/0.047175336386512 d0/
      data x(37)/0.904117256370475 d0/, a(37)/0.106939325995318 d0/
      data x(38)/0.769902674194305 d0/, a(38)/0.160078328543346 d0/
      data x(39)/0.587317954286617 d0/, a(39)/0.203167426723066 d0/
      data x(40)/0.367831498998180 d0/, a(40)/0.233492536538355 d0/
      data x(41)/0.125233408511469 d0/, a(41)/0.249147045813403 d0/
c     n=16
      data x(64)/0.989400934991650 d0/, a(64)/0.027152459411754 d0/
      data x(65)/0.944575023073233 d0/, a(65)/0.062253523938648 d0/
      data x(66)/0.865631202387832 d0/, a(66)/0.095158511682493 d0/
      data x(67)/0.755404408355003 d0/, a(67)/0.124628971255534 d0/
      data x(68)/0.617876244402644 d0/, a(68)/0.149595988816577 d0/
      data x(69)/0.458016777657227 d0/, a(69)/0.169156519395003 d0/
      data x(70)/0.281603550779259 d0/, a(70)/0.182603415044924 d0/
      data x(71)/0.095012509837637 d0/, a(71)/0.189450610455069 d0/      
c      n=20
      data x(72)/0.993128599185094 d0/, a(72)/0.017614007139152 d0/
      data x(73)/0.963971927277913 d0/, a(73)/0.040601429800386 d0/
      data x(74)/0.912234428251325 d0/, a(74)/0.062672048334109 d0/
      data x(75)/0.839116971822218 d0/, a(75)/0.083276741576704 d0/
      data x(76)/0.746331906460150 d0/, a(76)/0.101930119817240 d0/
      data x(77)/0.636053680726515 d0/, a(77)/0.118194531961518 d0/
      data x(78)/0.510867001950827 d0/, a(78)/0.131688638449176 d0/
      data x(79)/0.373706088715419 d0/, a(79)/0.142096109318382 d0/
      data x(80)/0.227785851141645 d0/, a(80)/0.149172986472603 d0/
      data x(81)/0.076526521133497 d0/, a(81)/0.152753387130725 d0/        
c       n=24       
        data x(82)/0.995187219997021 d0/, a(82)/0.012341229799987 d0/
        data x(83)/0.974728555971309 d0/, a(83)/0.028531388628933 d0/
        data x(84)/0.938274552002732 d0/, a(84)/0.044277438817419 d0/
        data x(85)/0.886415527004401 d0/, a(85)/0.059298584915436 d0/
        data x(86)/0.820001985973902 d0/, a(86)/0.073346481411080 d0/
        data x(87)/0.740124191578554 d0/, a(87)/0.086190161531953 d0/
        data x(88)/0.648093651936975 d0/, a(88)/0.097618652104113 d0/
        data x(89)/0.545421471388839 d0/, a(89)/0.107444270115965 d0/
        data x(90)/0.433793507626045 d0/, a(90)/0.115505668053725 d0/
        data x(91)/0.315042679696163 d0/, a(91)/0.121670472927803 d0/
        data x(92)/0.191118867473616 d0/, a(92)/0.125837456346828 d0/
        data x(93)/0.064056892862605 d0/, a(93)/0.127938195346752 d0/
c**** n=64
      data x(154)/0.999305041735772d0/, a(154)/0.001783280721696d0/
      data x(155)/0.996340116771955d0/, a(155)/0.004147033260562d0/
      data x(156)/0.991013371476744d0/, a(156)/0.006504457968978d0/
      data x(157)/0.983336253884625d0/, a(157)/0.008846759826363d0/
      data x(158)/0.973326827789910d0/, a(158)/0.011168139460131d0/
      data x(159)/0.961008799652053d0/, a(159)/0.013463047896718d0/
      data x(160)/0.946411374858402d0/, a(160)/0.015726030476024d0/
      data x(161)/0.929569172131939d0/, a(161)/0.017951715775697d0/
      data x(162)/0.910522137078502d0/, a(162)/0.020134823153530d0/
      data x(163)/0.889315445995114d0/, a(163)/0.022270173808383d0/
      data x(164)/0.865999398154092d0/, a(164)/0.024352702568710d0/
      data x(165)/0.840629296252580d0/, a(165)/0.026377469715054d0/
      data x(166)/0.813265315122797d0/, a(166)/0.028339672614259d0/
      data x(167)/0.783972358943341d0/, a(167)/0.030234657072402d0/
      data x(168)/0.752819907260531d0/, a(168)/0.032057928354851d0/
      data x(169)/0.719881850171610d0/, a(169)/0.033805161837141d0/
      data x(170)/0.685236313054233d0/, a(170)/0.035472213256882d0/
      data x(171)/0.648965471254657d0/, a(171)/0.037055128540240d0/
      data x(172)/0.611155355172393d0/, a(172)/0.038550153178615d0/
      data x(173)/0.571895646202634d0/, a(173)/0.039953741132720d0/
      data x(174)/0.531279464019894d0/, a(174)/0.041262563242623d0/
      data x(175)/0.489403145707052d0/, a(175)/0.042473515123653d0/
      data x(176)/0.446366017253464d0/, a(176)/0.043583724529323d0/
      data x(177)/0.402270157963991d0/, a(177)/0.044590558163756d0/
      data x(178)/0.357220158337668d0/, a(178)/0.045491627927418d0/
      data x(179)/0.311322871990210d0/, a(179)/0.046284796581314d0/
      data x(180)/0.264687162208767d0/, a(180)/0.046968182816210d0/
      data x(181)/0.217423643740007d0/, a(181)/0.047540165714830d0/
      data x(182)/0.169644420423992d0/, a(182)/0.047999388596458d0/
      data x(183)/0.121462819296120d0/, a(183)/0.048344762234802d0/
      data x(184)/0.072993121787799d0/, a(184)/0.048575467441503d0/
      data x(185)/0.024350292663424d0/, a(185)/0.048690957009139d0/        
c       n=96        
        data x(226)/0.999689503883230d0/, a(226)/0.000796792065552d0/
        data x(227)/0.998364375863181d0/, a(227)/0.001853960788946d0/
        data x(228)/0.995981842987209d0/, a(228)/0.002910731817934d0/
        data x(229)/0.992543900323762d0/, a(229)/0.003964554338444d0/
        data x(230)/0.988054126329623d0/, a(230)/0.005014202742927d0/
        data x(231)/0.982517263563014d0/, a(231)/0.006058545504235d0/
        data x(232)/0.975939174585136d0/, a(232)/0.007096470791153d0/
        data x(233)/0.968326828463264d0/, a(233)/0.008126876925698d0/
        data x(234)/0.959688291448742d0/, a(234)/0.009148671230783d0/
        data x(235)/0.950032717784437d0/, a(235)/0.010160770535008d0/
        data x(236)/0.939370339752755d0/, a(236)/0.011162102099838d0/
        data x(237)/0.927712456722308d0/, a(237)/0.012151604671088d0/
        data x(238)/0.915071423120898d0/, a(238)/0.013128229566961d0/
        data x(239)/0.901460635315852d0/, a(239)/0.014090941772314d0/
        data x(240)/0.886894517402420d0/, a(240)/0.015038721026994d0/
        data x(241)/0.871388505909296d0/, a(241)/0.015970562902562d0/
        data x(242)/0.854959033434601d0/, a(242)/0.016885479864245d0/
        data x(243)/0.837623511228187d0/, a(243)/0.017782502316045d0/
        data x(244)/0.819400310737931d0/, a(244)/0.018660679627411d0/
        data x(245)/0.800308744139140d0/, a(245)/0.019519081140145d0/
        data x(246)/0.780369043867433d0/, a(246)/0.020356797154333d0/
        data x(247)/0.759602341176647d0/, a(247)/0.021172939892191d0/
        data x(248)/0.738030643744400d0/, a(248)/0.021966644438744d0/
        data x(249)/0.715676812348967d0/, a(249)/0.022737069658329d0/
        data x(250)/0.692564536642171d0/, a(250)/0.023483399085926d0/
        data x(251)/0.668718310043916d0/, a(251)/0.024204841792364d0/
        data x(252)/0.644163403784967d0/, a(252)/0.024900633222483d0/
        data x(253)/0.618925840125468d0/, a(253)/0.025570036005349d0/
        data x(254)/0.593032364777572d0/, a(254)/0.026212340735672d0/
        data x(255)/0.566510418561397d0/, a(255)/0.026826866725591d0/
        data x(256)/0.539388108324357d0/, a(256)/0.027412962726029d0/
        data x(257)/0.511694177154667d0/, a(257)/0.027970007616848d0/
        data x(258)/0.483457973920596d0/, a(258)/0.028497411065085d0/
        data x(259)/0.454709422167743d0/, a(259)/0.028994614150555d0/
        data x(260)/0.425478988407300d0/, a(260)/0.029461089958167d0/ 
        data x(261)/0.395797649828908d0/, a(261)/0.029896344136328d0/ 
        data x(262)/0.365696861472313d0/, a(262)/0.030299915420827d0/ 
        data x(263)/0.335208522892625d0/, a(263)/0.030671376123669d0/ 
        data x(264)/0.304364944354496d0/, a(264)/0.031010332586313d0/ 
        data x(265)/0.273198812591049d0/, a(265)/0.031316425596861d0/ 
        data x(266)/0.241743156163840d0/, a(266)/0.031589330770727d0/ 
        data x(267)/0.210031310460567d0/, a(267)/0.031828758894411d0/ 
        data x(268)/0.178096882367618d0/, a(268)/0.032034456231992d0/ 
        data x(269)/0.145973714654896d0/, a(269)/0.032206204794030d0/
        data x(270)/0.113695850110665d0/, a(270)/0.032343822568575d0/
        data x(271)/0.081297495464425d0/, a(271)/0.032447163714064d0/
        data x(272)/0.048812985136049d0/, a(272)/0.032516118713868d0/
        data x(273)/0.016276744849602d0/, a(273)/0.032550614492363d0/
        wt=0.0d0 
        ct=0.0d0
        select case(n)
        case(8)
         do  i=1,n/2
         wt(i)=a(15+i)   
         ct(i)=-x(15+i)
         end do
         do  i=n,n/2,-1
         ct(i)=-ct(n+1-i)
         wt(i)=wt(n+1-i)
         end do  
      case(12)
         do  i=1,n/2
         wt(i)=a(35+i)   
         ct(i)=-x(35+i)
         end do
         do  i=n,n/2,-1
         ct(i)=-ct(n+1-i)
         wt(i)=wt(n+1-i)
         end do  
      case(16)
         do  i=1,n/2
         wt(i)=a(63+i)   
         ct(i)=-x(63+i)
         end do
         do  i=n,n/2,-1
         ct(i)=-ct(n+1-i)
         wt(i)=wt(n+1-i)
         end do
      case(20)
        do  i=1,n/2
        wt(i)=a(71+i)   
        ct(i)=-x(71+i)
        end do
        do  i=n,n/2,-1
        ct(i)=-ct(n+1-i)
        wt(i)=wt(n+1-i)
        end do
      case(24)
         do  i=1,n/2
         wt(i)=a(81+i)   
         ct(i)=-x(81+i)
         end do
         do  i=n,n/2,-1
         ct(i)=-ct(n+1-i)
         wt(i)=wt(n+1-i)
         end do
      case(64)
         do  i=1,n/2
         wt(i)=a(153+i)   
         ct(i)=-x(153+i)
         end do
         do  i=n,n/2,-1
         ct(i)=-ct(n+1-i)
         wt(i)=wt(n+1-i)
         end do
      case(96)
         do  i=1,n/2
         wt(i)=a(225+i)   
         ct(i)=-x(225+i)
         end do
         do  i=n,n/2,-1
         ct(i)=-ct(n+1-i)
         wt(i)=wt(n+1-i)
         end do
      end select
        do i=1,96
         ct(i)=xlb+(ct(i)+1.0d0)*(xub-xlb)/2.0d0
        end do
        wt=wt*(xub-xlb)/2.0d0
        return             
      end subroutine
            subroutine gset(ax,bx,n,z,w)
c
c
c        this code has been obtained from the CERN computer library
c        in the year of the lord 1972.
c
c
      implicit real*8 (a-h,o-z)
c
c     n-point gauss zeros and weights for the interval (ax,bx) are
c           stored in  arrays z and w respectively.
c
      dimension     a(273),x(273),ktab(96)
      dimension z(1),w(1)
c
c-----table of initial subscripts for n=2(1)16(4)96
      data ktab(2)/1/
      data ktab(3)/2/
      data ktab(4)/4/
      data ktab(5)/6/
      data ktab(6)/9/
      data ktab(7)/12/
      data ktab(8)/16/
      data ktab(9)/20/
      data ktab(10)/25/
      data ktab(11)/30/
      data ktab(12)/36/
      data ktab(13)/42/
      data ktab(14)/49/
      data ktab(15)/56/
      data ktab(16)/64/
      data ktab(20)/72/
      data ktab(24)/82/
      data ktab(28)/82/
      data ktab(32)/94/
      data ktab(36)/94/
      data ktab(40)/110/
      data ktab(44)/110/
      data ktab(48)/130/
      data ktab(52)/130/
      data ktab(56)/130/
      data ktab(60)/130/
      data ktab(64)/154/
      data ktab(68)/154/
      data ktab(72)/154/
      data ktab(76)/154/
      data ktab(80)/186/
      data ktab(84)/186/
      data ktab(88)/186/
      data ktab(92)/186/
      data ktab(96)/226/
c
c-----table of abscissae (x) and weights (a) for interval (-1,+1).
c
c**** n=2
      data x(1)/0.577350269189626  d0/, a(1)/1.000000000000000  d0/
c**** n=3
      data x(2)/0.774596669241483  d0/, a(2)/0.555555555555556  d0/
      data x(3)/0.000000000000000  d0/, a(3)/0.888888888888889  d0/
c**** n=4
      data x(4)/0.861136311594053  d0/, a(4)/0.347854845137454  d0/
      data x(5)/0.339981043584856  d0/, a(5)/0.652145154862546  d0/
c**** n=5
      data x(6)/0.906179845938664  d0/, a(6)/0.236926885056189  d0/
      data x(7)/0.538469310105683  d0/, a(7)/0.478628670499366  d0/
      data x(8)/0.000000000000000  d0/, a(8)/0.568888888888889  d0/
c**** n=6
      data x(9)/0.932469514203152  d0/, a(9)/0.171324492379170  d0/
      data x(10)/0.661209386466265 d0/, a(10)/0.360761573048139 d0/
      data x(11)/0.238619186083197 d0/, a(11)/0.467913934572691 d0/
c**** n=7
      data x(12)/0.949107912342759 d0/, a(12)/0.129484966168870 d0/
      data x(13)/0.741531185599394 d0/, a(13)/0.279705391489277 d0/
      data x(14)/0.405845151377397 d0/, a(14)/0.381830050505119 d0/
      data x(15)/0.000000000000000 d0/, a(15)/0.417959183673469 d0/
c**** n=8
      data x(16)/0.960289856497536 d0/, a(16)/0.101228536290376 d0/
      data x(17)/0.796666477413627 d0/, a(17)/0.222381034453374 d0/
      data x(18)/0.525532409916329 d0/, a(18)/0.313706645877887 d0/
      data x(19)/0.183434642495650 d0/, a(19)/0.362683783378362 d0/
c**** n=9
      data x(20)/0.968160239507626 d0/, a(20)/0.081274388361574 d0/
      data x(21)/0.836031107326636 d0/, a(21)/0.180648160694857 d0/
      data x(22)/0.613371432700590 d0/, a(22)/0.260610696402935 d0/
      data x(23)/0.324253423403809 d0/, a(23)/0.312347077040003 d0/
      data x(24)/0.000000000000000 d0/, a(24)/0.330239355001260 d0/
c**** n=10
      data x(25)/0.973906528517172 d0/, a(25)/0.066671344308688 d0/
      data x(26)/0.865063366688985 d0/, a(26)/0.149451349150581 d0/
      data x(27)/0.679409568299024 d0/, a(27)/0.219086362515982 d0/
      data x(28)/0.433395394129247 d0/, a(28)/0.269266719309996 d0/
      data x(29)/0.148874338981631 d0/, a(29)/0.295524224714753 d0/
c**** n=11
      data x(30)/0.978228658146057 d0/, a(30)/0.055668567116174 d0/
      data x(31)/0.887062599768095 d0/, a(31)/0.125580369464905 d0/
      data x(32)/0.730152005574049 d0/, a(32)/0.186290210927734 d0/
      data x(33)/0.519096129206812 d0/, a(33)/0.233193764591990 d0/
      data x(34)/0.269543155952345 d0/, a(34)/0.262804544510247 d0/
      data x(35)/0.000000000000000 d0/, a(35)/0.272925086777901 d0/
c**** n=12
      data x(36)/0.981560634246719 d0/, a(36)/0.047175336386512 d0/
      data x(37)/0.904117256370475 d0/, a(37)/0.106939325995318 d0/
      data x(38)/0.769902674194305 d0/, a(38)/0.160078328543346 d0/
      data x(39)/0.587317954286617 d0/, a(39)/0.203167426723066 d0/
      data x(40)/0.367831498998180 d0/, a(40)/0.233492536538355 d0/
      data x(41)/0.125233408511469 d0/, a(41)/0.249147045813403 d0/
c**** n=13
      data x(42)/0.984183054718588 d0/, a(42)/0.040484004765316 d0/
      data x(43)/0.917598399222978 d0/, a(43)/0.092121499837728 d0/
      data x(44)/0.801578090733310 d0/, a(44)/0.138873510219787 d0/
      data x(45)/0.642349339440340 d0/, a(45)/0.178145980761946 d0/
      data x(46)/0.448492751036447 d0/, a(46)/0.207816047536889 d0/
      data x(47)/0.230458315955135 d0/, a(47)/0.226283180262897 d0/
      data x(48)/0.000000000000000 d0/, a(48)/0.232551553230874 d0/
c**** n=14
      data x(49)/0.986283808696812 d0/, a(49)/0.035119460331752 d0/
      data x(50)/0.928434883663574 d0/, a(50)/0.080158087159760 d0/
      data x(51)/0.827201315069765 d0/, a(51)/0.121518570687903 d0/
      data x(52)/0.687292904811685 d0/, a(52)/0.157203167158194 d0/
      data x(53)/0.515248636358154 d0/, a(53)/0.185538397477938 d0/
      data x(54)/0.319112368927890 d0/, a(54)/0.205198463721296 d0/
      data x(55)/0.108054948707344 d0/, a(55)/0.215263853463158 d0/
c**** n=15
      data x(56)/0.987992518020485 d0/, a(56)/0.030753241996117 d0/
      data x(57)/0.937273392400706 d0/, a(57)/0.070366047488108 d0/
      data x(58)/0.848206583410427 d0/, a(58)/0.107159220467172 d0/
      data x(59)/0.724417731360170 d0/, a(59)/0.139570677926154 d0/
      data x(60)/0.570972172608539 d0/, a(60)/0.166269205816994 d0/
      data x(61)/0.394151347077563 d0/, a(61)/0.186161000015562 d0/
      data x(62)/0.201194093997435 d0/, a(62)/0.198431485327111 d0/
      data x(63)/0.000000000000000 d0/, a(63)/0.202578241925561 d0/
c**** n=16
      data x(64)/0.989400934991650 d0/, a(64)/0.027152459411754 d0/
      data x(65)/0.944575023073233 d0/, a(65)/0.062253523938648 d0/
      data x(66)/0.865631202387832 d0/, a(66)/0.095158511682493 d0/
      data x(67)/0.755404408355003 d0/, a(67)/0.124628971255534 d0/
      data x(68)/0.617876244402644 d0/, a(68)/0.149595988816577 d0/
      data x(69)/0.458016777657227 d0/, a(69)/0.169156519395003 d0/
      data x(70)/0.281603550779259 d0/, a(70)/0.182603415044924 d0/
      data x(71)/0.095012509837637 d0/, a(71)/0.189450610455069 d0/
c**** n=20
      data x(72)/0.993128599185094 d0/, a(72)/0.017614007139152 d0/
      data x(73)/0.963971927277913 d0/, a(73)/0.040601429800386 d0/
      data x(74)/0.912234428251325 d0/, a(74)/0.062672048334109 d0/
      data x(75)/0.839116971822218 d0/, a(75)/0.083276741576704 d0/
      data x(76)/0.746331906460150 d0/, a(76)/0.101930119817240 d0/
      data x(77)/0.636053680726515 d0/, a(77)/0.118194531961518 d0/
      data x(78)/0.510867001950827 d0/, a(78)/0.131688638449176 d0/
      data x(79)/0.373706088715419 d0/, a(79)/0.142096109318382 d0/
      data x(80)/0.227785851141645 d0/, a(80)/0.149172986472603 d0/
      data x(81)/0.076526521133497 d0/, a(81)/0.152753387130725 d0/
c**** n=24
      data x(82)/0.995187219997021 d0/, a(82)/0.012341229799987 d0/
      data x(83)/0.974728555971309 d0/, a(83)/0.028531388628933 d0/
      data x(84)/0.938274552002732 d0/, a(84)/0.044277438817419 d0/
      data x(85)/0.886415527004401 d0/, a(85)/0.059298584915436 d0/
      data x(86)/0.820001985973902 d0/, a(86)/0.073346481411080 d0/
      data x(87)/0.740124191578554 d0/, a(87)/0.086190161531953 d0/
      data x(88)/0.648093651936975 d0/, a(88)/0.097618652104113 d0/
      data x(89)/0.545421471388839 d0/, a(89)/0.107444270115965 d0/
      data x(90)/0.433793507626045 d0/, a(90)/0.115505668053725 d0/
      data x(91)/0.315042679696163 d0/, a(91)/0.121670472927803 d0/
      data x(92)/0.191118867473616 d0/, a(92)/0.125837456346828 d0/
      data x(93)/0.064056892862605 d0/, a(93)/0.127938195346752 d0/
c**** n=32
      data x(94)/0.997263861849481 d0/, a(94)/0.007018610009470 d0/
      data x(95)/0.985611511545268 d0/, a(95)/0.016274394730905 d0/
      data x(96)/0.964762255587506 d0/, a(96)/0.025392065309262 d0/
      data x(97)/0.934906075937739 d0/, a(97)/0.034273862913021 d0/
      data x(98)/0.896321155766052 d0/, a(98)/0.042835898022226 d0/
      data x(99)/0.849367613732569 d0/, a(99)/0.050998059262376 d0/
      data x(100)/0.794483795967942d0/, a(100)/0.058684093478535d0/
      data x(101)/0.732182118740289d0/, a(101)/0.065822222776361d0/
      data x(102)/0.663044266930215d0/, a(102)/0.072345794108848d0/
      data x(103)/0.587715757240762d0/, a(103)/0.078193895787070d0/
      data x(104)/0.506899908932229d0/, a(104)/0.083311924226946d0/
      data x(105)/0.421351276130635d0/, a(105)/0.087652093004403d0/
      data x(106)/0.331868602282127d0/, a(106)/0.091173878695763d0/
      data x(107)/0.239287362252137d0/, a(107)/0.093844399080804d0/
      data x(108)/0.144471961582796d0/, a(108)/0.095638720079274d0/
      data x(109)/0.048307665687738d0/, a(109)/0.096540088514727d0/
c**** n=40
      data x(110)/0.998237709710559d0/, a(110)/0.004521277098533d0/
      data x(111)/0.990726238699457d0/, a(111)/0.010498284531152d0/
      data x(112)/0.977259949983774d0/, a(112)/0.016421058381907d0/
      data x(113)/0.957916819213791d0/, a(113)/0.022245849194166d0/
      data x(114)/0.932812808278676d0/, a(114)/0.027937006980023d0/
      data x(115)/0.902098806968874d0/, a(115)/0.033460195282547d0/
      data x(116)/0.865959503212259d0/, a(116)/0.038782167974472d0/
      data x(117)/0.824612230833311d0/, a(117)/0.043870908185673d0/
      data x(118)/0.778305651426519d0/, a(118)/0.048695807635072d0/
      data x(119)/0.727318255189927d0/, a(119)/0.053227846983936d0/
      data x(120)/0.671956684614179d0/, a(120)/0.057439769099391d0/
      data x(121)/0.612553889667980d0/, a(121)/0.061306242492928d0/
      data x(122)/0.549467125095128d0/, a(122)/0.064804013456601d0/
      data x(123)/0.483075801686178d0/, a(123)/0.067912045815233d0/
      data x(124)/0.413779204371605d0/, a(124)/0.070611647391286d0/
      data x(125)/0.341994090825758d0/, a(125)/0.072886582395804d0/
      data x(126)/0.268152185007253d0/, a(126)/0.074723169057968d0/
      data x(127)/0.192697580701371d0/, a(127)/0.076110361900626d0/
      data x(128)/0.116084070675255d0/, a(128)/0.077039818164247d0/
      data x(129)/0.038772417506050d0/, a(129)/0.077505947978424d0/
c**** n=48
      data x(130)/0.998771007252426d0/, a(130)/0.003153346052305d0/
      data x(131)/0.993530172266350d0/, a(131)/0.007327553901276d0/
      data x(132)/0.984124583722826d0/, a(132)/0.011477234579234d0/
      data x(133)/0.970591592546247d0/, a(133)/0.015579315722943d0/
      data x(134)/0.952987703160430d0/, a(134)/0.019616160457355d0/
      data x(135)/0.931386690706554d0/, a(135)/0.023570760839324d0/
      data x(136)/0.905879136715569d0/, a(136)/0.027426509708356d0/
      data x(137)/0.876572020274247d0/, a(137)/0.031167227832798d0/
      data x(138)/0.843588261624393d0/, a(138)/0.034777222564770d0/
      data x(139)/0.807066204029442d0/, a(139)/0.038241351065830d0/
      data x(140)/0.767159032515740d0/, a(140)/0.041545082943464d0/
      data x(141)/0.724034130923814d0/, a(141)/0.044674560856694d0/
      data x(142)/0.677872379632663d0/, a(142)/0.047616658492490d0/
      data x(143)/0.628867396776513d0/, a(143)/0.050359035553854d0/
      data x(144)/0.577224726083972d0/, a(144)/0.052890189485193d0/
      data x(145)/0.523160974722233d0/, a(145)/0.055199503699984d0/
      data x(146)/0.466902904750958d0/, a(146)/0.057277292100403d0/
      data x(147)/0.408686481990716d0/, a(147)/0.059114839698395d0/
      data x(148)/0.348755886292160d0/, a(148)/0.060704439165893d0/
      data x(149)/0.287362487355455d0/, a(149)/0.062039423159892d0/
      data x(150)/0.224763790394689d0/, a(150)/0.063114192286254d0/
      data x(151)/0.161222356068891d0/, a(151)/0.063924238584648d0/
      data x(152)/0.097004699209462d0/, a(152)/0.064466164435950d0/
      data x(153)/0.032380170962869d0/, a(153)/0.064737696812683d0/
c**** n=64
      data x(154)/0.999305041735772d0/, a(154)/0.001783280721696d0/
      data x(155)/0.996340116771955d0/, a(155)/0.004147033260562d0/
      data x(156)/0.991013371476744d0/, a(156)/0.006504457968978d0/
      data x(157)/0.983336253884625d0/, a(157)/0.008846759826363d0/
      data x(158)/0.973326827789910d0/, a(158)/0.011168139460131d0/
      data x(159)/0.961008799652053d0/, a(159)/0.013463047896718d0/
      data x(160)/0.946411374858402d0/, a(160)/0.015726030476024d0/
      data x(161)/0.929569172131939d0/, a(161)/0.017951715775697d0/
      data x(162)/0.910522137078502d0/, a(162)/0.020134823153530d0/
      data x(163)/0.889315445995114d0/, a(163)/0.022270173808383d0/
      data x(164)/0.865999398154092d0/, a(164)/0.024352702568710d0/
      data x(165)/0.840629296252580d0/, a(165)/0.026377469715054d0/
      data x(166)/0.813265315122797d0/, a(166)/0.028339672614259d0/
      data x(167)/0.783972358943341d0/, a(167)/0.030234657072402d0/
      data x(168)/0.752819907260531d0/, a(168)/0.032057928354851d0/
      data x(169)/0.719881850171610d0/, a(169)/0.033805161837141d0/
      data x(170)/0.685236313054233d0/, a(170)/0.035472213256882d0/
      data x(171)/0.648965471254657d0/, a(171)/0.037055128540240d0/
      data x(172)/0.611155355172393d0/, a(172)/0.038550153178615d0/
      data x(173)/0.571895646202634d0/, a(173)/0.039953741132720d0/
      data x(174)/0.531279464019894d0/, a(174)/0.041262563242623d0/
      data x(175)/0.489403145707052d0/, a(175)/0.042473515123653d0/
      data x(176)/0.446366017253464d0/, a(176)/0.043583724529323d0/
      data x(177)/0.402270157963991d0/, a(177)/0.044590558163756d0/
      data x(178)/0.357220158337668d0/, a(178)/0.045491627927418d0/
      data x(179)/0.311322871990210d0/, a(179)/0.046284796581314d0/
      data x(180)/0.264687162208767d0/, a(180)/0.046968182816210d0/
      data x(181)/0.217423643740007d0/, a(181)/0.047540165714830d0/
      data x(182)/0.169644420423992d0/, a(182)/0.047999388596458d0/
      data x(183)/0.121462819296120d0/, a(183)/0.048344762234802d0/
      data x(184)/0.072993121787799d0/, a(184)/0.048575467441503d0/
      data x(185)/0.024350292663424d0/, a(185)/0.048690957009139d0/
c**** n=80
      data x(186)/0.999553822651630d0/, a(186)/0.001144950003186d0/
      data x(187)/0.997649864398237d0/, a(187)/0.002663533589512d0/
      data x(188)/0.994227540965688d0/, a(188)/0.004180313124694d0/
      data x(189)/0.989291302499755d0/, a(189)/0.005690922451403d0/
      data x(190)/0.982848572738629d0/, a(190)/0.007192904768117d0/
      data x(191)/0.974909140585727d0/, a(191)/0.008683945269260d0/
      data x(192)/0.965485089043799d0/, a(192)/0.010161766041103d0/
      data x(193)/0.954590766343634d0/, a(193)/0.011624114120797d0/
      data x(194)/0.942242761309872d0/, a(194)/0.013068761592401d0/
      data x(195)/0.928459877172445d0/, a(195)/0.014493508040509d0/
      data x(196)/0.913263102571757d0/, a(196)/0.015896183583725d0/
      data x(197)/0.896675579438770d0/, a(197)/0.017274652056269d0/
      data x(198)/0.878722567678213d0/, a(198)/0.018626814208299d0/
      data x(199)/0.859431406663111d0/, a(199)/0.019950610878141d0/
      data x(200)/0.838831473580255d0/, a(200)/0.021244026115782d0/
      data x(201)/0.816954138681463d0/, a(201)/0.022505090246332d0/
      data x(202)/0.793832717504605d0/, a(202)/0.023731882865930d0/
      data x(203)/0.769502420135041d0/, a(203)/0.024922535764115d0/
      data x(204)/0.744000297583597d0/, a(204)/0.026075235767565d0/
      data x(205)/0.717365185362099d0/, a(205)/0.027188227500486d0/
      data x(206)/0.689637644342027d0/, a(206)/0.028259816057276d0/
      data x(207)/0.660859898986119d0/, a(207)/0.029288369583267d0/
      data x(208)/0.631075773046871d0/, a(208)/0.030272321759557d0/
      data x(209)/0.600330622829751d0/, a(209)/0.031210174188114d0/
      data x(210)/0.568671268122709d0/, a(210)/0.032100498673487d0/
      data x(211)/0.536145920897131d0/, a(211)/0.032941939397645d0/
      data x(212)/0.502804111888784d0/, a(212)/0.033733214984611d0/
      data x(213)/0.468696615170544d0/, a(213)/0.034473120451753d0/
      data x(214)/0.433875370831756d0/, a(214)/0.035160529044747d0/
      data x(215)/0.398393405881969d0/, a(215)/0.035794393953416d0/
      data x(216)/0.362304753499487d0/, a(216)/0.036373749905835d0/
      data x(217)/0.325664370747701d0/, a(217)/0.036897714638276d0/
      data x(218)/0.288528054884511d0/, a(218)/0.037365490238730d0/
      data x(219)/0.250952358392272d0/, a(219)/0.037776364362001d0/
      data x(220)/0.212994502857666d0/, a(220)/0.038129711314477d0/
      data x(221)/0.174712291832646d0/, a(221)/0.038424993006959d0/
      data x(222)/0.136164022809143d0/, a(222)/0.038661759774076d0/
      data x(223)/0.097408398441584d0/, a(223)/0.038839651059051d0/
      data x(224)/0.058504437152420d0/, a(224)/0.038958395962769d0/
      data x(225)/0.019511383256793d0/, a(225)/0.039017813656306d0/
c**** n=96
      data x(226)/0.999689503883230d0/, a(226)/0.000796792065552d0/
      data x(227)/0.998364375863181d0/, a(227)/0.001853960788946d0/
      data x(228)/0.995981842987209d0/, a(228)/0.002910731817934d0/
      data x(229)/0.992543900323762d0/, a(229)/0.003964554338444d0/
      data x(230)/0.988054126329623d0/, a(230)/0.005014202742927d0/
      data x(231)/0.982517263563014d0/, a(231)/0.006058545504235d0/
      data x(232)/0.975939174585136d0/, a(232)/0.007096470791153d0/
      data x(233)/0.968326828463264d0/, a(233)/0.008126876925698d0/
      data x(234)/0.959688291448742d0/, a(234)/0.009148671230783d0/
      data x(235)/0.950032717784437d0/, a(235)/0.010160770535008d0/
      data x(236)/0.939370339752755d0/, a(236)/0.011162102099838d0/
      data x(237)/0.927712456722308d0/, a(237)/0.012151604671088d0/
      data x(238)/0.915071423120898d0/, a(238)/0.013128229566961d0/
      data x(239)/0.901460635315852d0/, a(239)/0.014090941772314d0/
      data x(240)/0.886894517402420d0/, a(240)/0.015038721026994d0/
      data x(241)/0.871388505909296d0/, a(241)/0.015970562902562d0/
      data x(242)/0.854959033434601d0/, a(242)/0.016885479864245d0/
      data x(243)/0.837623511228187d0/, a(243)/0.017782502316045d0/
      data x(244)/0.819400310737931d0/, a(244)/0.018660679627411d0/
      data x(245)/0.800308744139140d0/, a(245)/0.019519081140145d0/
      data x(246)/0.780369043867433d0/, a(246)/0.020356797154333d0/
      data x(247)/0.759602341176647d0/, a(247)/0.021172939892191d0/
      data x(248)/0.738030643744400d0/, a(248)/0.021966644438744d0/
      data x(249)/0.715676812348967d0/, a(249)/0.022737069658329d0/
      data x(250)/0.692564536642171d0/, a(250)/0.023483399085926d0/
      data x(251)/0.668718310043916d0/, a(251)/0.024204841792364d0/
      data x(252)/0.644163403784967d0/, a(252)/0.024900633222483d0/
      data x(253)/0.618925840125468d0/, a(253)/0.025570036005349d0/
      data x(254)/0.593032364777572d0/, a(254)/0.026212340735672d0/
      data x(255)/0.566510418561397d0/, a(255)/0.026826866725591d0/
      data x(256)/0.539388108324357d0/, a(256)/0.027412962726029d0/
      data x(257)/0.511694177154667d0/, a(257)/0.027970007616848d0/
      data x(258)/0.483457973920596d0/, a(258)/0.028497411065085d0/
      data x(259)/0.454709422167743d0/, a(259)/0.028994614150555d0/
      data x(260)/0.425478988407300d0/, a(260)/0.029461089958167d0/
      data x(261)/0.395797649828908d0/, a(261)/0.029896344136328d0/
      data x(262)/0.365696861472313d0/, a(262)/0.030299915420827d0/
      data x(263)/0.335208522892625d0/, a(263)/0.030671376123669d0/
      data x(264)/0.304364944354496d0/, a(264)/0.031010332586313d0/
      data x(265)/0.273198812591049d0/, a(265)/0.031316425596861d0/
      data x(266)/0.241743156163840d0/, a(266)/0.031589330770727d0/
      data x(267)/0.210031310460567d0/, a(267)/0.031828758894411d0/
      data x(268)/0.178096882367618d0/, a(268)/0.032034456231992d0/
      data x(269)/0.145973714654896d0/, a(269)/0.032206204794030d0/
      data x(270)/0.113695850110665d0/, a(270)/0.032343822568575d0/
      data x(271)/0.081297495464425d0/, a(271)/0.032447163714064d0/
      data x(272)/0.048812985136049d0/, a(272)/0.032516118713868d0/
      data x(273)/0.016276744849602d0/, a(273)/0.032550614492363d0/
c
c
c-----test n
      alpha=0.5d0*(ax+bx)
      beta=0.5d0*(bx-ax)
      if( n.lt.1 .or. n.gt.96 ) go to 100
      if(n.ne.1) go to 1
      z(1)=alpha
      w(1)=bx-ax
      return
c
    1 if (n.le.16) go to 3
      if (n.gt.24) go to 4
      n=4*(n/4)
      go to 3
    4 if (n.gt.48) go to 5
      n=8*(n/8)
      go to 3
    5 n=16*(n/16)
c
c----- set k equal to initial subscript and store results
    3 k=ktab(n)
      m=n/2
      do 2 j=1,m
      jtab=k-1+j
      wtemp=beta*a(jtab)
      delta=beta*x(jtab)
      z(j)=alpha-delta
      w(j)=wtemp
      jp=n+1-j
      z(jp)=alpha+delta
      w(jp)=wtemp
    2 continue
      if((n-m-m).eq.0) return
      z(m+1)=alpha
      jmid=k+m
      w(m+1)=beta*a(jmid)
      return
c
  100 zn=n
      write(6,200) zn
  200 format(/////' error in gset. n has the non-permissible value',
     1e11.3/' execution terminated.')
      stop
      end
c name:    dgelg
c        programmbibliothek rhrz bonn        02/02/81       dgelg
c                                            fortran iv     ibm 370/168
c
c purpose:
c
c to solve a general system of simultaneous linear equations.
c
c usage:   call dgelg(r,a,m,n,eps,ier)
c
c parameters:
c
c r:       double precision m by n right hand side matrix
c          (destroyed). on return r contains the solutions
c          of the equations.
c
c a:       double precision m by m coefficient matrix
c          (destroyed).
c
c m:       the number of equations in the system.
c
c n:       the number of right hand side vectors.
c
c eps:     single precision input constant which is used as
c          relative tolerance for test on loss of
c          significance.
c
c ier:     resulting error parameter coded as follows
c           ier=0  - no error,
c           ier=-1 - no result because of m less than 1 or
c                   pivot element at any elimination step
c                   equal to 0,
c           ier=k  - warning due to possible loss of signifi-
c                   cance indicated at elimination step k+1,
c                   where pivot element was less than or
c                   equal to the internal tolerance eps times
c                   absolutely greatest element of matrix a.
c
c remarks: (1) input matrices r and a are assumed to be stored
c              columnwise in m*n resp. m*m successive storage
c              locations. on return solution matrix r is stored
c              columnwise too.
c          (2) the procedure gives results if the number of equations m
c              is greater than 0 and pivot elements at all elimination
c              steps are different from 0. however warning ier=k - if
c              given indicates possible loss of significance. in case
c              of a well scaled matrix a and appropriate tolerance eps,
c              ier=k may be interpreted that matrix a has the rank k.
c              no warning is given in case m=1.
c
c method:
c
c solution is done by means of gauss-elimination with
c complete pivoting.
c
c programs required:
c          none
c
c access:
c
c load module:    sys3.fortlib(dgelg)
c source module:  sys3.symlib.fortran(dgelg)
c description:    sys3.infolib(dgelg)
c
c author:         ibm, ssp iii
c installation:   ibm 370/168, mvs-jes2, fortran iv (h ext. enh.)
c
c**********************************************************************
      subroutine dgelg(r,a,m,n,eps,ier)
c
c
      implicit real*8 (a-h,o-z)
      dimension a(1),r(1)
      real*4 eps
c
c
c
c
      if(m)23,23,1
c
c     search for greatest element in matrix a
    1 ier=0
      piv=0.d0
      mm=m*m
      nm=n*m
      do 3 l=1,mm
      tb=dabs(a(l))
      if(tb-piv)3,3,2
    2 piv=tb
      i=l
    3 continue
      tol=eps*piv
c     a(i) is pivot element. piv contains the absolute value of a(i).
c
c
c     start elimination loop
      lst=1
      do 17 k=1,m
c
c     test on singularity
      if(piv)23,23,4
    4 if(ier)7,5,7
    5 if(piv-tol)6,6,7
    6 ier=k-1
    7 pivi=1.d0/a(i)
      j=(i-1)/m
      i=i-j*m-k
      j=j+1-k
c     i+k is row-index, j+k column-index of pivot element
c
c     pivot row reduction and row interchange in right hand side r
      do 8 l=k,nm,m
      ll=l+i
      tb=pivi*r(ll)
      r(ll)=r(l)
    8 r(l)=tb
c
c     is elimination terminated
      if(k-m)9,18,18
c
c     column interchange in matrix a
    9 lend=lst+m-k
      if(j)12,12,10
   10 ii=j*m
      do 11 l=lst,lend
      tb=a(l)
      ll=l+ii
      a(l)=a(ll)
   11 a(ll)=tb
c
c     row interchange and pivot row reduction in matrix a
   12 do 13 l=lst,mm,m
      ll=l+i
      tb=pivi*a(ll)
      a(ll)=a(l)
   13 a(l)=tb
c
c     save column interchange information
      a(lst)=j
c
c     element reduction and next pivot search
      piv=0.d0
      lst=lst+1
      j=0
      do 16 ii=lst,lend
      pivi=-a(ii)
      ist=ii+m
      j=j+1
      do 15 l=ist,mm,m
      ll=l-j
      a(l)=a(l)+pivi*a(ll)
      tb=dabs(a(l))
      if(tb-piv)15,15,14
   14 piv=tb
      i=l
   15 continue
      do 16 l=k,nm,m
      ll=l+j
   16 r(ll)=r(ll)+pivi*r(l)
   17 lst=lst+m
c     end of elimination loop
c
c
c     back substitution and back interchange
   18 if(m-1)23,22,19
   19 ist=mm+m
      lst=m+1
      do 21 i=2,m
      ii=lst-i
      ist=ist-lst
      l=ist-m
      l=a(l)+.5d0
      do 21 j=ii,nm,m
      tb=r(j)
      ll=j
      do 20 k=ist,mm,m
      ll=ll+1
   20 tb=tb-a(k)*r(ll)
      k=j+l
      r(j)=r(k)
   21 r(k)=tb
   22 return
c
c
c     error return
   23 ier=-1
      return
      end