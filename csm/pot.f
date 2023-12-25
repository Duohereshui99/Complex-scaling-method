ccccccc
        module potentialfunctions
          use parameter
            contains
!!gauss potential function
            complex*16 function gausspot(r,v0,r0,a0)
                implicit none
                 real*8::v0,r0,a0
                 complex*16::r
                   if (a0.gt.1e-6) then
                     gausspot=V0*exp(-(r-r0)**2/a0**2)
                       else
                        write(*,*)'a too small in gausspot!'
                        stop
                   endif
                   return
            end function
ccccccc
        complex*16 function vpotnn(v0,a0,r0,r)
        real*8::v0,a0,r0
        complex*16::r
        vpotnn=-v0*(1d0+cosh(r0/a0))/(cosh(r/a0)+cosh(r0/a0))
        end function vpotnn
!
!!Coulomb potential: z12=z1*z2, radius parameter r0
!
        complex*16 function vpotcoul(z12,r0,r)
        real*8::z12,r0
        complex*16::r
        if (abs(r)>=r0) then
            vpotcoul=z12*e2/r
        else
            vpotcoul=z12*e2/2/r0*(3d0-r**2/r0**2)
        end if 
        end function vpotcoul
!
!!Centrifugal barrier potential: orbital angular momentum l, reduced mass mu
!
        complex*16 function vpotcent(mu,l,r)
        real*8::mu
        integer::l
        complex*16::r
        vpotcent=hbarc**2/2d0/mu*(l+0.5d0)**2/r**2
        end function vpotcent
!
!!alpha core potential in the cluster model
!
        complex*16 function vpots(v0,a0,r0,z12,mu,l,r)
        real*8::v0,a0,r0,z12,mu
        integer::l
        complex*16::r
        vpots=vpotnn(v0,a0,r0,r)+vpotcoul(z12,r0,r)
   !  &   +vpotcent(mu,l,r)
        end function vpots
ccccccc
        function potencc(r)
      implicit none
      complex*16,intent(in)                :: r
      complex*16                           :: potencc
      potencc=(-626.885D0*exp(-1.55D0*r)+1438.72d0*exp(-3.11D0*r))/r
      !potencc=-258.7D0*exp(-(r/2.0d0)*(r/2.0d0))
      end function
        end module