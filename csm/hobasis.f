!!hobasis (complex)
        module hobasis
            use parameter
            contains
!!norm of ho basis
            real*8 function normho(nu,n,l)
            implicit none
            real*8::nu
            integer::n,l
            normho=sqrt(sqrt(2*nu**3/pi)*2**(n+2*l+3)*fact(n)
     &      *nu**l/doublefact(2*n+2*l+1))
      !   if(abs(normho)<1e-6) then
      !      write(*,*)'nu,n,l,normho',nu,n,l,normho
      !      write(*,*)'fact(n)',fact(n)
      !      write(*,*)'fact(n+l)=',fact(n+l)
      !      write(*,*)'fact(2n+2l+1)=',fact(2*n+2*l+1)
      !      stop
      !   endif
            end function
!!hobasis
            complex*16 function ho3d(n,l,nu,r)            !!标准形式的三维BASIS
            implicit none                           !!修改后的nu取mu*omega/(2hbar)
            integer l,n
            real*8::norma,nu
            complex*16::r
            norma=normho(nu,n,l)
      !        if (norma<1e-6) then
      !           write(*,*)'ho3d: Norm=0!!!for  nu,n,l',nu,n,l
      !        endif
            ho3d=norma*r**l*exp(-nu*r**2)*
     &            general_laguerre(2*nu*r**2,n,l+0.5d0) 
            end function ho3d
!! fact            
            function fact(n)           
                  implicit none
                  integer n
                  real*8 fact, dgamma,x
                  x=dfloat(n+1)
                  fact=dgamma(x)
            end function fact
!!double fact
      real*8 function doublefact(n)!双阶乘
            implicit none
            integer::n,i
            real::s
            s=1.0
            if(mod(n,2)==0) then
                  do i=n,2,-2
                        s=s*i
                  end do
            else
                  do i=n,1,-2
                        s=s*i
                  end do
            end if
            doublefact=s
      end function
!!Laguerre function
      recursive function general_laguerre(x, k, alpha) result(s)        !!2*nu*r**2,n,l+0.5d0
      complex*16, intent(in) :: x
      integer, intent(in) :: k
      real*8, intent(in) :: alpha
      complex*16 :: s

      if (k == 0) then
         s = 1.0
      else if (k == 1) then
         s = 1.0 + alpha - x
      else
         s = ((2*k - 1 + alpha - x) * general_laguerre(x, k-1, alpha) 
     &    - (k - 1 + alpha) * general_laguerre(x, k-2, alpha)) / k
      end if
      end function general_laguerre
        end module
