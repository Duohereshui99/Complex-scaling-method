ccccccc
      module algorithm
        contains
!interpolation function for uniform grids (real functions)
!Y: value of function's index Y  
!F: function value array, size of the array
      function FFR4(Y,F,N)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 F(N),P,P1,P2,Q,X,FFR4
      REAL*8 Y
      PARAMETER(X=.16666666666667)
      P=Y
      I=P
ccccccc
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=P-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFR4=(-P2*F(I)+Q*F(I+3))*P*P1*X+(P1*F(I+1)-P*F(I+2))*Q*P2*.5
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFR4=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFR4=F(N)
      RETURN
      end function
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!gauss-legendre integral,   N:number of intergal index， x1,x2:intergal range   
!x,w:gauss position rr and corresponding weight rrw  similar to dx(n)
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE gauleg(N,x1,x2,X,W)
        IMPLICIT NONE
        INTEGER N
        REAL*8 x1,x2,X(N),W(N)
        REAL*8 z1,z,xm,xl,pp,p3,p2,p1,pi,tol
        INTEGER m,i,j

        pi=acos(-1.0)
        tol=1.E-12

        m=(n+1)/2
        xm=0.5*(x2+x1)
        xl=0.5*(x2-x1)

         DO 10 i=1,m
         z=cos(pi*(i-0.25)/(N+0.5))

 20      CONTINUE
         p1=1.0E0
         p2=0.0E0
         DO 30 j=1,N
          p3=p2
          p2=p1
          p1=((2*j-1)*z*p2-(j-1)*p3)/j
 30      CONTINUE
         pp=N*(z*p1-p2)/(z*z-1.0E0)
         z1=z
         z=z1-p1/pp
         IF( abs(z1-z) .GT. tol) GOTO 20 ! Scheifenende

         X(i) = xm - xl*z
         X(n+1-i) = xm + xl*z
         W(i) = 2.E0*xl/((1.0-z*z)*pp*pp)
         W(n+1-i) = W(i)
 10     CONTINUE
        END SUBROUTINE gauleg
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
! five points derivative formula for second derivative
! y: function value array
! d2y: second derivative array
! n:size of the array
! uniform grid 
!!complex type 
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine second_derivative(y,d2y,n,dx)
       implicit none
       integer,intent(in)::n
       complex*16,intent(in)::dx
       complex*16,dimension(1:n),intent(in)::y
       complex*16,dimension(1:n),intent(out)::d2y
       integer::i 
       d2y(1)=(35.d0/12.d0*y(1)-26.d0/3.d0*y(2)+19.d0/2.d0*y(3)
     &   -14.d0/3.d0*y(4)+11.d0/12.d0*y(5))/(dx**2)
       d2y(2)=(11.d0/12.d0*y(1)-5.d0/3.d0*y(2)+1.d0/2.d0*y(3)
     &  +1.d0/3.d0*y(4)-1.d0/12.d0*y(5))/(dx**2)
       d2y(n-1)=(-1.d0/12.d0*y(N-4)+1.d0/3.d0*y(N-3)+1.d0/2.d0*y(N-2)
     &  -5.d0/3.d0*y(N-1)+11.d0/12.d0*y(N))/(dx**2)
       d2y(n)=(11.d0/12.d0*y(N-4)-14.d0/3.d0*y(N-3)+19.d0/2.d0*y(N-2)
     &  -26.d0/3.d0*y(N-1)+35.d0/12.d0*y(N))/(dx**2)
       do i=3,n-2
       d2y(i)=(-y(i-2)+16.d0*y(i-1)-30.d0*y(i)+ 
     & 16.d0*y(i+1)-y(i+2))/(12.d0*dx**2)
          !  write(*,*) d2y(i)
       end do
       end subroutine second_derivative
ccccccc
c *** Calculate d^2u(r)/dr^2 using five points derivative formula
c     f(ndim)=function to make derivative
c     h      =step
c     j      =point for derivative
      function deriv1(f,h,ndim,j)
        implicit none
        integer ndim,j
        real*8 f(ndim),h,deriv1

        if ((j.eq.1).or.(j.eq.2)) then
           deriv1=(-f(j+2)+4d0*f(j+1)-3d0*f(j))/2d0/h
        else if (j.eq.ndim-1) then
           deriv1=(3d0*f(j)-4d0*f(j-1)+f(j-2))/2d0/h
        else if (j.eq.ndim) then
           deriv1=0 !!!CHECK
        else ! five points formula
           deriv1=(f(j-2)-8*f(j-1)+8*f(j+1)-f(j+2))/h/12.
         end if
      end function deriv1
ccccccc
c-----------------------------------------------------------------------
!complex interpolation function for uniform grids
      FUNCTION FFC(PP,F,N)
      COMPLEX*16 FFC,F(N)
      REAL*8 PP
      PARAMETER(X=.16666666666667)
      I=PP
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=PP-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFC=(-P2*F(I)+Q*F(I+3))*(P*P1*X)+(P1*F(I+1)-P*F(I+2))*(Q*P2*.5)
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFC=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFC=F(N)
      RETURN
      END function
ccccccc
!Lagrange interpolation function for arbitrary grids (complex version)
      complex*16 function interpolateLagrange(xi, x, y, n)
                complex*16, intent(in) :: xi
                complex*16, intent(in) :: x(:)
                complex*16, intent(in) :: y(:)
                integer, intent(in) :: n
                complex*16 :: result
                complex*16 :: term
                integer :: i, j

                result = 0.0

                do i = 1, n
                    term = y(i)
                    do j = 1, n
                        if (i /= j) then
                            term = term * (xi - x(j)) / (x(i) - x(j))
                        end if
                    end do
                    result = result + term
                end do
                interpolateLagrange= result
            end function interpolateLagrange
ccccccc
!************************************************************************
!*     REAL 4-point lagrange interpolation routine.
!*     interpolates thr FUNCTION value fival at point r from an
!*     array of points stored in fdis(ndm). this array is assumed
!*     to be defined such that the first element fdis(1) CONTAINS
!*     the FUNCTION value at r=xv(1) and xv(2 .. ndm) are monotonically
!*     increasing.
!************************************************************************
      FUNCTION cfival(r,xv,fdis,ndm,alpha)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 cfival,fdis(ndm),y1,y2,y3,y4
      DIMENSION xv(ndm)
      IF(r.GT.xv(ndm)) go to 9
      DO 5 k=1,ndm-2
 5    IF(r.LT.xv(k)) go to 6
      k=ndm-2
 6    nst=MAX(k-1,1)
      x1=xv(nst)
      x2=xv(nst+1)
      x3=xv(nst+2)
      x4=xv(nst+3)
      y1=fdis(nst+0)
      y2=fdis(nst+1)
      y3=fdis(nst+2)
      y4=fdis(nst+3)
      pii1=(x1-x2)*(x1-x3)*(x1-x4)
      pii2=(x2-x1)*(x2-x3)*(x2-x4)
      pii3=(x3-x1)*(x3-x2)*(x3-x4)
      pii4=(x4-x1)*(x4-x2)*(x4-x3)
      xd1=r-x1
      xd2=r-x2
      xd3=r-x3
      xd4=r-x4
      pi1=xd2*xd3*xd4
      pi2=xd1*xd3*xd4
      pi3=xd1*xd2*xd4
      pi4=xd1*xd2*xd3
      cfival=y1*pi1/pii1+y2*pi2/pii2+y3*pi3/pii3+y4*pi4/pii4
      RETURN
 9    cfival=fdis(ndm) * EXP(alpha*(xv(ndm)-r))
      RETURN
      END function
      end module