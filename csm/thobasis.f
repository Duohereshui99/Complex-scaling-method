        module thobasis
            use parameter
            use hobasis
            contains
ccccccc       
      complex*16 function LSTFUN(gamma,m,r) !LST变换,相当于公式里面的s(r)
      implicit none
      real*8::gamma,m
      complex*16::r
      LSTFUN=(1/((1.d0/r)**m+
     & (1.d0/gamma/sqrt(r))**m))**(1.d0/m)
      end function
ccccccc
      complex*16 function D1LSTFUN(gamma,m,r)     !s(r)的一阶导数的解析表达式
      implicit none
      real*8::gamma,m
      complex*16::x,y       !计算的中间量
      complex*16::r
      x=(1.d0/r)**(m+1.d0)+1.d0/(2.0*r)*(1.d0/(gamma*sqrt(r)))**m
      y=(1.d0/r)**m+(1.d0/(gamma*sqrt(r)))**m
      D1LSTFUN=LSTFUN(gamma,m,r)*x/y
      end function
ccccccc
      complex*16 function THOFUNC(n,l,alpha,gamma,m,r)    !alpha是HOBASIS的无量纲参数
      implicit none
      integer::n,l
      real*8::alpha,gamma,m
      complex*16::r
      THOFUNC=sqrt(D1LSTFUN(gamma,m,r))
     & *ho3d(n,l,alpha,LSTFUN(gamma,m,r))
     & *LSTFUN(gamma,m,r)/r !再乘s(r)
      end function
        end module
