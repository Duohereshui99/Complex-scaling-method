ccccccc
        module mesh 
            implicit none
            integer::n_int,n_diff   !integral and differential mesh
            real*8::hcm             !mesh size, n_diff*hcm=range
        end module
ccccccc
        module system
            implicit none
            real*8::z1,z2           !charge num 
            real*8::z12
            real*8::mass_1,mass_2   !mass num
            real*8::mu              !reduced mass
            integer::L              !angular momentum num
        end module
ccccccc
        module lstpar
            implicit none
            real*8::gamma,m         !lst parameter
            real*8::b,alpha         
            integer::n_basis        !real basis num=n_basis+1
            !n_basis ranges from 0 to ...
        end module
ccccccc
        module potential
            implicit none
            real*8::v0 
            real*8::r0 
            real*8::a0
        end module
ccccccc
        module complexscaling
            implicit none
            real*8::theta  !!complex scaling angle (degree)
        end module
ccccccc