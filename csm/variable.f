ccccccc
        module variable
            implicit none
            real*8,allocatable::rr(:)    !real r coordinate first
            real*8,allocatable::rrw(:)  !integral weight
            complex*16,allocatable::r(:) !complex r coordinate (uniformed mesh)
            complex*16,allocatable::rc(:) !complex r coordinate
            complex*16,allocatable::vpot1(:)    
            complex*16,allocatable::vpot(:) !complex potential
            complex*16,allocatable::H(:,:)  !complex hamiltonian
            complex*16,allocatable::psi(:,:),d2psi(:,:)
            complex*16,allocatable::psi1(:),d2psi1(:)
            complex*16,allocatable::w(:)        !!complex eigenvalue
            complex*16,allocatable::vl(:,:),vr(:,:) !!complex eigenvector (left&right)
            complex*16,allocatable::wf(:)           !!!test wf
        end module