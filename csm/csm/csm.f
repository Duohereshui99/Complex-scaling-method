        program main
            use parameter
            use hobasis
            use thobasis
            use potentialfunctions
            use algorithm
            use mesh
            use system 
            use lstpar
            use potential
            use variable
            use complexscaling
            use eigenvalue
ccccccc
            implicit none
ccccccc
            integer::i,j,k
            real*8::t1,t2
            complex*16::s
ccccccc
            namelist /systems/ z1,z2,mass_1,mass_2,L
            namelist /meshs/ n_int,n_diff,hcm
            namelist /lst/ gamma,m,b,n_basis
            namelist /pots/ v0,r0,a0
            namelist /cplxscaling/ theta
ccccccc
            call get_info()
            call cpu_time(t1)
ccccccc
            read(5,nml=systems)
            read(5,nml=lst)
            read(5,nml=pots)
            read(5,nml=meshs)
            read(5,nml=cplxscaling)
ccccccc
            allocate(r(1:n_diff))               !!complex uniformed coordinate r
            allocate(rr(1:n_int))          
            allocate(rrw(1:n_int))      
            allocate(rc(1:n_int))               !!complex gauss coordinate rc
ccccccc
            allocate(vpot1(1:n_diff))
            allocate(vpot(1:n_int))             !!complex v
            allocate(H(0:n_basis,0:n_basis))    !!complex H
            allocate(w(0:n_basis))              !!complex eigenvalue
            allocate(vl(0:n_basis,0:n_basis))   !!complex eigenvector
            allocate(vr(0:n_basis,0:n_basis))   !!complex eigenvector
ccccccc
            allocate(psi1(1:n_diff),d2psi1(1:n_diff)) 
            allocate(psi(1:n_int,0:n_basis),d2psi(1:n_int,0:n_basis))
            allocate(psi_1(1:n_diff,0:n_basis))
            allocate(wf(1:n_diff,0:n_basis))
!! n column
ccccccc
            allocate(wavefunction(1:n_diff))
            allocate(d2wavefunction(1:n_diff))
ccccccc
            mu=amu*mass_1*mass_2/(mass_1+mass_2)
            z12=z1*z2
            alpha=1d0/2d0/b**2
ccccccc
 77         format(A,F15.7)
 777        format(A,I3)
            write(*,77)  'omega=',hbarc/mu/b**2
            write(*,777) 'nbasis=',n_basis
ccccccc
            call gauleg(n_int,0d0,hcm*n_diff,rr,rrw)     !!real rr (gauss points)
ccccccc
            rc=rr*exp(ii*theta*pi/180d0)                 !!rotated gauss points rc
ccccccc
            do i=1,n_diff
                r(i)=i*hcm                      
            end do
ccccccc
            r=r*exp(ii*theta*pi/180d0)                  !!rotated r (uniformed mesh) 
ccccccc
            do i=1,n_diff
                vpot1(i)=gausspot(r(i),v0,r0,a0)!potencc(r(i))!gausspot(r(i),v0,r0,a0)        !!rotated V
            end do
ccccccc
            do i=1,n_int
                vpot(i)=FFC(rr(i)/hcm,vpot1,n_diff)
            end do
ccccccc
        do i=0,n_basis
            do j=1,n_diff           !!real basis, not rotated!
                s=j*hcm   
                psi1(j)=THOFUNC(i,L,alpha,gamma,m,s)*s
            end do
            psi_1(:,i)=psi1
            call second_derivative(psi1,d2psi1,n_diff,hcm*exp(ii*0d0))
            do j=1,n_int
                psi(j,i)=FFC(rr(j)/hcm,psi1,n_diff)
                d2psi(j,i)=FFC(rr(j)/hcm,d2psi1,n_diff)
            end do
        end do
ccccccc
        do i=0,n_basis              !complex H matrix
            do j=0,n_basis
                do k=1,n_int
                    H(i,j)=H(i,j)+psi(k,i)*vpot(k)*psi(k,j)*rrw(k)      !!rotated V
     &      +exp(-2*ii*theta*pi/180)*(-hbarc**2/2/mu)*rrw(k)*d2psi(k,i)
     &        *psi(k,j)
     &      +exp(-2*ii*theta*pi/180)*(hbarc**2/2/mu)*(l+1d0)*l/rr(k)**2
     &        *rrw(k)*psi(k,i)*psi(k,j)
                end do
            end do 
        end do
ccccccc
         do i=0,n_basis
            write(20,*) H(i,:)
         end do
ccccccc
         call ZGEEVS(n_basis+1,H,w,vl,vr)
ccccccc
         do i=0,n_basis     !!original complex eigenvalues with index correlated to eigenvectors
            write(21,*) w(i)
         end do
ccccccc
         do i=0,n_basis     !!ith eigenvector,jth integral point
                do j=1,n_diff
                    do k=0,n_basis
                        wf(j,i)=wf(j,i)+psi_1(j,k)*vr(k,i)
                    end do
                end do
         end do
ccccccc
!! multiply by complex scaling factor exp(i\theta/2)
         wf=exp(ii*theta*pi/360d0)*wf
ccccccc
!!norm of wf**2, square norm of wf
         do i=0,n_basis
            s=0
            do j=1,n_diff
                s=s+wf(j,i)**2*hcm
            end do
            wf(:,i)=wf(:,i)/sqrt(s)
         end do
ccccccc
         do i=0,n_basis
            do j=1,n_diff
                write(22,*) abs(r(j)),real(wf(j,i))
            end do
            write(22,*) '&'
         end do
ccccccc
       !  k=14                      !!uniformed mesh for kth's eigenwf's expectation value
!         do k=0,n_basis
!          wavefunction(:)=wf(:,k)        
!          call second_derivative(wavefunction,d2wavefunction,n_diff,hcm*exp(ii*0d0))
!          s=0
! ccccccc
!             do i=1,n_diff
!                 s=s+wavefunction(i)**2*vpot1(i)*hcm      !!rotated V
!      &      +exp(-2*ii*theta*pi/180)*(-hbarc**2/2/mu)*hcm*d2wavefunction(i)
!      &        *wavefunction(i)
!      &      +exp(-2*ii*theta*pi/180)*(hbarc**2/2/mu)*(l+1d0)*l/abs(r(i))**2
!      &        *hcm*wavefunction(i)**2
!             end do
!             write(*,*) k,'expectation value=',s
!         end do
!             do i=1,n_diff
!                 write(23,*) abs(r(i)),real(wavefunction(i))
!             end do
! !!bubbling sort
         do i=0,n_basis
            do j=i+1,n_basis
                if(real(w(j))<real(w(i))) then
                    s=w(i)                          !!sort eigenvalues
                    w(i)=w(j)
                    w(j)=s
                end if
            end do
         end do
ccccccc
 200      format('***********complex eigenvalues***********')
 201      format('=========================================')
 202      format(I3,3X,A,F15.7,1X,A,1X,F15.7,A,3X,A)
 203      format(A,F7.3,1X,A)
ccccccc
         write(*,203) 'theta=',theta,'degree'
         write(*,200)
         do i=0,n_basis                             !!actually theta*0.9>-(1/2)arg(E)
            if(theta*pi/180d0*0.9d0>-atan2(aimag(w(i)),real(w(i)))/2d0.and.aimag(w(i))<0.and.real(w(i)).gt.0)  then 
                write(*,202) i,'Er=',real(w(i)),'+',aimag(w(i)) ,'i','resonance'
 !    &            'Gamma=',-2*aimag(w(i)),'MeV',
 !    &            't_half=',hbarc*log(2d0)/(-2*aimag(w(i)))*ratio,'s'
            else 
                write(*,202) i,'E=',real(w(i)),'+',aimag(w(i)),'i'
            end if
         end do
ccccccc 
         write(*,201)
ccccccc

            deallocate(r,rc,rr,rrw)
            deallocate(vpot,vpot1)
            deallocate(H,w,vl,vr)
            deallocate(psi,psi1,d2psi,d2psi1)
            deallocate(wf,psi_1)
            deallocate(wavefunction,d2wavefunction)
ccccccc            
            call cpu_time(t2)
            write(*,*) 'running time=',t2-t1
ccccccc
         contains
ccccccc
            subroutine get_info()
            ! 使用预处理器检查宏是否被定义
#ifdef BASE
        print *, 'Base directory: ', BASE
#endif

#ifdef VERDATE
        print *, 'Version date: ', VERDATE
#endif

#ifdef VERREV
        print *, 'Version revision: ', VERREV
#endif

#ifdef COMPDATE
        print *, 'Compilation date: ', COMPDATE
#endif
            end subroutine

        end program