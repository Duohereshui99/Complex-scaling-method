ccccccc
        module eigenvalue
!complex eigenvalue solutions from lapack        
!n:dimension of matrix,a:matrix
!w:eigenvalues,vl,vr: left & right eigenvectors   
        contains
ccccccc
        subroutine ZGEEVS(n,a,w,vl,vr)             
            implicit none 
            integer::n
            integer::info,lda,ldvr,ldvl,lwork
            complex*16::a(n,n),w(n),vl(n,n),vr(n,n)
            complex*16::work(2*n),rwork(2*n)
        !    allocate(work(lwork),rwork(lwork))
            ldvr=n
            ldvl=n
            lda=n
            lwork=2*n
            call zgeev('N','V',n,A,lda,w,vl,ldvl,vr,ldvr,
     &       work,lwork,rwork,info)
           ! write(*,*) 'info=',info
           ! deallocate(work,rwork)
        end subroutine

        end module