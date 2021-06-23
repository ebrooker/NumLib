MODULE proctorOptimize 
USE kindSettings, ONLY : rkp
USE multiVarFuncIntrfc, ONLY : f_intrfc, gradf_intrfc
USE optimization
IMPLICIT NONE
PRIVATE


    PUBLIC :: proctorOpt

    CHARACTER(LEN=50) :: line="!!************************************************"
    CHARACTER(LEN=2 ) :: excl="!!"

    TYPE proctorOpt

        PROCEDURE(f_intrfc    ), POINTER, NOPASS :: f
        PROCEDURE(gradf_intrfc), POINTER, NOPASS :: gradf
        REAL(rkp)                                :: tol   = 1e-6
        REAL(rkp)                                :: norm  = HUGE(1.0_rkp)
        REAL(rkp)                                :: alph0 = 1.0
        REAL(rkp)                                :: alpha = HUGE(1.0_rkp)
        REAL(rkp)                                :: gamma = 0.5
        REAL(rkp)                                :: rho   = 0.5
        REAL(rkp)                                :: L2
        REAL(rkp), ALLOCATABLE                   :: x0(:,:)
        REAL(rkp), ALLOCATABLE                   :: x(:,:)
        REAL(rkp), ALLOCATABLE                   :: p(:,:)
        REAL(rkp), ALLOCATABLE                   :: B(:,:)
        REAL(rkp), ALLOCATABLE                   :: iB(:,:)

        INTEGER                                  :: m    = 2
        INTEGER                                  :: n    = 0
        INTEGER                                  :: nmax = 5000

        CONTAINS
            PROCEDURE :: init => constructor
            PROCEDURE :: test_lineSearch
            PROCEDURE :: test_BFGS
            PROCEDURE :: test_steepDescent
            PROCEDURE :: output
            PROCEDURE :: set_f
            PROCEDURE :: set_x0
            FINAL     :: destructor

    END TYPE proctorOpt

CONTAINS

    SUBROUTINE destructor(this)
        TYPE(proctorOpt) :: this

        WRITE(*,*) "[proctorOpt] Deallocating dynamic arrays..."
        IF (ALLOCATED(this%x0)) DEALLOCATE(this%x0)
        IF (ALLOCATED(this%x )) DEALLOCATE(this%x )
        IF (ALLOCATED(this%p )) DEALLOCATE(this%p )
        IF (ALLOCATED(this%B )) DEALLOCATE(this%B )
        IF (ALLOCATED(this%iB)) DEALLOCATE(this%iB)

        WRITE(*,*) "[proctorOpt] Nullifying procedure pointers..."
        IF (ASSOCIATED(this%f    )) this%   f  => NULL()
        IF (ASSOCIATED(this%gradf)) this%gradf => NULL()

        WRITE(*,*) "[proctorOpt] Class instance destroyed..."
    END SUBROUTINE destructor

    SUBROUTINE constructor(this,fin,gradfin,m,xi,B,iB)
        CLASS(proctorOpt)       :: this
        PROCEDURE(f_intrfc    ) :: fin
        PROCEDURE(gradf_intrfc) :: gradfin
        INTEGER,     INTENT(IN) :: m
        REAL(rkp),   INTENT(IN) :: xi(:,:),B(:,:),iB(:,:)
        INTEGER                 :: i

        CALL this%set_f(fin,gradfin)

        this%m = m
        IF (ALLOCATED(this%x0)) DEALLOCATE(this%x0)
        IF (ALLOCATED(this%x )) DEALLOCATE(this%x )
        IF (ALLOCATED(this%p )) DEALLOCATE(this%p )
        IF (ALLOCATED(this%B )) DEALLOCATE(this%B )
        IF (ALLOCATED(this%iB)) DEALLOCATE(this%iB)

        ALLOCATE( this%x0 (m,1), source=xi      )
        ALLOCATE( this%x  (m,1), source=xi      )
        ALLOCATE( this%p  (m,1), source=0.0_rkp )
        ALLOCATE( this%B  (m,m), source=B       )
        ALLOCATE( this%iB (m,m), source=iB      )

    END SUBROUTINE constructor


    SUBROUTINE set_f(this,fin,gradfin)
        CLASS(proctorOpt)       :: this
        PROCEDURE(f_intrfc    ) :: fin
        PROCEDURE(gradf_intrfc) :: gradfin

        IF (ASSOCIATED(this%f    )) this%f     => NULL()
        this%f  => fin
        IF (ASSOCIATED(this%gradf)) this%gradf => NULL()
        this%gradf => gradfin

    END SUBROUTINE set_f

    SUBROUTINE set_x0(this,xi)
        CLASS(proctorOpt)     :: this
        REAL(rkp), INTENT(IN) :: xi(:,:)

        this%x0 = xi
        this%x  = this%x0

    END SUBROUTINE set_x0

    SUBROUTINE test_lineSearch(this)
        CLASS(proctorOpt) :: this
        CALL this%test_steepDescent()
        CALL this%test_BFGS()

    END SUBROUTINE test_lineSearch

    SUBROUTINE test_steepDescent(this)
        CLASS(proctorOpt) :: this

        this%L2 = HUGE(1.0_rkp)
        this%n  = 0
        DO WHILE (this%L2 >= this%tol .and. this%n < this%nmax)
            
            CALL steepDescent(this%gradf,this%x,this%p)

            CALL backtrack(this%f,this%gradf,this%x,this%p,         &
                           this%alph0,this%rho,this%gamma,this%alpha)
            
            this%x0 = this%x
            this%x  = this%x + this%alpha*this%p
            this%n  = this%n + 1
            this%L2 = SQRT(SUM(ABS(this%x-this%x0)**2))
        
        ENDDO
        write(*,*) "Steepest Descent SOLUTION: "
        write(*,*) "n,x,L2 = ", this%n,this%x(:,1),this%L2
        
    END SUBROUTINE test_steepDescent


    SUBROUTINE test_BFGS(this)
        CLASS(proctorOpt) :: this
        LOGICAL           :: flag=.true.

        this%L2 = HUGE(1.0_rkp)
        this%n  = 0
        DO WHILE (this%L2 >= this%tol .and. this%n < this%nmax)
            
            CALL BFGS(this%gradf,flag,this%x,this%x0,    &
                                 this%p,this%B,this%iB)

            CALL backtrack(this%f,this%gradf,this%x,this%p,         &
                           this%alph0,this%rho,this%gamma,this%alpha)
            
            this%x0 = this%x
            this%x  = this%x + this%alpha*this%p
            this%n  = this%n + 1
            this%L2 = SQRT(SUM(ABS(this%x-this%x0)**2))

            if (this%n.gt.0) flag=.false.
        
        ENDDO
        write(*,*) "BFGS SOLUTION:"
        write(*,*) "n,x,L2 = ", this%n,this%x(:,1),this%L2
        
    END SUBROUTINE test_BFGS

    SUBROUTINE output(this)
        CLASS(proctorOpt) :: this


    END SUBROUTINE output

END MODULE proctorOptimize