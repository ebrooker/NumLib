!!
!!
!!
MODULE proctorClass
USE kindSettings, ONLY : rkp
USE functionInterface, ONLY : f_1d_intrfc
USE differentiation
IMPLICIT NONE
PRIVATE

    PUBLIC :: proctor_1d

    CHARACTER(LEN=50) :: line="!!************************************************"
    CHARACTER(LEN=2 ) :: excl="!!"

    TYPE proctor_1d

        PROCEDURE(f_1d_intrfc), POINTER, NOPASS :: f
        PROCEDURE(f_1d_intrfc), POINTER, NOPASS :: df
        REAL(rkp)                               :: x = 0.0
        REAL(rkp)                               :: h = 0.1
        REAL(rkp)                               :: y = 0.0
        REAL(rkp)                               :: deriv(3,3) = 0.0
        REAL(rkp)                               :: error(3,3) = 0.0
        INTEGER                                 :: ideriv
        CHARACTER(LEN=3)                        :: cderiv

        CONTAINS
            PROCEDURE :: init => constructor
            PROCEDURE :: dtest
            PROCEDURE :: compute
            PROCEDURE :: output
            PROCEDURE :: rextest => richExtrapTest
            PROCEDURE :: set_f
            FINAL     :: destructor

    END TYPE proctor_1d

CONTAINS

    SUBROUTINE destructor(this)
        TYPE(proctor_1d) :: this

        ! WRITE(*,*) "[PROCTOR_1D] Deallocating dynamic arrays..."
        ! IF (ALLOCATED(this%x    )) DEALLOCATE(this%x    )
        ! IF (ALLOCATED(this%h    )) DEALLOCATE(this%h    )
        ! IF (ALLOCATED(this%deriv)) DEALLOCATE(this%deriv)
        ! IF (ALLOCATED(this%error)) DEALLOCATE(this%error)

        WRITE(*,*) "[PROCTOR_1D] Nullifying procedure pointers..."
        IF (ASSOCIATED(this%f )) this%f  => NULL()
        IF (ASSOCIATED(this%df)) this%df => NULL()

        WRITE(*,*) "[PROCTOR_1D] Class instance destroyed..."
    END SUBROUTINE destructor

    SUBROUTINE constructor(this,fin,dfin,ideriv,xin,hin)
        CLASS(proctor_1d)               :: this
        PROCEDURE(f_1d_intrfc)          :: fin
        PROCEDURE(f_1d_intrfc)          :: dfin
        REAL(rkp), OPTIONAL, INTENT(IN) :: xin
        REAL(rkp), OPTIONAL, INTENT(IN) :: hin
        INTEGER,             INTENT(IN) :: ideriv

        CALL this%set_f(fin,dfin,ideriv)
        IF (PRESENT(xin)) this%x = xin
        IF (PRESENT(hin)) this%h = hin

    END SUBROUTINE constructor


    SUBROUTINE set_f(this,fin,dfin,ideriv)
        CLASS(proctor_1d)      :: this
        PROCEDURE(f_1d_intrfc) :: fin
        PROCEDURE(f_1d_intrfc) :: dfin
        INTEGER, INTENT(IN)    :: ideriv

        IF (ASSOCIATED(this%f )) this%f  => NULL()
        this%f  => fin
        IF (ASSOCIATED(this%df)) this%df => NULL()
        this%df => dfin

        SELECT CASE (ideriv)
            CASE (1)
                this%cderiv="1st"
            CASE (2)
                this%cderiv="2nd"
            CASE (3)
                this%cderiv="3rd"
            CASE(4)
                this%cderiv="4th"
            CASE default
                this%cderiv="1st"
        END SELECT
        this%ideriv = ideriv

    END SUBROUTINE set_f


    SUBROUTINE dtest(this)
        CLASS(proctor_1d) :: this
        CALL this%compute()
        CALL this%output()
        CALL this%rextest()

    END SUBROUTINE dtest


    SUBROUTINE compute(this)
        CLASS(proctor_1d) :: this
        this%deriv(1,1) = cntrDiff(this%f,this%x,this%h,deriv=this%ideriv        )
        this%deriv(1,2) = cntrDiff(this%f,this%x,this%h,deriv=this%ideriv,order=2)
        this%deriv(1,3) = cntrDiff(this%f,this%x,this%h,deriv=this%ideriv,order=4)
        this%deriv(2,1) = fwrdDiff(this%f,this%x,this%h,deriv=this%ideriv        )
        this%deriv(2,2) = fwrdDiff(this%f,this%x,this%h,deriv=this%ideriv,order=1)
        this%deriv(2,3) = fwrdDiff(this%f,this%x,this%h,deriv=this%ideriv,order=2)
        this%deriv(3,1) = bwrdDiff(this%f,this%x,this%h,deriv=this%ideriv        )
        this%deriv(3,2) = bwrdDiff(this%f,this%x,this%h,deriv=this%ideriv,order=1)
        this%deriv(3,3) = bwrdDiff(this%f,this%x,this%h,deriv=this%ideriv,order=2)
    END SUBROUTINE compute


    SUBROUTINE output(this)
        CLASS(proctor_1d) :: this
        write(*,*) ""
        write(*,*) line
        write(*,*) excl
        write(*,*) "!! Computing ",this%cderiv," derivative"
        write(*,*) excl
        write(*,*) "!! cntrDiff default =", this%deriv(1,1)
        write(*,*) "!! cntrDiff 2nd ord =", this%deriv(1,2)
        write(*,*) "!! cntrDiff 4th ord =", this%deriv(1,3)
        write(*,*) "!! ",this%cderiv," derivative   =", this%df(this%x)
        write(*,*) excl
        write(*,*) "!! fwrdDiff default =", this%deriv(2,1)
        write(*,*) "!! fwrdDiff 1st ord =", this%deriv(2,2)
        write(*,*) "!! fwrdDiff 2nd ord =", this%deriv(2,3)
        write(*,*) "!! ",this%cderiv," derivative   =", this%df(this%x)
        write(*,*) excl
        write(*,*) "!! bwrdDiff default =", this%deriv(3,1)
        write(*,*) "!! bwrdDiff 1st ord =", this%deriv(3,2)
        write(*,*) "!! bwrdDiff 2nd ord =", this%deriv(3,3)
        write(*,*) "!! ",this%cderiv," derivative   =", this%df(this%x)
        write(*,*) excl
        write(*,*) line
        write(*,*) ""

    END SUBROUTINE output


    SUBROUTINE richExtrapTest(this)
        CLASS(proctor_1d) :: this
        REAL(rkp)         :: yh(3),yh2(3)

        this%deriv(1,1) = cntrDiff(this%f,this%x,0.5*this%h,deriv=this%ideriv,order=2)
        this%deriv(2,1) = fwrdDiff(this%f,this%x,0.5*this%h,deriv=this%ideriv,order=1)        
        this%deriv(3,1) = bwrdDiff(this%f,this%x,0.5*this%h,deriv=this%ideriv,order=1)        

        write(*,*) ""
        write(*,*) line
        write(*,*) excl
        write(*,*) "!! Computing ",this%cderiv," derivative"
        write(*,*) excl
        write(*,*) "!! cntrDiff 2nd ord h/2 =", this%deriv(1,1)
        write(*,*) "!! cntrDiff 2nd ord h   =", this%deriv(1,2)
        write(*,*) "!! cntrDiff 4th ord h   =", this%deriv(1,3)
        write(*,*) "!! cntrDiff richExtrap  =", richExtrap(this%deriv(1,2),this%deriv(1,1))
        write(*,*) "!! ",this%cderiv," derivative       =", this%df(this%x)
        write(*,*) excl
        write(*,*) "!! fwrdDiff 1st ord h/2 =", this%deriv(2,1)
        write(*,*) "!! fwrdDiff 1st ord h   =", this%deriv(2,2)
        write(*,*) "!! fwrdDiff 2nd ord h   =", this%deriv(2,3)
        write(*,*) "!! fwrdDiff richExtrap  =", richExtrap(this%deriv(2,2),this%deriv(2,1))
        write(*,*) "!! ",this%cderiv," derivative       =", this%df(this%x)
        write(*,*) excl
        write(*,*) "!! bwrdDiff 1st ord h/2 =", this%deriv(3,1)
        write(*,*) "!! bwrdDiff 1st ord h   =", this%deriv(3,2)
        write(*,*) "!! bwrdDiff 2nd ord h   =", this%deriv(3,3)
        write(*,*) "!! bwrdDiff richExtrap  =", richExtrap(this%deriv(3,2),this%deriv(3,1))
        write(*,*) "!! ",this%cderiv," derivative       =", this%df(this%x)
        write(*,*) excl
        write(*,*) line
        write(*,*) ""

    END SUBROUTINE richExtrapTest


END MODULE proctorClass