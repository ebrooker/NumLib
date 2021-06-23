!!*********************************************
!!
!! Author: Ezra Brooker
!! 2021
!! 
!! src/Utils/parameterType.f90
!!
!! This module defines a parameter derived type
!! that is used to store a key and value pair
!! of any type that is defined or extended to.
!! 
!! For now, this supports CHARACTER, LOGICAL,
!! REAL, and INTEGER, but derived types can be
!! added for comprehension.
!!
!!*********************************************
MODULE parameter_type
IMPLICIT NONE
PRIVATE

    PUBLIC :: param_t


    TYPE param_t
        CHARACTER(LEN=:), ALLOCATABLE :: key
        CLASS(*),         ALLOCATABLE :: value

        CONTAINS
            PROCEDURE :: getKey
            PROCEDURE :: getValue
            PROCEDURE :: addKeyVal

    END TYPE param_t

CONTAINS

    FUNCTION getKey(this)
        CLASS(param_t),   INTENT(IN)  :: this
        CHARACTER(LEN=:), ALLOCATABLE :: getKey

        getKey = this%key

    END FUNCTION getKey

    FUNCTION getValue(this)
        CLASS(param_t), INTENT(IN)  :: this
        CLASS(*),       ALLOCATABLE :: getValue

        getValue = this%value

    END FUNCTION getValue

    SUBROUTINE addKeyVal(this,key,value)
        CLASS(param_t),    INTENT(INOUT) :: this
        CHARACTER(LEN=*),  INTENT(IN)    :: key
        CLASS(*),          INTENT(IN)    :: value

        this%key   = key
        this%value = value

    END SUBROUTINE addKeyVal

END MODULE parameter_type

