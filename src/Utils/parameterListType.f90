!!*********************************************
!!
!! Author: Ezra Brooker
!! 2021
!! 
!! src/Utils/parameterListType.f90
!!
!! This module defines a parameter list derived
!! type object that one can use to pass an N
!! number of parameters to a function, provided
!! the function knows it will be receiving the
!! derived type list. This is a crude way to
!! achieve the **kwargs style function input of
!! Python functions.
!!
!! The primary procedures one needs to use with
!! an instance of this parameter list type are
!! the init, addEntry, getEntry procedures
!!
!!*********************************************
MODULE parameter_list_type
USE kindSettings
USE parameter_type
IMPLICIT NONE
PRIVATE

    PUBLIC :: param_list_t

    TYPE param_list_t
        PRIVATE
        INTEGER(ilkp)              :: length         !! length of list
        INTEGER(ilkp)              :: lastIndex = 0  !! last index filled
        TYPE(param_t), ALLOCATABLE :: params(:)      !! list of param_t types

    CONTAINS

        !! List initializer
        PROCEDURE :: init

        !! Add key-value entry routines
        PROCEDURE, PRIVATE :: addInt
        PROCEDURE, PRIVATE :: addReal
        PROCEDURE, PRIVATE :: addBool
        PROCEDURE, PRIVATE :: addChar
        GENERIC            :: addEntry => addInt,addReal,addBool,addChar

        !! Get key-value routines
        PROCEDURE, PRIVATE :: getInt
        PROCEDURE, PRIVATE :: getReal
        PROCEDURE, PRIVATE :: getBool
        PROCEDURE, PRIVATE :: getChar
        GENERIC            :: getEntry => getInt,getReal,getBool,getChar

        !! Some utility functions
        PROCEDURE          :: getLength
        PROCEDURE          :: getLastIndex
        PROCEDURE          :: printEntries
        PROCEDURE, PRIVATE :: isKeyInList

    END TYPE param_list_t


CONTAINS


    SUBROUTINE init(this,length)
        !! Initialize the list with list length
        CLASS(param_list_t)            :: this
        INTEGER(ilkp),      INTENT(IN) :: length

        this%length = length

        IF (ALLOCATED(this%params)) DEALLOCATE(this%params)
        ALLOCATE(this%params(this%length))

    END SUBROUTINE init

    FUNCTION getLength(this)
        !! Get the length of this list
        CLASS(param_list_t), INTENT(IN) :: this
        INTEGER(ilkp)                   :: getLength

        getLength = this%length

    END FUNCTION getLength

    FUNCTION getLastIndex(this)
        !! Get the last index of the list filled
        CLASS(param_list_t), INTENT(IN) :: this
        INTEGER(ilkp)                   :: getLastIndex

        getLastIndex = this%lastIndex

    END FUNCTION getLastIndex

    SUBROUTINE printEntries(this)
        !! Print all entries in the list
        CLASS(param_list_t), INTENT(IN) :: this
        INTEGER(ilkp)                   :: i
        CHARACTER(LEN=:), ALLOCATABLE   :: key

        DO i=1,this%lastIndex
      
            key = this%params(i)%getKey()
      
            SELECT TYPE(value => this%params(i)%getValue())

            TYPE IS(INTEGER(ilkp))
                WRITE(*,*) " key = ", key, ", value = ", value
      
            TYPE IS(REAL(rkp))
                WRITE(*,*) " key = ", key, ", value = ", value
      
            TYPE IS(CHARACTER(LEN=*))
                WRITE(*,*) " key = ", key, ", value = ", value
      
            TYPE IS(LOGICAL)
                WRITE(*,*) " key = ", key, ", value = ", value
      
            CLASS DEFAULT
                WRITE(*,*) " key = ", key, ", value = UNKOWN DATA TYPE"
      
            END SELECT

        END DO

    END SUBROUTINE printEntries


    SUBROUTINE isKeyInList(this,key,inlist,index)
        !! Check if key is in list and return index
        !! if it is
        CLASS(param_list_t), INTENT(IN)  :: this
        CHARACTER(LEN=*),    INTENT(IN)  :: key
        LOGICAL,             INTENT(OUT) :: inlist
        INTEGER(ilkp),       INTENT(OUT) :: index
        INTEGER(ilkp)                    :: i

        inlist = .false.
        index  = 0 
        DO i = 1, this%lastIndex
            IF (this%params(i)%getKey() == key) THEN
                inlist = .true.
                index = i
            END IF
        END DO

    END SUBROUTINE isKeyInList


    !!********************************!!
    !! Define type-dependent routines !!
    !! for adding key-value pairs     !!
    !!********************************!!

    SUBROUTINE addInt(this,key,value)
        !! Add integer parameter
        CLASS(param_list_t), INTENT(INOUT) :: this
        CHARACTER(LEN=*),    INTENT(IN)    :: key
        INTEGER(ilkp)                      :: value
        CLASS(*),            ALLOCATABLE   :: newval
        INTEGER(ilkp)                      :: ierr
        INTEGER(ilkp)                      :: index
        LOGICAL                            :: inlist

        ALLOCATE(newval, source=value)
        CALL this%isKeyInList(key,inlist,index) ! Check if key in list

        IF (inlist) THEN
            ! overwrite existing key-value entry
            CALL this%params(index)%addKeyVal(key,newval)
        
        ELSE
            ! add new key-value entry
            this%lastIndex = this%lastIndex + 1
            CALL this%params(this%lastIndex)%addKeyVal(key,newval)
        
        END IF

    END SUBROUTINE addInt

    SUBROUTINE addReal(this,key,value)
        !! Add real parameter
        CLASS(param_list_t), INTENT(INOUT) :: this
        CHARACTER(LEN=*),    INTENT(IN)    :: key
        REAL(rkp)                          :: value
        CLASS(*),            ALLOCATABLE   :: newval
        INTEGER(ilkp)                      :: ierr
        INTEGER(ilkp)                      :: index
        LOGICAL                            :: inlist

        ALLOCATE(newval, source=value)
        CALL this%isKeyInList(key,inlist,index)

        IF (inlist) THEN
        
            CALL this%params(index)%addKeyVal(key,newval)
        
        ELSE
            this%lastIndex = this%lastIndex + 1
            CALL this%params(this%lastIndex)%addKeyVal(key,newval)
        
        END IF

    END SUBROUTINE addReal

    SUBROUTINE addBool(this,key,value)
        !! Add boolean parameter
        CLASS(param_list_t), INTENT(INOUT) :: this
        CHARACTER(LEN=*),    INTENT(IN)    :: key
        LOGICAL                            :: value
        CLASS(*),            ALLOCATABLE   :: newval
        INTEGER(ilkp)                      :: ierr
        INTEGER(ilkp)                      :: index
        LOGICAL                            :: inlist

        ALLOCATE(newval, source=value)
        CALL this%isKeyInList(key,inlist,index)

        IF (inlist) THEN
        
            CALL this%params(index)%addKeyVal(key,newval)
        
        ELSE
            this%lastIndex = this%lastIndex + 1
            CALL this%params(this%lastIndex)%addKeyVal(key,newval)
        
        END IF

    END SUBROUTINE addBool

    SUBROUTINE addChar(this,key,value)
        !! Add string parameter
        CLASS(param_list_t), INTENT(INOUT) :: this
        CHARACTER(LEN=*),    INTENT(IN)    :: key
        CHARACTER(LEN=*),    INTENT(IN)    :: value
        CLASS(*),            ALLOCATABLE   :: newval
        INTEGER(ilkp)                      :: ierr
        INTEGER(ilkp)                      :: index
        LOGICAL                            :: inlist

        ALLOCATE(newval, source=value)
        CALL this%isKeyInList(key,inlist,index)

        IF (inlist) THEN
        
            CALL this%params(index)%addKeyVal(key,newval)
        
        ELSE
            this%lastIndex = this%lastIndex + 1
            CALL this%params(this%lastIndex)%addKeyVal(key,newval)
        
        END IF

    END SUBROUTINE addChar


    !!********************************!!
    !! Define type-dependent routines !!
    !! for fetching key-values        !!
    !!********************************!!

    SUBROUTINE getInt(this,key,value)
        CLASS(param_list_t), INTENT(IN)  :: this
        CHARACTER(LEN=*),    INTENT(IN)  :: key
        INTEGER(ilkp),       INTENT(OUT) :: value
        INTEGER(ilkp)                    :: ierr
        INTEGER(ilkp)                    :: index
        LOGICAL                          :: inlist

        CALL this%isKeyInList(key,inlist,index)
        IF (inlist) THEN
            SELECT TYPE (temp => this%params(index)%getValue())
            TYPE IS (INTEGER(ilkp))
                value = temp
            CLASS DEFAULT
                stop "[getEntry ERROR] Entry VALUE not type INTEGER"
            END SELECT
        ELSE
            stop "[getEntry ERROR] Entry KEY not in LIST"
        ENDIF
    END SUBROUTINE getInt


    SUBROUTINE getReal(this,key,value)
        CLASS(param_list_t), INTENT(IN)  :: this
        CHARACTER(LEN=*),    INTENT(IN)  :: key
        REAL(rkp),           INTENT(OUT) :: value
        INTEGER(ilkp)                    :: ierr
        INTEGER(ilkp)                    :: index
        LOGICAL                          :: inlist

        CALL this%isKeyInList(key,inlist,index)
        IF (inlist) THEN
            SELECT TYPE (temp => this%params(index)%getValue())
            TYPE IS (REAL(rkp))
                value = temp
            CLASS DEFAULT
                stop "[getEntry ERROR] Entry value not type REAL"
            END SELECT
        ELSE
            stop "[getEntry ERROR] Entry KEY not in LIST"
        ENDIF
    END SUBROUTINE getReal


    SUBROUTINE getChar(this,key,value)
        CLASS(param_list_t),           INTENT(IN)  :: this
        CHARACTER(LEN=*),              INTENT(IN)  :: key
        CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: value
        INTEGER(ilkp)                              :: ierr
        INTEGER(ilkp)                              :: index
        LOGICAL                                    :: inlist

        CALL this%isKeyInList(key,inlist,index)
        IF (inlist) THEN
            SELECT TYPE (temp => this%params(index)%getValue())
            TYPE IS (CHARACTER(LEN=*))
                value = temp
            CLASS DEFAULT
                stop "[getEntry ERROR] Entry value not type REAL"
            END SELECT
        ELSE
            stop "[getEntry ERROR] Entry KEY not in LIST"
        ENDIF
    END SUBROUTINE getChar


    SUBROUTINE getBool(this,key,value)
        CLASS(param_list_t), INTENT(IN)  :: this
        CHARACTER(LEN=*),    INTENT(IN)  :: key
        LOGICAL,             INTENT(OUT) :: value
        INTEGER(ilkp)                    :: ierr
        INTEGER(ilkp)                    :: index
        LOGICAL                          :: inlist

        CALL this%isKeyInList(key,inlist,index)
        IF (inlist) THEN
            SELECT TYPE (temp => this%params(index)%getValue())
            TYPE IS (LOGICAL)
                value = temp
            CLASS DEFAULT
                stop "[getEntry ERROR] Entry value not type LOGICAL"
            END SELECT
        ELSE
            stop "[getEntry ERROR] Entry KEY not in LIST"
        ENDIF
    END SUBROUTINE getBool

END MODULE parameter_list_type