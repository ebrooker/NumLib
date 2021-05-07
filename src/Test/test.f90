PROGRAM test
USE kindSettings, ONLY : rkp, ilkp, iskp
USE numericalDifferentiation
USE funcs
USE derivatives
IMPLICIT NONE

	REAL(rkp) :: x, h, ext(2), y1(3), y2(3)

	h = 0.01
	x = 2.0
	ext = 1.0

    write(*,*) ""
    write(*,*) "!!************************************************"
    write(*,*) "!!"
    write(*,*) "!! 1st derivative test #1"
    write(*,*) "!!"
    write(*,*) "!! cntrDiff default =", cntrDiff_d1(f1,x,h        )
    write(*,*) "!! cntrDiff 2nd ord =", cntrDiff_d1(f1,x,h,order=2)
    write(*,*) "!! cntrDiff 4th ord =", cntrDiff_d1(f1,x,h,order=4)
    write(*,*) "!! 1st derivative   =", df1(x)
    write(*,*) "!!"
    write(*,*) "!! fwrdDiff default =", fwrdDiff_d1(f1,x,h        )
    write(*,*) "!! fwrdDiff 1st ord =", fwrdDiff_d1(f1,x,h,order=1)
    write(*,*) "!! fwrdDiff 2nd ord =", fwrdDiff_d1(f1,x,h,order=2)
    write(*,*) "!! 1st derivative   =", df1(x)
    write(*,*) "!!"
    write(*,*) "!! bwrdDiff default =", bwrdDiff_d1(f1,x,h        )
    write(*,*) "!! bwrdDiff 1st ord =", bwrdDiff_d1(f1,x,h,order=1)
    write(*,*) "!! bwrdDiff 2nd ord =", bwrdDiff_d1(f1,x,h,order=2)
    write(*,*) "!! 1st derivative   =", df1(x)
    write(*,*) "!!"
    write(*,*) "!!************************************************"
    write(*,*) ""

    write(*,*) ""
    write(*,*) "!!************************************************"
    write(*,*) "!!"
    write(*,*) "!! 2nd derivative test #1"
    write(*,*) "!!"
    write(*,*) "!! cntrDiff default =", cntrDiff_d2(f1,x,h        )
    write(*,*) "!! cntrDiff 2nd ord =", cntrDiff_d2(f1,x,h,order=2)
    write(*,*) "!! cntrDiff 4th ord =", cntrDiff_d2(f1,x,h,order=4)
    write(*,*) "!! 2nd derivative   =", d2f1(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!! fwrdDiff default =", fwrdDiff_d2(f1,x,h        )
    write(*,*) "!! fwrdDiff 1st ord =", fwrdDiff_d2(f1,x,h,order=1)
    write(*,*) "!! fwrdDiff 2nd ord =", fwrdDiff_d2(f1,x,h,order=2)
    write(*,*) "!! 2nd derivative   =", d2f1(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!! bwrdDiff default =", bwrdDiff_d2(f1,x,h        )
    write(*,*) "!! bwrdDiff 1st ord =", bwrdDiff_d2(f1,x,h,order=1)
    write(*,*) "!! bwrdDiff 2nd ord =", bwrdDiff_d2(f1,x,h,order=2)
    write(*,*) "!! 2nd derivative   =", d2f1(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!!************************************************"
    write(*,*) ""

    write(*,*) ""
    write(*,*) "!!************************************************"
    write(*,*) "!!"
    write(*,*) "!! 1st derivative test #2"
    write(*,*) "!!"
    write(*,*) "!! cntrDiff default =", cntrDiff_d1(f2,x,h,        args=ext)
    write(*,*) "!! cntrDiff 2nd ord =", cntrDiff_d1(f2,x,h,order=2,args=ext)
    write(*,*) "!! cntrDiff 4th ord =", cntrDiff_d1(f2,x,h,order=4,args=ext)
    write(*,*) "!! 1st derivative   =", df2(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!! fwrdDiff default =", fwrdDiff_d1(f2,x,h,        args=ext)
    write(*,*) "!! fwrdDiff 1st ord =", fwrdDiff_d1(f2,x,h,order=1,args=ext)
    write(*,*) "!! fwrdDiff 2nd ord =", fwrdDiff_d1(f2,x,h,order=2,args=ext)
    write(*,*) "!! 1st derivative   =", df2(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!! bwrdDiff default =", bwrdDiff_d1(f2,x,h,        args=ext)
    write(*,*) "!! bwrdDiff 1st ord =", bwrdDiff_d1(f2,x,h,order=1,args=ext)
    write(*,*) "!! bwrdDiff 2nd ord =", bwrdDiff_d1(f2,x,h,order=2,args=ext)
    write(*,*) "!! 1st derivative   =", df2(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!!************************************************"
    write(*,*) ""


    write(*,*) ""
    write(*,*) "!!************************************************"
    write(*,*) "!!"
    write(*,*) "!! 2nd derivative test #2"
    write(*,*) "!!"
    write(*,*) "!! cntrDiff default =", cntrDiff_d2(f2,x,h,        args=ext)
    write(*,*) "!! cntrDiff 2nd ord =", cntrDiff_d2(f2,x,h,order=2,args=ext)
    write(*,*) "!! cntrDiff 4th ord =", cntrDiff_d2(f2,x,h,order=4,args=ext)
    write(*,*) "!! 2nd derivative   =", d2f2(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!! fwrdDiff default =", fwrdDiff_d2(f2,x,h,        args=ext)
    write(*,*) "!! fwrdDiff 1st ord =", fwrdDiff_d2(f2,x,h,order=1,args=ext)
    write(*,*) "!! fwrdDiff 2nd ord =", fwrdDiff_d2(f2,x,h,order=2,args=ext)
    write(*,*) "!! 2nd derivative   =", d2f2(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!! bwrdDiff default =", bwrdDiff_d2(f2,x,h,        args=ext)
    write(*,*) "!! bwrdDiff 1st ord =", bwrdDiff_d2(f2,x,h,order=1,args=ext)
    write(*,*) "!! bwrdDiff 2nd ord =", bwrdDiff_d2(f2,x,h,order=2,args=ext)
    write(*,*) "!! 2nd derivative   =", d2f2(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!!************************************************"
    write(*,*) ""


    write(*,*) ""
    write(*,*) "!!************************************************"
    write(*,*) "!!"
    write(*,*) "!! 3rd derivative test"
    write(*,*) "!!"
    write(*,*) "!! cntrDiff default =", cntrDiff_d3(f2,x,h,        args=ext)
    write(*,*) "!! cntrDiff 2nd ord =", cntrDiff_d3(f2,x,h,order=2,args=ext)
    write(*,*) "!! cntrDiff 4th ord =", cntrDiff_d3(f2,x,h,order=4,args=ext)
    write(*,*) "!! 3rd derivative   =", d3f2(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!! fwrdDiff default =", fwrdDiff_d3(f2,x,h,        args=ext)
    write(*,*) "!! fwrdDiff 1st ord =", fwrdDiff_d3(f2,x,h,order=1,args=ext)
    write(*,*) "!! fwrdDiff 2nd ord =", fwrdDiff_d3(f2,x,h,order=2,args=ext)
    write(*,*) "!! 3rd derivative   =", d3f2(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!! bwrdDiff default =", bwrdDiff_d3(f2,x,h,        args=ext)
    write(*,*) "!! bwrdDiff 1st ord =", bwrdDiff_d3(f2,x,h,order=1,args=ext)
    write(*,*) "!! bwrdDiff 2nd ord =", bwrdDiff_d3(f2,x,h,order=2,args=ext)
    write(*,*) "!! 3rd derivative   =", d3f2(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!!************************************************"
    write(*,*) ""


    write(*,*) ""
    write(*,*) "!!************************************************"
    write(*,*) "!!"
    write(*,*) "!! 4th derivative test"
    write(*,*) "!!"
    write(*,*) "!! cntrDiff default =", cntrDiff_d4(f2,x,h,        args=ext)
    write(*,*) "!! cntrDiff 2nd ord =", cntrDiff_d4(f2,x,h,order=2,args=ext)
    write(*,*) "!! cntrDiff 4th ord =", cntrDiff_d4(f2,x,h,order=4,args=ext)
    write(*,*) "!! 4th derivative   =", d4f2(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!! fwrdDiff default =", fwrdDiff_d4(f2,x,h,        args=ext)
    write(*,*) "!! fwrdDiff 1st ord =", fwrdDiff_d4(f2,x,h,order=1,args=ext)
    write(*,*) "!! fwrdDiff 2nd ord =", fwrdDiff_d4(f2,x,h,order=2,args=ext)
    write(*,*) "!! 4th derivative   =", d4f2(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!! bwrdDiff default =", bwrdDiff_d4(f2,x,h,        args=ext)
    write(*,*) "!! bwrdDiff 1st ord =", bwrdDiff_d4(f2,x,h,order=1,args=ext)
    write(*,*) "!! bwrdDiff 2nd ord =", bwrdDiff_d4(f2,x,h,order=2,args=ext)
    write(*,*) "!! 4th derivative   =", d4f2(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!!************************************************"
    write(*,*) ""

    y1(1) = cntrDiff_d4(f2,x,    h,order=2,args=ext)
    y2(1) = cntrDiff_d4(f2,x,0.5*h,order=2,args=ext)

    y1(2) = fwrdDiff_d4(f2,x,    h,order=1,args=ext)
    y2(2) = fwrdDiff_d4(f2,x,0.5*h,order=1,args=ext)
    
    y1(3) = bwrdDiff_d4(f2,x,    h,order=1,args=ext)
    y2(3) = bwrdDiff_d4(f2,x,0.5*h,order=1,args=ext)

    write(*,*) ""
    write(*,*) "!!************************************************"
    write(*,*) "!!"
    write(*,*) "!! 4th derivative test w/ Richardson Extrapolation"
    write(*,*) "!!"
    write(*,*) "!! cntrDiff 2nd ord    =", y1(1)
    write(*,*) "!! cntrDiff 2nd ord    =", y2(1)
    write(*,*) "!! cntrDiff 4th ord    =", cntrDiff_d4(f2,x,    h,order=4,args=ext)
    write(*,*) "!! cntrDiff richExtrap =", richExtrap(y1(1),y2(1))
    write(*,*) "!! 4th derivative      =", d4f2(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!! fwrdDiff 1st ord    =", y1(2) 
    write(*,*) "!! fwrdDiff 1st ord    =", y2(2) 
    write(*,*) "!! fwrdDiff 2nd ord    =", fwrdDiff_d4(f2,x,    h,order=2,args=ext)
    write(*,*) "!! fwrdDiff richExtrap =", richExtrap(y1(2),y2(2))
    write(*,*) "!! 4th derivative      =", d4f2(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!! bwrdDiff 1st ord    =", y1(3)
    write(*,*) "!! bwrdDiff 1st ord    =", y2(3)
    write(*,*) "!! bwrdDiff 2nd ord    =", bwrdDiff_d4(f2,x,    h,order=2,args=ext)
    write(*,*) "!! bwrdDiff richExtrap =", richExtrap(y1(3),y2(3))
    write(*,*) "!! 4th derivative      =", d4f2(x,args=ext)
    write(*,*) "!!"
    write(*,*) "!!************************************************"
    write(*,*) ""


END PROGRAM test