C FILE: MBX_FORTRAN_TEST.F90
	  SUBROUTINE CALC_MAX9(ret, d)
	  implicit none
C
C     Calc background
C
      REAL*8 arr(9*2+(9-2)*2)											! g = gaussian result
	  REAL*8 ret
	  REAL*8 d(9,9)														! d = data = pixel values
Cf2py intent(in) d														! intent(in) means input
Cf2py intent(out) ret													! intent(out) means only this will be returned
Cf2py depend(d) ret
	  
      arr(1:9) = d(:,1)
	  arr(10:18) = d(:,9)
	  arr(19:25) = d(2:8,1)
	  arr(26:32) = d(2:8,9)
	  
	  ret = sum(arr)/32
	  
      END
	  
	  SUBROUTINE NORM5(ret, d)
	  implicit none
C
C     Norm for 5 variables
C
      REAL*8 d(5)											! g = gaussian result
	  REAL*8 ret													! d = data = pixel values
	  INTEGER i
Cf2py intent(in) d														! intent(in) means input
Cf2py intent(out) ret													! intent(out) means only this will be returned
Cf2py depend(d) ret
	  
      ret = 0
	  
	  do i = 1,5
		if (d(i) .gt. ret) then
		ret = d(i)
		endif
		if (-d(i) .gt. ret) then
		ret = -d(i)
		endif
	  enddo
	  
      END
	  
	  SUBROUTINE NORM5_2(ret, d)
	  implicit none
C
C     Norm for 5 variables
C
      REAL*8 d(5)											! g = gaussian result
	  REAL*8 ret													! d = data = pixel values
	  INTEGER i
Cf2py intent(in) d														! intent(in) means input
Cf2py intent(out) ret													! intent(out) means only this will be returned
Cf2py depend(d) ret
	  
      ret = 0
	  
	  do i = 1,5
		if (abs(d(i)) .gt. ret) then
		ret = abs(d(i))
		endif
	  enddo
	  
      END
	  
	  SUBROUTINE NORM5_3(ret, d)
	  implicit none
C
C     Norm for 5 variables
C
      REAL*8 d(5)											! g = gaussian result
	  REAL*8 ret													! d = data = pixel values
	  INTEGER i
Cf2py intent(in) d														! intent(in) means input
Cf2py intent(out) ret													! intent(out) means only this will be returned
Cf2py depend(d) ret
	  
      ret = maxval(abs(d))
	  
      END
	  
	  SUBROUTINE MAX9(ret, d)
	  implicit none
C
C     Norm for 5 variables
C
      REAL*8 d(9,9)											! g = gaussian result
	  REAL*8 ret													! d = data = pixel values
Cf2py intent(in) d														! intent(in) means input
Cf2py intent(out) ret													! intent(out) means only this will be returned
Cf2py depend(d) ret
	  	  
	  ret = maxval(d)
	  
      END
	  
	  SUBROUTINE MAX7(ret, d)
	  implicit none
C
C     Norm for 5 variables
C
      REAL*8 d(7,7)											! g = gaussian result
	  REAL*8 ret													! d = data = pixel values
Cf2py intent(in) d														! intent(in) means input
Cf2py intent(out) ret													! intent(out) means only this will be returned
Cf2py depend(d) ret
	  	  
	  ret = maxval(d)
	  
      END
	  	  
C END FILE MBX_FORTRAN_TEST.F90