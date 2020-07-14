C FILE: MBX_FORTRAN_TEST.F90
	  SUBROUTINE CALC_BG9(ret, d)
	  implicit none
C
C     Calc background
C
      REAL*8 arr(9*2+(9-2)*2)											! arr = temporary array
	  REAL*8 ret														! ret = return
	  REAL*8 d(9,9)														! d = data = pixel values
Cf2py intent(in) d														! intent(in) means input
Cf2py intent(out) ret													! intent(out) means only this will be returned
Cf2py depend(d) ret
	  
      arr(1:9) = d(:,1)
	  arr(10:18) = d(:,9)
	  arr(19:25) = d(2:8,1)
	  arr(26:32) = d(2:8,9)												! copy side of array
	  
	  ret = sum(arr)/32													! calculate mean
	  
      END
	  
	  SUBROUTINE CALC_BG7(ret, d)
	  implicit none
C
C     Calc background
C
      REAL*8 arr(7*2+(7-2)*2)											! arr = temporary array
	  REAL*8 ret														! ret = return
	  REAL*8 d(7,7)														! d = data = pixel values
Cf2py intent(in) d														! intent(in) means input
Cf2py intent(out) ret													! intent(out) means only this will be returned
Cf2py depend(d) ret
	  
      arr(1:7) = d(:,1)
	  arr(8:14) = d(:,7)
	  arr(15:19) = d(2:6,1)
	  arr(20:24) = d(2:6,7)												! copy side of array
	  
	  ret = sum(arr)/24													! calculate mean
	  
      END
	  
	  SUBROUTINE NORM5(ret, d)
	  implicit none
C
C     Norm for 5 variables
C
      REAL*8 d(5)														! d = data = 5 variables
	  REAL*8 ret														! ret = return
	  INTEGER i
Cf2py intent(in) d														! intent(in) means input
Cf2py intent(out) ret													! intent(out) means only this will be returned
Cf2py depend(d) ret
	  
      ret = 0
	  
	  do i = 1,5														! loop over all and find absolute largest value
		if (abs(d(i)) .gt. ret) then
		ret = abs(d(i))
		endif
	  enddo
	  
      END
	  
	  SUBROUTINE NORM6(ret, d)
	  implicit none
C
C     Norm for 6 variables
C
      REAL*8 d(6)														! d = data = 6 variables
	  REAL*8 ret														! ret = return
	  INTEGER i
Cf2py intent(in) d														! intent(in) means input
Cf2py intent(out) ret													! intent(out) means only this will be returned
Cf2py depend(d) ret
	  
      ret = 0
	  
	  do i = 1,6														! loop over all and find absolute largest value
		if (abs(d(i)) .gt. ret) then
		ret = abs(d(i))
		endif
	  enddo
	  
      END
	  
	  SUBROUTINE MAX9(ret, d)
	  implicit none
C
C     Maximum of 9x9
C
      REAL*8 d(9,9)														! d = data = pixel values
	  REAL*8 ret														! ret = return
Cf2py intent(in) d														! intent(in) means input
Cf2py intent(out) ret													! intent(out) means only this will be returned
Cf2py depend(d) ret
	  	  
	  ret = maxval(d)
	  
      END
	  
	  SUBROUTINE MAX7(ret, d)
	  implicit none
C
C     Maximum of 7x7
C
      REAL*8 d(7,7)														! d = data = pixel values
	  REAL*8 ret														! ret = return
Cf2py intent(in) d														! intent(in) means input
Cf2py intent(out) ret													! intent(out) means only this will be returned
Cf2py depend(d) ret
	  	  
	  ret = maxval(d)
	  
      END
	  	  
C END FILE MBX_FORTRAN_TEST.F90