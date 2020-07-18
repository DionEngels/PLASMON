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
	  
	  SUBROUTINE MIN9(ret, d)
	  implicit none
C
C     Minimum of 9x9
C
      REAL*8 d(9,9)														! d = data = pixel values
	  REAL*8 ret														! ret = return
Cf2py intent(in) d														! intent(in) means input
Cf2py intent(out) ret													! intent(out) means only this will be returned
Cf2py depend(d) ret
	  	  
	  ret = minval(d)
	  
      END
	  
	  SUBROUTINE MIN7(ret, d)
	  implicit none
C
C     Minimum of 7x7
C
      REAL*8 d(7,7)														! d = data = pixel values
	  REAL*8 ret														! ret = return
Cf2py intent(in) d														! intent(in) means input
Cf2py intent(out) ret													! intent(out) means only this will be returned
Cf2py depend(d) ret
	  	  
	  ret = minval(d)
	  
      END
	  
	  SUBROUTINE FFT9(x_re, x_im, y_re, y_im, d)
	  implicit none
C
C     FFT for 9x9 ROI
C
      REAL*8 d(9,9)														! d = data = pixel values
	  REAL*8 x_re, y_re, y_im, x_im										! real and imaginary parts
	  REAL*8 fit_o(9), fit_cos(9), fit_sin(9)							! omega, sin and cosine
	  REAL*8 loop_x_im(9), loop_x_re(9), loop_y_im(9), loop_y_re(9)		! loop variables
	  real*8 pi
	  INTEGER s, i, j
Cf2py intent(in) d														! intent(in) means input
Cf2py intent(out) x_re, y_re, y_im, x_im								! intent(out) means only this will be returned
Cf2py depend(d) x_re, y_re, y_im, x_im
	  	  
	  s = 9
	  pi = 2.d0 * asin(1.d0)											! define pi
	  
	  do i = 1, s
	  fit_o(i) = i*2*pi/s
	  fit_cos(i) = cos(fit_o(i))
	  fit_sin(i) = sin(fit_o(i))
	  loop_x_re(i) = 0													! predefine
	  loop_x_im(i) = 0
	  loop_y_re(i) = 0
	  loop_y_im(i) = 0
	  enddo
	  
	  do i=1, s
	  do j=1, s
	  loop_x_re(j) = loop_x_re(j) + fit_cos(j)*d(i,j)	
	  loop_x_im(j) = loop_x_im(j) - fit_sin(j)*d(i,j)
	  loop_y_re(i) = loop_y_re(i) + fit_cos(i)*d(i,j)
	  loop_y_im(i) = loop_y_im(i) - fit_sin(i)*d(i,j)
	  enddo
	  enddo
	  
	  x_re = sum(loop_x_re)												! sum and return
	  x_im = sum(loop_x_im)
	  y_re = sum(loop_y_re)
	  y_im = sum(loop_y_im)
	  
      END
	  
	  SUBROUTINE FFT7(x_re, x_im, y_re, y_im, d)
	  implicit none
C
C     FFT for 7x7 ROI
C
      REAL*8 d(7,7)											! g = gaussian result
	  REAL*8 x_re, y_re, y_im, x_im							! d = data = pixel values
	  REAL*8 fit_o(7), fit_cos(7), fit_sin(7)
	  REAL*8 loop_x_im(7), loop_x_re(7), loop_y_im(7), loop_y_re(7)
	  real*8 pi
	  INTEGER s, i, j
Cf2py intent(in) d											! intent(in) means input
Cf2py intent(out) x_re, y_re, y_im, x_im					! intent(out) means only this will be returned
Cf2py depend(d) x_re, y_re, y_im, x_im
	  	  
	  s = 7
	  pi = 2.d0 * asin(1.d0)
	  
	  do i = 1, s
	  fit_o(i) = i*2*pi/s
	  fit_cos(i) = cos(fit_o(i))
	  fit_sin(i) = sin(fit_o(i))
	  loop_x_re(i) = 0
	  loop_x_im(i) = 0
	  loop_y_re(i) = 0
	  loop_y_im(i) = 0
	  enddo
	  
	  do i=1, s
	  do j=1, s
	  loop_x_re(j) = loop_x_re(j) + fit_cos(j)*d(i,j)	
	  loop_x_im(j) = loop_x_im(j) - fit_sin(j)*d(i,j)
	  loop_y_re(i) = loop_y_re(i) + fit_cos(i)*d(i,j)
	  loop_y_im(i) = loop_y_im(i) - fit_sin(i)*d(i,j)
	  enddo
	  enddo
	  
	  x_re = sum(loop_x_re)
	  x_im = sum(loop_x_im)
	  y_re = sum(loop_y_re)
	  y_im = sum(loop_y_im)
	  
      END
	  	  
C END FILE MBX_FORTRAN_TEST.F90