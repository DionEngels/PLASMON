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
	  
	  SUBROUTINE FFT9(x_re, x_im, y_re, y_im, d)
	  implicit none
C
C     FFT for 9x9 ROI
C
      REAL*8 d(9,9)											! g = gaussian result
	  REAL*8 x_re, y_re, y_im, x_im							! d = data = pixel values
	  REAL*8 fit_o(9), fit_cos(9), fit_sin(9)
	  REAL*8 loop_x_im(9), loop_x_re(9), loop_y_im(9), loop_y_re(9)
	  real*8 pi
	  INTEGER s, i, j
Cf2py intent(in) d											! intent(in) means input
Cf2py intent(out) x_re, y_re, y_im, x_im					! intent(out) means only this will be returned
Cf2py depend(d) x_re, y_re, y_im, x_im
	  	  
	  s = 9
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