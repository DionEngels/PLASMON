C FILE: MBX_FORTRAN.F90
	  SUBROUTINE GAUSSIAN(g, h, c_x, c_y, w_x, w_y, s, d)
	  implicit none
C
C     Make a gaussian without background
C
      REAL*8 h, c_x, c_y, w_x, w_y										! height, center (x and y), and widths
	  REAL*8 expo														! temporary save location of exponent
	  INTEGER s, i, j													! s = roi_size
      REAL*8 g(s*s)														! g = gaussian result
	  REAL*8 d(9,9)														! d = data = pixel values
Cf2py intent(in) h, c_x, c_y, w_x, w_y, s, d							! intent(in) means input
Cf2py intent(out) g														! intent(out) means only this will be returned
Cf2py depend(c_x) g
	  
      do i = 0, s - 1
         do j = 0, s - 1
			expo = h*exp(-(((c_x-i)/w_x)**2+((c_y-j)/w_y)**2)/2)		! calculate exponential
            g(i*s+j+1) = expo-d(i+1,j+1)								! subtract pixel value and store
         enddo
      enddo
      END
	  
	  SUBROUTINE GS_BG(g, h, c_x, c_y, w_x, w_y, b, s, d)
	  implicit none
C
C     Make a gaussian with background
C
      REAL*8 h, c_x, c_y, w_x, w_y, b									! height, center (x and y), widths, and background
	  REAL*8 expo														! temporary save location of exponent
	  INTEGER s, i, j													! s = roi_size
      REAL*8 g(s*s)														! g = gaussian result
	  REAL*8 d(9,9)														! d = data = pixel values
Cf2py intent(in) h, c_x, c_y, w_x, w_y, s, d, b							! intent(in) means input
Cf2py intent(out) g														! intent(out) means only this will be returned
Cf2py depend(c_x) g
      do i = 0, s - 1
         do j = 0, s - 1
			expo = h*exp(-(((c_x-i)/w_x)**2+((c_y-j)/w_y)**2)/2)		! calculate exponential
            g(i*s+j+1) = expo-d(i+1,j+1)+b								! subtract pixel value and add background and store
         enddo
      enddo
      END
	  
	  SUBROUTINE DENSE_DIF(dif, x, stp, comp, s, s2, d)
	  implicit none
C
C     Dense difference calculation
C	  
	  INTEGER s, s2 													! s = num of params (5 in this case), s2 = ROI size
	  REAL*8 x(5), stp													! x = params, stp = EPS^(1/3)
	  REAL*8 d(9,9)														! d = data = pixel values
	  REAL*8 dif(s2*s2, 5)												! dif = dense difference result
	  REAL*8 h(5), comp(5)												! comp = compare (all ones), h = h values (calculated first)
	  REAL*8 x1(5), x2(5), dx											! other params and difference
	  REAL*8 f2(s2*s2), f1(s2*s2), df(s2*s2)							! other function values and difference
	  INTEGER i
	  
Cf2py intent(in) x, stp, d, comp, s, s2
Cf2py intent(out) dif
Cf2py depend(x) dif	 
	  
	  h = ABS(MAX(comp, x)*stp)
	  
	  do i =1, s														! based on scipy.optimize.least_squares (Python)
         x1 = x
		 x1(i) = x1(i) - h(i)
		 x2 = x
		 x2(i) = x2(i) + h(i)
		 dx = 2*h(i)
		 call GAUSSIAN(f1, x1(1), x1(2), x1(3), x1(4), x1(5), s2, d)
		 call GAUSSIAN(f2, x2(1), x2(2), x2(3), x2(4), x2(5), s2, d)
		 df = f2 - f1
		 dif(:,i) = df/dx
	  enddo
	  
	  END
	  
	  SUBROUTINE DENSE_DIF_BG(dif, x, stp, comp, s, s2, d)
	  implicit none
C
C     Dense difference calculation
C	  
	  INTEGER s, s2 													! s = num of params (6 in this case), s2 = ROI size
	  REAL*8 x(6), stp													! x = params, stp = EPS^(1/3)
	  REAL*8 d(9,9)														! d = data = pixel values
	  REAL*8 dif(s2*s2, 6)												! dif = dense difference result
	  REAL*8 h(6), comp(6)												! comp = compare (all ones), h = h values (calculated first)
	  REAL*8 x1(6), x2(6), dx											! other params and difference
	  REAL*8 f2(s2*s2), f1(s2*s2), df(s2*s2)							! other function values and difference
	  INTEGER i
	  
Cf2py intent(in) x, stp, d, comp
Cf2py intent(out) dif
Cf2py depend(x) dif	 
	  
	  h = ABS(MAX(comp, x)*stp)
	  
	  do i =1, s														! based on scipy.optimize.least_squares (Python)
         x1 = x
		 x1(i) = x1(i) - h(i)
		 x2 = x
		 x2(i) = x2(i) + h(i)
		 dx = 2*h(i)
		 call GS_BG(f1, x1(1), x1(2), x1(3), x1(4), x1(5), x1(6), s2, d)
		 call GS_BG(f2, x2(1), x2(2), x2(3), x2(4), x2(5), x2(6), s2, d)
		 df = f2 - f1
		 dif(:,i) = df/dx
	  enddo
	  
	  END
	  
	  SUBROUTINE GAUSSIAN7(g, h, c_x, c_y, w_x, w_y, s, d)
	  implicit none
C
C     Make a gaussian without background, ROI7
C
      REAL*8 h, c_x, c_y, w_x, w_y										! height, center (x and y), and widths
	  REAL*8 expo														! temporary save location of exponent
	  INTEGER s, i, j													! s = roi_size
      REAL*8 g(s*s)														! g = gaussian result
	  REAL*8 d(7,7)														! d = data = pixel values
Cf2py intent(in) h, c_x, c_y, w_x, w_y, s, d							! intent(in) means input
Cf2py intent(out) g														! intent(out) means only this will be returned
Cf2py depend(c_x) g
	  
      do i = 0, s - 1
         do j = 0, s - 1
			expo = h*exp(-(((c_x-i)/w_x)**2+((c_y-j)/w_y)**2)/2)		! calculate exponential
            g(i*s+j+1) = expo-d(i+1,j+1)								! subtract pixel value and store
         enddo
      enddo
      END
	  
	  SUBROUTINE GS_BG7(g, h, c_x, c_y, w_x, w_y, b, s, d)
	  implicit none
C
C     Make a gaussian with background
C
      REAL*8 h, c_x, c_y, w_x, w_y, b									! height, center (x and y), widths, and background
	  REAL*8 expo														! temporary save location of exponent
	  INTEGER s, i, j													! s = roi_size
      REAL*8 g(s*s)														! g = gaussian result
	  REAL*8 d(7,7)														! d = data = pixel values
Cf2py intent(in) h, c_x, c_y, w_x, w_y, s, d, b							! intent(in) means input
Cf2py intent(out) g														! intent(out) means only this will be returned
Cf2py depend(c_x) g
      do i = 0, s - 1
         do j = 0, s - 1
			expo = h*exp(-(((c_x-i)/w_x)**2+((c_y-j)/w_y)**2)/2)		! calculate exponential
            g(i*s+j+1) = expo-d(i+1,j+1)+b								! subtract pixel value and add background and store
         enddo
      enddo
      END
	  
	  SUBROUTINE DENSE_DIF7(dif, x, stp, comp, s, s2, d)
	  implicit none
C
C     Dense difference calculation
C	  
	  INTEGER s, s2 													! s = num of params (5 in this case), s2 = ROI size
	  REAL*8 x(5), stp													! x = params, stp = EPS^(1/3)
	  REAL*8 d(7,7)														! d = data = pixel values
	  REAL*8 dif(s2*s2, 5)												! dif = dense difference result
	  REAL*8 h(5), comp(5)												! comp = compare (all ones), h = h values (calculated first)
	  REAL*8 x1(5), x2(5), dx											! other params and difference
	  REAL*8 f2(s2*s2), f1(s2*s2), df(s2*s2)							! other function values and difference
	  INTEGER i
	  
Cf2py intent(in) x, stp, d, comp, s, s2
Cf2py intent(out) dif
Cf2py depend(x) dif	 
	  
	  h = ABS(MAX(comp, x)*stp)
	  
	  do i =1, s														! based on scipy.optimize.least_squares (Python)
         x1 = x
		 x1(i) = x1(i) - h(i)
		 x2 = x
		 x2(i) = x2(i) + h(i)
		 dx = 2*h(i)
		 call GAUSSIAN7(f1, x1(1), x1(2), x1(3), x1(4), x1(5), s2, d)
		 call GAUSSIAN7(f2, x2(1), x2(2), x2(3), x2(4), x2(5), s2, d)
		 df = f2 - f1
		 dif(:,i) = df/dx
	  enddo
	  
	  END
	  
	  SUBROUTINE DENSE_DIF_BG7(dif, x, stp, comp, s, s2, d)
	  implicit none
C
C     Dense difference calculation
C	  
	  INTEGER s, s2 													! s = num of params (6 in this case), s2 = ROI size
	  REAL*8 x(6), stp													! x = params, stp = EPS^(1/3)
	  REAL*8 d(7,7)														! d = data = pixel values
	  REAL*8 dif(s2*s2, 6)												! dif = dense difference result
	  REAL*8 h(6), comp(6)												! comp = compare (all ones), h = h values (calculated first)
	  REAL*8 x1(6), x2(6), dx											! other params and difference
	  REAL*8 f2(s2*s2), f1(s2*s2), df(s2*s2)							! other function values and difference
	  INTEGER i
	  
Cf2py intent(in) x, stp, d, comp
Cf2py intent(out) dif
Cf2py depend(x) dif	 
	  
	  h = ABS(MAX(comp, x)*stp)
	  
	  do i =1, s														! based on scipy.optimize.least_squares (Python)
         x1 = x
		 x1(i) = x1(i) - h(i)
		 x2 = x
		 x2(i) = x2(i) + h(i)
		 dx = 2*h(i)
		 call GS_BG7(f1, x1(1), x1(2), x1(3), x1(4), x1(5), x1(6), s2, d)
		 call GS_BG7(f2, x2(1), x2(2), x2(3), x2(4), x2(5), x2(6), s2, d)
		 df = f2 - f1
		 dif(:,i) = df/dx
	  enddo
	  
	  END
	  	  
C END FILE MBX_FORTRAN.F90