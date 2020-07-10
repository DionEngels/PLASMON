C FILE: GAUSSIAN_FULL.F90
	  SUBROUTINE GAUSSIAN(g, h, c_x, c_y, w_x, w_y, s, d)
	  implicit none
C
C     Make a gaussian without background
C
      REAL*8 h, c_x, c_y, w_x, w_y
	  REAL*8 expo
	  INTEGER s, i, j
      REAL*8 g(s*s), d(9,9)
Cf2py intent(in) h, c_x, c_y, w_x, w_y, s, d
Cf2py intent(out) g
Cf2py depend(c_x) g
      do i = 0, s - 1
         do j = 0, s - 1
			expo = h*exp(-(((c_x-i)/w_x)**2+((c_y-j)/w_y)**2)/2)
            g(i*s+j+1) = expo-d(i+1,j+1)
         enddo
      enddo
      END
	  
	  SUBROUTINE GAUSSIAN_BACK(g, h, c_x, c_y, w_x, w_y, b, s, d)
	  implicit none
C
C     Make a gaussian with background
C
      REAL*8 h, c_x, c_y, w_x, w_y, b
	  REAL*8 expo
	  INTEGER s, i, j
      REAL*8 g(s*s), d(9,9)
Cf2py intent(in) h, c_x, c_y, w_x, w_y, s, d, b
Cf2py intent(out) g
Cf2py depend(c_x) g
      do i = 0, s - 1
         do j = 0, s - 1
			expo = h*exp(-(((c_x-i)/w_x)**2+((c_y-j)/w_y)**2)/2)
            g(i*s+j+1) = expo-d(i+1,j+1)+b
         enddo
      enddo
      END
	  
C END FILE GAUSSIAN_FULL.F90