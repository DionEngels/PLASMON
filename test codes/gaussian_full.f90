C FILE: GAUSSIAN_FULL.F90
      SUBROUTINE GAUSSIAN(g, h, c_x, c_y, w_x, w_y, s)
	  implicit none
C
C     Make a gaussian
C
      REAL*8 h, c_x, c_y, w_x, w_y
	  INTEGER s, i, j
      REAL*8 g(s*s)
Cf2py intent(in) h, c_x, c_y, w_x, w_y, s
Cf2py intent(out) g
Cf2py depend(c_x) g
      do i = 0, s - 1
         do j = 0, s - 1
            g(i*s+j+1) = h*exp(-(((c_x-i)/w_x)**2+((c_y-j)/w_y)**2)/2)
         enddo
      enddo
      END
	  
	  SUBROUTINE GAUSSIAN_BACKGROUND(g, h, c_x, c_y, w_x, w_y, b, s)
	  implicit none
C
C     Make a gaussian with background
C
      REAL*8 h, c_x, c_y, w_x, w_y, b
	  INTEGER s, i, j
      REAL*8 g(s*s)
Cf2py intent(in) h, c_x, c_y, w_x, w_y, s
Cf2py intent(out) g
Cf2py depend(c_x) g
      do i = 0, s - 1
         do j = 0, s - 1
            g(i*s+j+1) = h*exp(-(((c_x-i)/w_x)**2+((c_y-j)/w_y)**2)/2)+b
         enddo
      enddo
      END
	  
	  SUBROUTINE GAUSSIAN_VEC(g, h, c_x, c_y, w_x, w_y, s_int, s)
C
C     Make a gaussian without background
C
      REAL*8 h, c_x, c_y, w_x, w_y
	  integer s_int
	  REAL*8 s(2, s_int*s_int)
Cf2py intent(in) h, c_x, c_y, w_x, w_y, s
	  REAL*8 g(s_int*s_int)
Cf2py intent(out) g
Cf2py depend(s) g
	  
	  g = h*exp(-(((c_x-s(1,:))/w_x)**2+((c_y-s(2,:))/w_y)**2)/2)
  
      END
	  
C END FILE GAUSSIAN_FULL.F90