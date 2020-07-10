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
	  
	  SUBROUTINE DENSE_DIF(dif, x, stp, comp, d)
	  implicit none
C
C     Dense difference calculation
C	  
	  INTEGER s, s_roi
	  REAL*8 x(5), stp
	  REAL*8 d(9,9)
	  REAL*8 dif(81, 5)
	  REAL*8 h(5), comp(5)
	  REAL*8 x1(5), x2(5), dx
	  REAL*8 f2(81), f1(81), df(81)
	  INTEGER i
	  
Cf2py intent(in) x, stp, d, comp
Cf2py intent(out) dif
Cf2py depend(x) dif	 

	  s = 5
	  s_roi = 9
	  
	  h = ABS(MAX(comp, x)*stp)
	  
	  !open(1, file = 'output2.dat', status='new')
	  !write(1,'(A)') 'node1 node2 node3 node4 node5'
	  
	  do i =1, s
         x1 = x
		 x1(i) = x1(i) - h(i)
		 x2 = x
		 x2(i) = x2(i) + h(i)
		 dx = 2*h(i)
		 call GAUSSIAN(f1, x1(1), x1(2), x1(3), x1(4), x1(5), s_roi, d)
		 call GAUSSIAN(f2, x2(1), x2(2), x2(3), x2(4), x2(5), s_roi, d)
		 df = f2 - f1
		 dif(:,i) = df/dx
		 !write(1,*) df
	  enddo
	  
	  END
	  	  
C END FILE GAUSSIAN_FULL.F90