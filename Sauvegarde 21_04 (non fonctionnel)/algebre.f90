MODULE mod_algebre

USE numerics

IMPLICIT NONE 

CONTAINS 

   SUBROUTINE gauss_seidel (n, A, X, B)
      INTEGER, INTENT(IN)  ::  n
      INTEGER :: i
      INTEGER :: j
      INTEGER :: k = 0
      REAL(rp), DIMENSION(n,n),  INTENT(IN)     ::  A
      REAL(rp), DIMENSION(n),    INTENT(INOUT)  ::  X
      REAL(rp), DIMENSION(n),    INTENT(IN)     ::  B
		REAL :: ERR = 1
		REAL :: S
		REAL :: R
		ERR = 1
		k = 0

		do while (ERR > 10E-8 .and. k < 10000)
		k = k+1
			ERR = 0._rp
			do i = 1,n
				S = dot_product(A(i,:),X(:))
				R = (B(i) - S)/A(i,i)
				ERR = ERR + R**2
				X(i) = X(i) + R
			end do
			
		end do

!
      RETURN
      END SUBROUTINE gauss_seidel


END MODULE mod_algebre 
