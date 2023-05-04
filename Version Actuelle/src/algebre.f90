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
				ERR = sqrt(ERR)
			end do
			
		end do

!
      RETURN
      END SUBROUTINE gauss_seidel


! Fonction de la décomposition de Cholesky d'une matrice A
! ! @ Variable d'entrées 
! 		n  : Rang de la matrice 
!		A : La matrice à décomposer
! @ Sortie :
!		Cholesky : La solution du problème (Contenant les deux matrices triangulaires)

	FUNCTION Cholesky(n,A)

	IMPLICIT NONE 

	INTEGER :: n
	REAL(rp), DIMENSION(n,n) :: A
	
	REAL(rp), DIMENSION(n,n) :: Cholesky
	REAL(rp), DIMENSION(n,n) :: L
	REAL(rp) :: temp
	INTEGER :: i, j, k
	Cholesky = 0._rp
	L = 0._rp
	DO i = 1, n 
		L(i,i) = sqrt(A(i,i) - dot_product(L(i,1:i-1),L(i,1:i-1)))
		DO j = i+1,n
			L(j,i) = (A(j,i) - dot_product(L(i,1:j-1),L(j,1:j-1)))/L(i,i)
			L(i,j) = L(j,i) !-- Création de la parite inférieure
		END DO
	END DO 
	
	Cholesky = L	
	RETURN 
END FUNCTION Cholesky	
	
	
		
! Fonction qui résoud le système triangulaire supérieur AX = B
! ! @ Variable d'entrées 
! 				n  : Le rang de la matrice 
!				A : Une matrice triangulaire supérieure
!				B : Le second membre du problème
! @ Sortie :
! 				X : La solution du problème		

	SUBROUTINE sys_trig_sup(n,A,X,B)
	
	IMPLICIT NONE 

	INTEGER, INTENT(IN)  ::  n
	INTEGER :: i
	REAL(rp), DIMENSION(n,n),  INTENT(IN)     ::  A
	REAL(rp), DIMENSION(n),  INTENT(INOUT)  :: X
	REAL(rp), DIMENSION(n),    INTENT(IN)     ::  B	
	
	DO i = n,1,-1
		X(i) = (B(i) - dot_product(A(i+1:n,i),X(i+1:n)))/A(i,i)
	END DO 


	END SUBROUTINE sys_trig_sup
	
	
		
! Fonction qui résoud le système triangulaire inférieur AX = B
! ! @ Variable d'entrées 
! 				n  : Le rang de la matrice 
!				A : Une matrice triangulaire inférieure
!				B : Le second membre du problème
! @ Sortie :
! 				X : La solution du problème			
	
	
	SUBROUTINE sys_trig_inf(n,A,X,B)
	
	IMPLICIT NONE 

	INTEGER, INTENT(IN)  ::  n
	INTEGER :: i
	INTEGER :: j
	REAL :: S
	REAL(rp), DIMENSION(n,n),  INTENT(IN)     ::  A
	REAL(rp), DIMENSION(n),    INTENT(INOUT)  ::  X
	REAL(rp), DIMENSION(n),    INTENT(IN)     ::  B	
	
	
	DO i = 1,n
		X(i) = (B(i) - dot_product(A(i,1:j-1),X(:)))/A(i,i)
	END DO 	
	
	END SUBROUTINE sys_trig_inf
	
! Fonction qui résoud le système matriciel AX = B, à l'aide de la méthode d'une décomposition de Cholesky
! ! @ Variable d'entrées 
! 				n  : Le rang de la matrice A
!				A : La matrice du problème
!				B : Le second membre du problème
! @ Sortie :
! 				X : La solution du problème

      SUBROUTINE res_Chol(n, A, X, B)
      
      IMPLICIT NONE 
      
      INTEGER, INTENT(IN)  ::  n
      INTEGER :: i
      INTEGER :: j
      REAL(rp), DIMENSION(n,n),  INTENT(IN)     ::  A
      REAL(rp), DIMENSION(n),    INTENT(INOUT)  ::  X
      REAL(rp), DIMENSION(n),    INTENT(IN)     ::  B
 
	REAL(rp), DIMENSION(n,n) :: LtL 
	REAL(rp), DIMENSION(n) :: Ly
 	Ly = 0._rp
      LtL = Cholesky(n,A)

      CALL sys_trig_inf(n,LtL,Ly,B)
      
      CALL sys_trig_sup(n,LtL,X,Ly)

!
      RETURN
      END SUBROUTINE res_Chol
      
	
	

END MODULE mod_algebre 
