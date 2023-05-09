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
		X(i) = (B(i) - dot_product(A(i,1:i-1),X(:)))/A(i,i)
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
      
	
	
	
	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!		OPERATION SUR DES MATRICES CREUSES 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





	FUNCTION Cholesky_creux(n,Ne,vois, pos_diag, A_Creux)

	IMPLICIT NONE 

	INTEGER :: n
	INTEGER :: Ne
	INTEGER, DIMENSION(n) :: vois
	INTEGER, DIMENSION(Ne+1) :: pos_diag 
	REAL(rp), DIMENSION(n) :: A_Creux
	
	
	REAL(rp), DIMENSION(n) :: Cholesky_creux
	REAL(rp), DIMENSION(n) :: L
	REAL(rp) :: temp
	INTEGER :: i, j, k, m, p, compteur 
	
	Cholesky_creux = 0._rp
	L = 0._rp
	DO i = 1, Ne 
		DO j = pos_diag(i) + 1, pos_diag(i+1) - 1
			IF (vois(j) < i) THEN  
				L(pos_diag(i)) = A_creux(pos_diag(i)) - L(j)**2
			END IF 
		END DO 
		
		DO j = pos_diag(i)+1,pos_diag(i+1)-1
			IF(vois(j) > i) THEN 
				DO k = pos_diag(vois(j))+1, pos_diag(vois(j)+1) -1
					IF (vois(k) == i) THEN  	
						DO m = pos_diag(j) + 1,pos_diag(i+1) - 1
							IF(vois(m) < vois(j)-1) THEN 
								p = pos_diag(vois(i))+1
								compteur = 0
								write(*,*) vois(k),vois(p),vois(m),i
								DO WHILE(p /= pos_diag(vois(i)+1) .AND. compteur == 0)  
									IF (vois(p) == vois(m)) THEN 
										L(k) = A_creux(k) - L(p)*L(m)/L(pos_diag(i))
										compteur = 1
									END IF 
									p = p+1
								END DO 
							END IF 
						END DO 
						L(j) = L(k)
					END IF   
				END DO 
			END IF 	 
		END DO 
		
	!	DO j = i+1,n
	!		L(j,i) = (A(j,i) - dot_product(L(i,1:j-1),L(j,1:j-1)))/L(i,i)
	!		L(i,j) = L(j,i) !-- Création de la parite inférieure
	!	END DO
	END DO 
	
	Cholesky_Creux = L	
	RETURN 
END FUNCTION Cholesky_Creux
END MODULE mod_algebre 
