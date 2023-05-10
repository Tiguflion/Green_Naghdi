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





		!-- Fonction qui résouds le problème matriciel AU = L , à l'aide de la matrice A_creux, en utilisant la méthode de Gauss-Seidel : 
! ! @ Variable d'entrées 
 
!				N_coeff : Taille du vecteur des poids A_creux
!				Ns : Nombre de points dans le maillage 
!				A_Creux : Vecteur des poids
!				A_cellule : Table des coefficients non nul de la matrice des poids 
!				Coeff_Diag : Tableau des coefficients diagonaux (= Pos_Creux)
!				Sec_mem : Second membre L
!	@ Variable de sortie
!				U : Solution du problème matriciel AU = L


      SUBROUTINE gauss_seidel_creux (N_Coeff,Ns, A_Creux,Vois,Coeff_Diag,U,B)
      INTEGER, INTENT(IN)  ::  Ns, N_Coeff
      INTEGER :: i
      INTEGER :: j
      INTEGER :: k 
      REAL(rp), DIMENSION(N_Coeff),  INTENT(IN)     ::  A_Creux
      INTEGER, DIMENSION(N_Coeff),  INTENT(IN)     ::  Vois
      INTEGER, DIMENSION(Ns+1), INTENT(IN) :: Coeff_Diag
      REAL(rp), DIMENSION(Ns),    INTENT(IN)::  B
      REAL(rp), DIMENSION(Ns),    INTENT(OUT)  ::  U
		REAL(rp) :: ERR
		REAL(rp) :: S
		REAL(rp) :: R
		k = 0
		ERR = 1._rp
		!U(:) = 0._rp
		do while (ERR > 10E-10 .and. k < 20000)
			k = k+1
			ERR = 0._rp
			do i = 1,Ns
				R = B(i)/A_creux(Coeff_diag(i))
				do j = 1, (Coeff_diag(i+1) - Coeff_diag(i))
					R = R - A_Creux(j-1+Coeff_diag(i))*U(Vois(j-1 + Coeff_diag(i)))/A_creux(coeff_diag(i))
				end do 			
			
				! Code original 
				!S = 0
				
				!do j = 1, (Coeff_diag(i+1) - Coeff_diag(i))
				!	S = S+A_Creux(j-1+Coeff_diag(i))*U(Vois(j-1 + Coeff_diag(i)))
				!end do 
				
				!R = (B(i) - S)/A_Creux(Coeff_diag(i))
				ERR = ERR + R*R
				U(i) = U(i) + R
				ERR = sqrt(ERR)
				
			end do
		end do

!
      RETURN
      END SUBROUTINE gauss_seidel_creux
!




END MODULE mod_algebre 
