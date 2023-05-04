MODULE mod_correction

USE numerics 
USE mod_algebre
USE mod_cond
IMPLICIT NONE 

CONTAINS 

SUBROUTINE Correction(W,Ne)

	IMPLICIT NONE 
	
	INTEGER :: Ne ! Nombre d'éléments à corriger 
	REAL(rp), DIMENSION(:,:) :: W
	REAL(rp), DIMENSION(:,:),  ALLOCATABLE :: sol
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: A
	REAL(rp),DIMENSION(:),ALLOCATABLE :: U ! Vitesse horizontale 
	REAL(rp),DIMENSION(:),ALLOCATABLE :: omega ! Vitesse verticale
	REAL(rp),DIMENSION(:),ALLOCATABLE :: sigma ! Ecart-type 
	REAL(rp), DIMENSION(:), ALLOCATABLE :: B ! Second membre
	REAL(rp), DIMENSION(:), ALLOCATABLE :: A_diag !- Coefficients diagonale de la matrice A
	REAL(rp), DIMENSION(:), ALLOCATABLE :: um2 
	REAL(rp), DIMENSION(:), ALLOCATABLE :: up2
	
	REAL(rp), DIMENSION(Ne) :: vec_debug

	INTEGER i,j

	CALL Init_Vec_Cor(Ne,sol,A,B,U,A_diag,um2,up2,omega,sigma)


	A = Diag(A_diag,Ne,0) + Diag(up2,Ne,-2) + Diag(up2,Ne, 2)
	

	CALL Cond_Bord_coeff(Ne,A,B) !--- Modification des coefficients liés aux conditions aux bords

	
	!DO i = 1,Nx
	!	IF(B(i) /= 0._rp) THEN 
	!		WRITE(*,*) i, B(i)
	!	END IF 
	!END DO 
	CALL res_Chol(Ne, A, U, B)	
	
	

	!CALL gauss_seidel (Nx, A, U, B) !Faire un Choleksy
	!WRITE(*,*) 'B :'
	
	!CALL u_theo(X,date,U)	
	!CALL w_theo(X,date,omega)
	!CALL sig_theo(X,date,sigma)
	
	!vec_debug = matmul(A,U)
	!DO i = 1,Nx
	!	IF(abs(B(i)-vec_debug(i)) > 1E-5) then 
	!		WRITE(*,*) i, B(i) - vec_debug(i), B(i), vec_debug(i)
	!	END IF 
	!END DO 
	!READ(*,*)
	
	!WRITE(*,*) norm2(matmul(A,U)-B)
	!READ(*,*)
	CALL omega_sigma(Ne,U,omega,sigma)

	CALL bord_omega_sigma(Ne,U,omega,sigma)
	
	W(:,:) = assemble_correction(Ne,U,omega,sigma)
	
	!WRITE(*,*) sum(W(1,:))

END SUBROUTINE Correction 



SUBROUTINE Init_Vec_Cor(Ne,sol,A,B,U,A_diag,um2,up2,omega,sigma)

IMPLICIT NONE 
	INTEGER,			       INTENT(IN)    :: Ne
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: sol
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: A
	REAL(rp), DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: U
	REAL(rp), DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: B
	REAL(rp), DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: A_diag
	REAL(rp), DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: um2
	REAL(rp), DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: up2	
	REAL(rp), DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: omega
	REAL(rp), DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: sigma
	
	
	
	INTEGER :: i 
	
	REAL :: cl, cr
	
! Faire une condition avec une variable choix qui vérifie si les variables sont allouées

	ALLOCATE(sol(4,Ne))
	ALLOCATE(A(Ne,Ne))
	ALLOCATE(B(Ne),A_diag(Ne),U(Ne))
	ALLOCATE(um2(Ne),up2(ne))
	ALLOCATE(omega(Ne),sigma(Ne))

	sol(:,:) = 0._rp
	A(:,:) = 0._rp
	B(:) = 0._rp
	A_diag(:) = 0._rp
	um2(:) = 0._rp
	up2(:) = 0._rp
	omega(:) = 0._rp
	sigma(:) = 0._rp
	U(:) = 0._rp
	
	DO i = 2,Ne-1
		A_diag(i) = W(1,i) + 1._rp/(12._rp*dx**2)*(W(1,i+1)**3 + W(1,i-1)**3)
		up2(i) =-1._rp/(12._rp*dx**2)*(W(1,i+1)**3)
		cl = W(1,i-1)*(sqrt(3._rp)*W(3,i-1) + W(4,i-1))
		cr = W(1,i+1)*(sqrt(3._rp)*W(3,i+1) + W(4,i+1))
		B(i) = W(2,i) + (cr-cl)/(4._rp*sqrt(3._rp)*dx)
	END DO 
	
END SUBROUTINE Init_Vec_Cor



SUBROUTINE Cond_bord_Coeff(Ne,A,B) 

	IMPLICIT NONE 
	
	INTEGER,		  INTENT(IN) 	:: Ne
	REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: A
	REAL(rp), DIMENSION(:),   INTENT(INOUT) :: B
	
	
	REAL(rp) :: c1,c2,c3,cN2,cN1,cN
	INTEGER :: i,j
	
	!=========================================================================
	!
	!			CONDITIONS AU BORD GAUCHE 
	!
	!=========================================================================
	
	A(1,1) = W(1,1) + (W(1,2)**3 + W(1,1)**3)/(12._rp*dx**2)
	A(1,2) = 	       (W(1,1)**3)       /(12._rp*dx**2) !- retirer le "-"
	A(1,3) = 	      - (W(1,2)**3)	 /(12._rp*dx**2)
	A(2,2) = W(1,2) + (W(1,3)**3 + W(1,1)**3)/(12._rp*dx**2)
	A(2,4) =              -(W(1,3)**3)	 /(12._rp*dx**2)
	
	A(2,1) = A(1,2)
	A(3,1) = A(1,3)
	A(4,2) = A(2,4)
	

	c1 = W(1,1)*(sqrt(3._rp)*W(3,1) + W(4,1))
	c2 = W(1,2)*(sqrt(3._rp)*W(3,2) + W(4,2))
	c3 = W(1,3)*(sqrt(3._rp)*W(3,3) + W(4,3))
	
	B(1) = W(2,1) + (c2 - c1)/(4._rp*sqrt(3._rp)*dx)
	B(2) = W(2,2) + (c3 - c1)/(4._rp*sqrt(3._rp)*dx)

	!=========================================================================
	!
	!			CONDITIONS AU BORD GAUCHE 
	!
	!=========================================================================
	
	A(Ne,Ne)    = W(1,Ne)    + (W(1,Ne)**3 + W(1,Ne-1)**3)  /(12._rp*dx**2)
	A(Ne,Ne-1)  =		    W(1,Ne)**3	 	        /(12._rp*dx**2) !-- Retirer le "-" à cause des conditions aux bords 
	A(Ne,Ne-2)  =		 -W(1,Ne-1)**3 			/(12._rp*dx**2)
	A(Ne-1,Ne-1) = W(1,Ne-1) + (W(1,Ne)**3 + W(1,Ne-2)**3)  /(12._rp*dx**2)
	A(Ne-1,Ne-3) = 		-(W(1,Ne - 2)**3)	        /(12._rp*dx**2)	
	
	
	A(Ne-1,Ne) = A(Ne,Ne-1) 
	A(Ne-2,Ne) =  A(Ne,Ne-2) 
	A(Ne-3,Ne-1) = A(Ne-1,Ne-3)	


	cN =  W(1,Ne)  *(sqrt(3._rp)*W(3,Ne  ) + W(4,Ne  ))
	cN1 = W(1,Ne-1)*(sqrt(3._rp)*W(3,Ne-1) + W(4,Ne-1))
	cN2 = W(1,Ne-2)*(sqrt(3._rp)*W(3,Ne-2) + W(4,Ne-2))


	
	B(Ne)   = W(2,Ne)     + (cN - cN1)/(4._rp*sqrt(3._rp)*dx)
	B(Ne-1) = W(2,Ne - 1) + (cN - cN2)/(4._rp*sqrt(3._rp)*dx)

END SUBROUTINE Cond_bord_Coeff




SUBROUTINE omega_sigma(Ne,U,omega,sigma)

IMPLICIT NONE 

INTEGER,		INTENT(IN) :: Ne
REAL(rp), DIMENSION(:), INTENT(IN) :: U
REAL(rp), DIMENSION(:), INTENT(INOUT) :: omega
REAL(rp), DIMENSION(:), INTENT(INOUT) :: sigma

INTEGER :: i

DO i = 2,Ne-1
	omega(i) = -W(1,i)*(U(i+1) - U(i-1))/(4._rp*dx) !-- Gradiant de bathymétrie ici, négligé car on est sur fond plat
	sigma(i) = -W(1,i)*(U(i+1) - U(i-1))/(4._rp*dx*sqrt(3._rp))
END DO 

END SUBROUTINE omega_sigma




SUBROUTINE bord_omega_sigma(Ne,U,omega,sigma)

IMPLICIT NONE 

INTEGER,		INTENT(IN) :: Ne
REAL(rp), DIMENSION(:), INTENT(IN) :: U
REAL(rp), DIMENSION(:), INTENT(INOUT) :: omega
REAL(rp), DIMENSION(:), INTENT(INOUT) :: sigma

INTEGER :: i

omega(1) = -W(1,1)*(U(2) + U(1))/(4._rp*dx)
sigma(1) = -W(1,1)*(U(2) + U(1))/(4._rp*dx*sqrt(3._rp))

omega(Ne) = W(1,Ne)*(U(Ne) + U(Ne-1))/(4._rp*dx)
sigma(Ne) = W(1,Ne)*(U(Ne) + U(Ne-1))/(4._rp*dx*sqrt(3._rp))

END SUBROUTINE bord_omega_sigma





FUNCTION Diag(vec,Ne,k) 

IMPLICIT NONE 
	INTEGER :: Ne
	REAL(rp), DIMENSION(:) :: vec
	INTEGER :: k
	REAL(rp), DIMENSION(Ne,Ne) :: Diag
	INTEGER :: i
	
	Diag(:,:) = 0._rp
	IF(k <= 0) THEN
		DO i = 1, Ne+k
			Diag(i,i-k) = vec(i)
		END DO 
	ELSE IF(k > 0) THEN
		DO i = 1, Ne - k
			Diag(i+k,i) = vec(i)
		END DO 
	END IF 

END FUNCTION Diag
	



	
FUNCTION assemble_correction(Ne,U,omega,sigma)

IMPLICIT NONE 

	INTEGER :: Ne
	REAL(rp),DIMENSION(Ne) :: U
	REAL(rp),DIMENSION(Ne) :: omega
	REAL(rp),DIMENSION(Ne) :: sigma
	
	REAL(rp), DIMENSION(4,Ne) :: assemble_correction
	
	REAL(rp), DIMENSION(Ne) :: h
	h(:) = 0._rp
	h(:) = W(1,:)
	
	assemble_correction(1,:) = h(:)
	assemble_correction(2,:) = h(:)*U(:)
	assemble_correction(3,:) = h(:)*omega(:)
	assemble_correction(4,:) = h(:)*sigma(:)
	
END FUNCTION assemble_correction
	




	
		
	
	
END MODULE mod_correction
