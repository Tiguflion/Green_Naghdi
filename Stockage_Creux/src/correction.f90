MODULE mod_correction

USE numerics 
USE mod_algebre
USE mod_operateur 
USE mod_cond
IMPLICIT NONE 

CONTAINS 

SUBROUTINE Correction(W,ind_min,ind_max,Ne)

	IMPLICIT NONE 
	
	INTEGER :: ind_min, ind_max 
	INTEGER :: Nb_Vois
	INTEGER :: Ne ! Nombre d'éléments à corriger 
	INTEGER , DIMENSION(:)  , ALLOCATABLE :: Pos_diag !-- Position des coefficients diagonaux dans A_Creux
	INTEGER , DIMENSION(:)  , ALLOCATABLE :: Vois
	REAL(rp), DIMENSION(:,:)	      :: W
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: A
	REAL(rp), DIMENSION(:)  , ALLOCATABLE :: A_Creux
	REAL(rp), DIMENSION(:)  , ALLOCATABLE :: U ! Vitesse horizontale 
	REAL(rp), DIMENSION(:)  , ALLOCATABLE :: omega ! Vitesse verticale
	REAL(rp), DIMENSION(:)  , ALLOCATABLE :: sigma ! Ecart-type 
	REAL(rp), DIMENSION(:)  , ALLOCATABLE :: B ! Second membre
	REAL(rp), DIMENSION(:)  , ALLOCATABLE :: A_diag !- Coefficients diagonale de la matrice A
	
	INTEGER, DIMENSION(:)	, ALLOCATABLE :: vois_debug 
	REAL(rp), DIMENSION(:)  , ALLOCATABLE :: vec_debug
	REAL(rp), DIMENSION(Ne,Ne) :: Mat_debug 

	INTEGER i,j,k
	


	CALL Init_Vec_Cor(Ne,B,U,omega,sigma)
	
	CALL init_creux(Ne,Nb_vois,Vois,Pos_diag,A_creux)
	
	CALL rempli_A_creux(Ne,ind_min,ind_max,pos_diag,vois,A_creux)

	CALL rempli_B(Ne,ind_min,ind_max,dx,B)
		
	CALL Cond_Bord_coeff(Ne,B) !--- Modification des coefficients liés aux conditions aux bords (à retravailler)
	
	CALL Cond_bord_Creux(Ne,W(1,ind_min:ind_max),pos_diag,Nb_vois,vois,A_creux)
		
	CALL gauss_seidel_creux(Nb_vois,Ne, A_Creux,Vois,pos_diag,U,B)
	
	CALL omega_sigma(Ne,U,omega,sigma)

	CALL bord_omega_sigma(Ne,U,omega,sigma)
	
	W(:,:) = assemble_correction(Ne,U,omega,sigma)

END SUBROUTINE Correction 



SUBROUTINE Init_Vec_Cor(Ne,B,U,omega,sigma)

IMPLICIT NONE 
	INTEGER,			       INTENT(IN)    :: Ne
	REAL(rp), DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: U
	REAL(rp), DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: B	
	REAL(rp), DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: omega
	REAL(rp), DIMENSION(:),   ALLOCATABLE, INTENT(INOUT) :: sigma	
	
	INTEGER :: i 
	
	
	ALLOCATE(B(Ne),U(Ne))
	ALLOCATE(omega(Ne),sigma(Ne))
	B(:) = 0._rp
	omega(:) = 0._rp
	sigma(:) = 0._rp
	U(:) = 0._rp
	

END SUBROUTINE Init_Vec_Cor



SUBROUTINE init_Creux(Ne,Nb_vois,Vois,Pos_diag,A_creux)

	INTEGER, INTENT(IN) :: Ne
	INTEGER, INTENT(INOUT) :: Nb_vois
	INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Vois 
	INTEGER, DIMENSION(:), ALLOCATABLE , INTENT(INOUT):: Pos_diag 
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: A_creux 
	
	INTEGER :: i
	
	ALLOCATE(Pos_diag(Ne+1))
	
	DO i = 1,Ne+1
		Pos_diag(i) = i
	END DO 
	
	
	Nb_Vois = Ne+1 !-- Prise en compte des éléments diagonaux 
	
	DO i = 1,Ne 
	!--Pour chaques opérateurs, on ne compte pas l'élément diagonal, car il est déjà présent dans Nb_Vois
		CALL Vois_DaDX(i,Ne,Nb_Vois,pos_diag(:)) 
	END DO 
	
	!-- Ajouter une partie qui compte les voisins sur une condition aux bords 
	
	ALLOCATE(A_Creux(Nb_Vois))
	ALLOCATE(Vois(Nb_Vois))
	
	Vois(:) = 0
	A_Creux(:) = 0._rp

	
	DO i = 1,Ne+1
		Vois(pos_diag(i)) = i
	END DO 

END SUBROUTINE Init_Creux



SUBROUTINE rempli_A_creux(Ne,ind_min,ind_max,pos_diag,vois,A_creux)
	INTEGER, INTENT(IN) :: Ne
	INTEGER, INTENT(IN) :: ind_min, ind_max
	INTEGER, DIMENSION(:), INTENT(INOUT) :: pos_diag
	INTEGER, DIMENSION(:), INTENT(INOUT) :: vois
	REAL(rp), DIMENSION(:), INTENT(INOUT) :: A_creux 
	
	INTEGER :: i
	!--- Si on utilise plusieurs opérateurs, ajouter une fonction de tri en dessous de diag_creux
	DO i = ind_min+2,ind_max-2
		CALL DaDX_Creux(i,Ne,-(W(1,ind_min:ind_max)**3/3._rp),pos_diag,vois,A_creux)
	END DO 
	
	CALL Diag_Creux(W(1,ind_min:ind_max),Ne, 0, Pos_diag,Vois, A_creux)
	
	
END SUBROUTINE rempli_A_creux 
	

SUBROUTINE rempli_B(Ne,ind_min,ind_max,dx,B)
	INTEGER, INTENT(IN) :: Ne
	INTEGER, INTENT(IN) :: ind_min
	INTEGER, INTENT(IN) :: ind_max
	REAL(rp), INTENT(IN) :: dx
	REAL(rp), DIMENSION(:), INTENT(INOUT) :: B
	
	
B = W(2,ind_min:ind_max) + Da(W(1,ind_min:ind_max)*(W(4,ind_min:ind_max) + sqrt(3._rp)*W(3,ind_min:ind_max))/(2*sqrt(3._rp)),dx,Ne)


END SUBROUTINE rempli_B



SUBROUTINE Cond_bord_Coeff(Ne,B) 

	IMPLICIT NONE 
	
	INTEGER,		  INTENT(IN) 	:: Ne
	REAL(rp), DIMENSION(:),   INTENT(INOUT) :: B
	
	
	REAL(rp) :: c1,c2,c3,cN2,cN1,cN
	INTEGER :: i,j

	c1 = W(1,1)*(sqrt(3._rp)*W(3,1) + W(4,1))
	c2 = W(1,2)*(sqrt(3._rp)*W(3,2) + W(4,2))
	c3 = W(1,3)*(sqrt(3._rp)*W(3,3) + W(4,3))
	
	B(1) = W(2,1) + (c2 - c1)/(4._rp*sqrt(3._rp)*dx)
	B(2) = W(2,2) + (c3 - c1)/(4._rp*sqrt(3._rp)*dx)

	cN =  W(1,Ne)  *(sqrt(3._rp)*W(3,Ne  ) + W(4,Ne  ))
	cN1 = W(1,Ne-1)*(sqrt(3._rp)*W(3,Ne-1) + W(4,Ne-1))
	cN2 = W(1,Ne-2)*(sqrt(3._rp)*W(3,Ne-2) + W(4,Ne-2))


	
	B(Ne)   = W(2,Ne)     + (cN - cN1)/(4._rp*sqrt(3._rp)*dx)
	B(Ne-1) = W(2,Ne - 1) + (cN - cN2)/(4._rp*sqrt(3._rp)*dx)

END SUBROUTINE Cond_bord_Coeff




SUBROUTINE Cond_bord_Creux(Ne,Vec,pos_diag,Nb_vois,vois,A_creux)
	INTEGER, INTENT(IN) :: Ne
	INTEGER, INTENT(IN) :: Nb_vois
	INTEGER, DIMENSION(:),  INTENT(IN)    :: pos_diag
	REAL(rp), DIMENSION(:), INTENT(IN)    :: Vec
	INTEGER, DIMENSION(:),  INTENT(INOUT) :: vois
	REAL(rp), DIMENSION(:), INTENT(INOUT) :: A_creux
	
	INTEGER :: i,j,k 
	
	
	A_creux(1) = Vec(1) + (Vec(2)**3 + Vec(1)**3)/(12._rp*dx**2)
	A_creux(2) = 	       (Vec(1)**3)       /(12._rp*dx**2) !- retirer le "-"
	A_creux(3) = 	      - (Vec(2)**3)	 /(12._rp*dx**2)
	A_creux(4) = Vec(2) + (Vec(3)**3 + Vec(1)**3)/(12._rp*dx**2)
	A_creux(5) =          - (Vec(1)**3)	 /(12._rp*dx**2)
	A_creux(6) =          - (Vec(3)**3)      /(12._rp*dx**2)
	A_creux(pos_diag(2) + 1) = A_Creux(pos_diag(1) + 1)

	vois(2) = 2
	vois(3) = 3
	vois(5) = 1
	vois(6) = 4

	
	A_Creux(Pos_diag(Ne))    = Vec(Ne)    + (Vec(Ne)**3 + Vec(Ne-1)**3)  /(12._rp*dx**2)
	A_Creux(Pos_diag(Ne)+2) =         Vec(Ne)**3	 	        /(12._rp*dx**2) !-- Retirer le "-" à cause des conditions aux bords 
	A_Creux(Pos_diag(Ne)+1)  =		 -(Vec(Ne-1)**3) 			/(12._rp*dx**2)
	A_Creux(Pos_diag(Ne-1) ) = Vec(Ne-1) + (Vec(Ne)**3 + Vec(Ne-2)**3)  /(12._rp*dx**2)
	A_Creux(Pos_diag(Ne-1)+1) = 		-(Vec(Ne - 2)**3)	        /(12._rp*dx**2)	
	A_Creux(Pos_diag(Ne-1)+2) = A_creux(Pos_diag(Ne)+2)
	!
	
	vois(Pos_diag(Ne)+1) = Ne-2
	vois(Pos_diag(Ne)+2) = Ne-1
	vois(Pos_diag(Ne-1)+2) = Ne
	vois(Pos_diag(Ne-1)+1) = Ne-3
	
	
END SUBROUTINE Cond_bord_creux






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
