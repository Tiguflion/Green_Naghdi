MODULE mod_operateur

USE numerics


!===================================================================================================================
!
! Afin de comprendre de que font les fonctions, on utilise les accronymes suivants : 
!
! - D signifie "Dérivé"
! - a Correspond à la variable que l'on mets en indice (exemple la fonction Da(u) corresponds à la fonction qui calcule la dérivée de u)
!
! -a1,a2, ... ont le même fonctionnement que a, elles sont utilisées dans le cadre où l'on utilise plusieurs variables. 
!
!
! Dans le cas d'un problème matriciel AX = B, il apparait également :
! 
! - X, l'inconnue du problème, (il apparait généralement dans l'énoncé des fonctionq qui fabriquent A)
!
! 
! Les fonctions qui n'utilisent pas ce jeu de variable, seront commentées
!
!
! Lexique des variables : 
! - Ne, le nombre d'éléments des vecteurs 
! - dx, l'ecart entre deux points 
!
!=================================================================================================================== 



CONTAINS 










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!					PARTIE 1 : FONCTIONS SUR DES VECTEURS
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











FUNCTION Da(a,dx,Ne)   ! ATTENTION : CETTE FONCTION NE PRENDS PAS EN COMPTE LES CONDITIONS AUX BORDS
	INTEGER :: Ne
	REAL(rp):: dx
	REAL(rp), DIMENSION(Ne) :: a
	
	
	REAL(rp), DIMENSION(Ne) :: Da
	
	INTEGER :: i
	
	Da = 0._rp

!Gradiant centré
	
	Da(2:Ne-1) = (a(3:Ne) - a(1:Ne-2))/(2._rp*dx)

END FUNCTION Da










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!					PARTIE 2 : FONCTIONS SUR DES MATRICES
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










!==============================================================
! Fonction qui remplie une diagonale de matrice 
!
!Entrée : Vec, un vecteur de valeur
!	  Ne, Le nombre d'éléments de la matrice 
!	  k, la position de la diagonale (positif si le décalage est la droite et négatif à droite, et vaut 0 s'il n'y a pas de décalage)	
!==============================================================

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


FUNCTION DaDX(Vec,dx,Ne)
	INTEGER  :: Ne
	REAL(rp) :: dx ! Dans le cadre ou l'on a un écart constant entre les points d'observations
	REAL(rp), DIMENSION(Ne) :: Vec
	
	REAL(rp), DIMENSION(Ne,Ne) :: DaDX
	
	DaDX(:,:) = 0._rp
	
	DaDX(2:Ne-1,2:Ne-1) = -Diag( (vec(1:Ne-2) + vec(3:Ne))/(4._rp*dx**2) ,Ne-2,0) &
	&+ Diag( vec(3:Ne)/(4._rp*dx**2),Ne-2,-2) + Diag( vec(3:Ne)/(4._rp*dx**2)  , Ne-2,2)


END FUNCTION DaDX










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!					PARTIE 2 : FONCTIONS SUR DES MATRICES
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










!==============================================================================================================
!
!					FONCTION DE RECHERCHE DE VOISINS EN FONCTION DES OPERATEURS
!
!==============================================================================================================

SUBROUTINE Vois_DaDX(i,Ne,Nb_Vois,pos_diag)

	INTEGER 	      , INTENT(IN)    :: i
	INTEGER		      , INTENT(IN)    :: Ne
	INTEGER, DIMENSION(:) , INTENT(INOUT) :: pos_diag
	INTEGER 	      , INTENT(INOUT) :: Nb_vois

	INTEGER :: j
	
	Nb_Vois = Nb_vois+2
	DO j = i+1,Ne+1
		pos_diag(j) = pos_diag(j)+2
	END DO 
	
END SUBROUTINE Vois_DaDX






!==============================================================================================================
!
!					FONCTION D'OPERATEURS CREUX 
!
!==============================================================================================================


SUBROUTINE Diag_creux(Vec,Ne, k, Pos_diag, Vois, A_creux)
	INTEGER, INTENT(IN) :: Ne
	INTEGER, INTENT(IN) :: k
	INTEGER,  DIMENSION(:), INTENT(IN)    :: Pos_diag
	REAL(rp), DIMENSION(:), INTENT(IN)    :: vec
	INTEGER,  DIMENSION(:), INTENT(INOUT) :: Vois 
	REAL(rp), DIMENSION(:), INTENT(INOUT) :: A_creux 
	
	INTEGER :: i,j
	
	DO i = 1,Ne
		DO j = Pos_diag(i), Pos_diag(i+1) - 1
			IF(i+k == Vois(j) .or. (Vois(j) == 0 .and. i+k > 0)) then 
				Vois(j) = i+k
				A_Creux(j) = A_creux(j) + Vec(i)
				EXIT
			END IF 
		END DO 
	END DO 


END SUBROUTINE Diag_Creux






SUBROUTINE DaDX_Creux(i,Ne,Vec,pos_diag,vois,A_creux)

	INTEGER, DIMENSION(:),  INTENT(IN) :: pos_diag 
	INTEGER,	        INTENT(IN) :: i
	INTEGER, 	        INTENT(IN) :: Ne
	REAL(rp), DIMENSION(:),	INTENT(IN) :: Vec
	INTEGER, DIMENSION(:),  INTENT(INOUT) :: vois
	REAL(rp), DIMENSION(:), INTENT(INOUT) :: A_creux
	
	
	REAL(rp), DIMENSION(3) :: Coeff
	INTEGER, DIMENSION(3) :: Pos
	
	INTEGER k,j


	Coeff(:) = 0._rp
	
	Pos(1) = i
	Pos(2) = i-2
	Pos(3) = i+2
	
	Coeff(1) = -(Vec(i+1) + Vec(i-1))/(4._rp*dx**2)
	Coeff(2) = Vec(i-1)/(4._rp*dx**2)
	Coeff(3) = Vec(i+1)/(4._rp*dx**2)
	
	!A_Creux(pos_diag(i)) = A_creux(pos_diag(i)) + Coeff(1)


	!----- Shémas type de remplissage creux
	k = 1
	DO j = Pos_diag(i), Pos_diag(i+1)-1
		if((vois(j) == Pos(k) .or. vois(j) == 0) .and. k <= 3) then 
			if(vois(j) == 0) then 
				A_Creux(j) = Coeff(k)
			else 
				A_creux(j) = A_creux(j) + Coeff(k)
			end if 
			vois(j) = pos(k)
			k = k+1
		end if 
	END DO 
 

END SUBROUTINE DaDX_Creux

END MODULE mod_operateur 
