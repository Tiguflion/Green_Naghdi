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
	&+ Diag( vec(1:Ne-2)/(4._rp*dx**2),Ne-2,-2) + Diag( vec(3:Ne)/(4._rp*dx**2)  , Ne-2,2)


END FUNCTION DaDX


END MODULE mod_operateur 
