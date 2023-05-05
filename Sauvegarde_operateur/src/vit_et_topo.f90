MODULE mod_vit_top

USE numerics

CONTAINS 

!=========================================================================
!
!	FONCTION DE CALCUL DE LA VITESSE A PARTIR DE W
!
!=========================================================================

FUNCTION vitesse(W_1,W_2)
	REAL(rp) :: W_1
	REAL(rp) :: W_2
	REAL(rp) :: vitesse
	
	IF(abs(W_1) < 1.E-10) THEN 
		vitesse = 0.0_rp
	ELSE
		vitesse = W_2/W_1
	END IF 
END FUNCTION vitesse

!=========================================================================
!
!	FONCTION DE CALCUL DE LA BATHYMETRIE DU DOMAINE 
!
!=========================================================================
SUBROUTINE topo(X,z)
	REAL(rp), INTENT(IN) :: X
	REAL(rp), INTENT(OUT) :: z

!	Ajouter plusieurs types de topographies (bosse, creux, ...)
		z = 0.0_rp	
END SUBROUTINE topo


FUNCTION sigma(W,i) !-- C'éer un tableau h_temp qui répertorie h en t-dt
	REAL(rp), DIMENSION(:,:) :: W
	INTEGER :: i
	
	REAL(rp) :: sigma 
	REAL(rp) :: ul,ur
	
	

sigma = vitesse(W(1,i),W(4,i))

END FUNCTION sigma 	


FUNCTION omega(W,i)
	REAL(rp),DIMENSION(:,:) :: W
	INTEGER :: i
	
	REAL(rp) :: omega 
	REAL(rp) :: ul,ur	

omega = vitesse(W(1,i),W(3,i))

END FUNCTION omega


FUNCTION Grad_k(u,i)
	REAL(rp), DIMENSION(:) :: u
	INTEGER :: i 
	
	REAL(rp) :: Grad_k
	
	IF(i == 1) THEN 
		Grad_k = (u(i+1) + u(i))/(2*dx)
	ELSE IF(i == Nx) THEN 
		Grad_k = (u(i+1)+u(i))/(2*dx)
	ELSE 
		Grad_k = (u(i+1) + u(i-1))/(2*dx)
	END IF 
	
END FUNCTION Grad_k



END MODULE mod_vit_top
