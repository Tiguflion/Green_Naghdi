MODULE mod_calcul

USE numerics
USE mod_cond 
USE mod_flux
USE mod_vit_top

IMPLICIT NONE 


CONTAINS 



!==================================================
!	FONCTION QUI CALCULE LE PAS DE TEMPS
!
! @Entrée : W(Nx,2), Le vecteur des valeurs de la solution approchée
!
! @Sortie : dt, Le pas de temps

!==================================================
 
SUBROUTINE calcul_pas_de_temps(W,dt) 
      REAL(rp), DIMENSION(:,:), INTENT(in)  :: W
      REAL(rp),                 INTENT(out) :: dt

      INTEGER :: i
      REAL(rp)    :: hl, hr
      REAL(rp)    :: ul, cl, ur, cr
      REAL(rp)    :: aux

      aux=0.0_rp

      do i=1,Nx-1

         hl = W(1,i)
         hr = W(1,i+1)

         ul = vitesse(W(1,i),W(2,i))
         ur = vitesse(W(1,i+1),W(2,i+1))

         cl = sqrt(g*hl)
         cr = sqrt(g*hr)

         aux = max(1.E-9_rp,max(abs(ur)+cr,abs(ul)+cl,aux))
      end do
     Lambda = aux
      dt = .5_rp*dx/aux
      dt = min(dt,Tfin - date)
	

END SUBROUTINE calcul_pas_de_temps


SUBROUTINE Calcul_Shallow_Water()
INTEGER :: i
REAL(rp), DIMENSION(4) :: SolR,solL

CALL flux_godunov_WB(W,X,Flux,Nx,g)
	compt = 0
	h_temp(:) = W(1,:)
	DO i = 2,Nx-1
		CALL second_membre_godunov(X(i-1),X(i),W(:,i-1),W(:,i),SolL)
		CALL second_membre_godunov(X(i),X(i+1),W(:,i),W(:,i+1),SolR)
		
	     W(1,i) = W(1,i) - dt/dx*(Flux(1,i)-Flux(1,i-1)) 
	     W(2,i) = W(2,i) - dt/dx*(Flux(2,i)-Flux(2,i-1)) - dt/2.0_rp*g*(SolR(2) + SolL(2))
	     W(3,i) = W(3,i) - dt/dx*(Flux(3,i)-Flux(3,i-1)) 
	     W(4,i) = W(4,i) - dt/dx*(Flux(4,i)-Flux(4,i-1))
	     
	     !Pour éviter les problème de très faible Flux qui peut faire apparaitre du NaN dans la transition sec mouillée. 
	     if(W(1,i) < 1.E-8_rp) then 
	     W(1,i) = 0._rp	     
	     end if 	   
	       
	     if(abs(W(2,i)) < 1.E-8_rp) then 
	     W(2,i) = 0._rp
	     end if 
	    
	    !Calcul de la ligne de niveau  
	   ! IF(W(1,i) == 0._rp .and. compt == 0 .and. date > compt2*Tfin/Nobs) THEN
	    ! 	compt = 1
	     !	compt2 = compt2 + 1
	    ! 	CALL topo(X(i),evo_nu(compt2))
	     !	
	     !	WRITE(34,*) date,evo_nu(compt2)-H   	
	    !END IF 

	END DO 
END SUBROUTINE Calcul_Shallow_Water


SUBROUTINE Calcul_friction()
INTEGER :: i
	CALL flux_godunov_WB_friction(W,X,Flux,Nx,g)

      DO i=2,Nx-1
         W(1,i) = W(1,i) - dt/dx*(lambda*(W(1,i) - Flux(1,i) + W(1,i) -  Flux(1,i-1)))
         W(2,i) = W(2,i) - dt/dx*(lambda*(W(2,i) - Flux(2,i) + W(2,i) -  Flux(2,i-1)))
      END DO 

END SUBROUTINE Calcul_friction


END MODULE mod_calcul 
