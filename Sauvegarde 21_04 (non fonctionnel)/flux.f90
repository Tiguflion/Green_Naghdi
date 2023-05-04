MODULE mod_flux 

USE numerics
USE mod_cond 
USE mod_vit_top
IMPLICIT NONE 

CONTAINS 


!======================================================================
!	FONCTION QUI CALCULE LE FLUX DE GODUNOV WELL BALANCED DANS TOUT 
!	LE DOMAINE POUR UN CAS SANS FRICTION
!
! @Entrée : W(Nx,2), Le vecteur des valeurs de la solution approchée
!			X(Nx), Le vecteur de discretisation du domaine 
!			g, La constante de gravité (g = 9.81 m.s^{-2})
!			Nx, Le nombre de volumes d'approximations
!
! @Sortie : Flux, Le flux de Godunov Well-Balanced 
!======================================================================

SUBROUTINE flux_godunov_WB(W,X,Flux,ns,g)

      IMPLICIT NONE

      INTEGER,              INTENT(IN)  :: ns
      REAL(rp), DIMENSION(:,:), INTENT(IN)  :: W
      REAL(rp), DIMENSION(:,:), INTENT(OUT) :: Flux
      REAL(rp),                 INTENT(IN)  :: g
      REAL(rp), DIMENSION(:)  , INTENT(IN)  :: X

      INTEGER :: i
      REAL(rp)    :: hl, hr, pil, pir,tl,tr,  sl, sr 
      REAL(rp)    ::  ul, cl, ur, cr, wr, wl,sigl, sigr
      REAL(rp)    :: zl,zr, lambda 

      DO i=1,ns-1

	CALL topo(X(i),zl)
	CALL topo(X(i+1),zr)
	
         hl = W(1,i)
         hr = W(1,i+1)

         ul = vitesse(W(1,i),W(2,i))
         ur = vitesse(W(1,i+1),W(2,i+1))      

         pil = hl*ul**2 + g*0.5*hl**2
         pir = hr*ur**2 + g*0.5*hr**2
         
         tl = hl*ul*wl
         tr = hr*ur*wr


         wl = omega(W(:,:),i)
         wr = omega(W(:,:),i+1)
         
         sigl = sigma(W(:,:),i)
         sigr = sigma(W(:,:),i+1)
         
         sl = hl*sigl*ul
         sr = hr*sigr*ur
         
         lambda = max(1.E-9_rp,max(0._rp,ul,ur) + sqrt(g*max(hr,hl,0._rp)))
         

         Flux(1,i) = 0.5*(hl*ul+hr*ur) - 0.5*lambda*( (hr+zr) - (hl+zl))
         Flux(2,i) = 0.5*(pil + pir) - 0.5*lambda*(hr*ur-hl*ul)
         Flux(3,i) = 0.5*(tl + tr) - 0.5*lambda*(hr*wr - hl*wl)
         Flux(4,i) = 0.5*(sl + sr) - 0.5*lambda*(hr*sigr - hl*sigl) 

         !-- Filtre pour éviter d'avoir des NaN dans le cas de vagues  
         
         IF(abs(Flux(1,i)) < 1.E-10) THEN
         	Flux(1,i) = 0.0_rp
         END IF
         IF(abs(Flux(2,i)) < 1.E-10) THEN
         	Flux(2,i) = 0.0_rp
         END IF 
      END DO 


END SUBROUTINE flux_godunov_WB


SUBROUTINE flux_godunov_WB_friction(W,X,Flux,ns,g)

      IMPLICIT NONE

      INTEGER,              INTENT(IN)  :: ns
      REAL(rp), DIMENSION(:,:), INTENT(IN)  :: W
      REAL(rp), DIMENSION(:,:), INTENT(OUT) :: Flux
      REAL(rp),                 INTENT(IN)  :: g
      REAL(rp), DIMENSION(:)  , INTENT(IN)  :: X

      INTEGER :: i
      REAL(rp)    :: hl, hr, pil, pir,  tl , tr
      REAL(rp)    :: ul, cl, wl, ur, cr, wr
 

      DO i=1,Nx-1
		!-- Essayer de voir un schémas qui ajoute de la topographie dans la résolution 
         hl = W(1,i)
         hr = W(1,i+1)

         ul = vitesse(W(1,i),W(2,i+1))
         ur = vitesse(W(1,i+1),W(2,i+1))

         wl = vitesse(W(1,i),W(3,i+1))
         wr = vitesse(W(1,i+1),W(3,i+1))

         pil = hl*ul**2 + g*0.5*hl**2
         pir = hr*ur**2 + g*0.5*hr**2
   
         Flux(1,i) = 0.5*(hl+hr) - 1./(2*lambda)*( (hr*ur) - (hl*ul))

         Flux(2,i)=0.5*( (hl*ul) + (hr*ur)) -1./(2*lambda*alp*dt)*(1-exp(-alp*dt))*((pir-pil)&
         & + alp*dx*(hl*ul+hr*ur)*(1-exp(-alp*dt))*.5)

	 !-- Filtre pour éviter d'avoir des NaN dans le cas de vagues 
	 
         IF(Flux(1,i) < 1.E-10) THEN
         	Flux(1,i) = 0.0_rp
         END IF
         IF(Flux(2,i) < 1.E-10) THEN
         	Flux(2,i) = 0.0_rp
         END IF 
         
      END DO 


END SUBROUTINE flux_godunov_WB_friction


!=========================================================================
!	FONCTION QUI CALCULE LE SECOND MEMBRE DE L'EQUATION DE SHALLOW WATER
!	POUR UN SCHEMAS DE GODUNOV WELL BALANCED SANS FRICTION
!
!@Entrées : xl, le point d'approximation gauche de l'intervalle
!			xr, le point d'approximation droit de l'intervalle
!			W_l, les valeurs d'approximation sur l'intervalle précédent
!			W_r, les valeurs d'approximation sur l'intervalle suivant 
!
!@Sorties : sol, le vecteur du second membre de l'approximation
!
!=========================================================================

SUBROUTINE second_membre_godunov(xl,xr,W_l,W_r,sol)
	REAL(rp), INTENT(IN) :: xl,xr
	REAL(rp), DIMENSION(4), INTENT(IN) :: W_l,W_r
	REAL(rp), DIMENSION(4), INTENT(OUT) :: sol
	REAL(rp) :: sol1,sol2
	CALL topo(xl,sol1)
	CALL topo(xr,sol2)
	sol(1) = 0.0_rp
	sol(2) = .5_rp*(W_l(1)+W_r(1))*(sol2-sol1)/dx
	sol(3) = 0.0_rp
	sol(4) = 0.0_rp
END SUBROUTINE second_membre_godunov


END MODULE mod_flux
