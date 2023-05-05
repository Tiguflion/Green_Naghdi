MODULE mod_cond 

USE numerics
USE mod_vit_top


IMPLICIT NONE 

CONTAINS 


!===========================================================================
!
!	FONCTION  DES CONDITIONS INITIAL DU PROBLÈME
!
! @Entrées : X(Nx), Le vecteur de discretisation du domaine 
!
! @Sorties : W(4,Nx), Le Vecteur des solutions approchées
!			 H, La hauteur d'eau + fond dans le domaine (h_0 dans l'article)
!
!===========================================================================

SUBROUTINE Cond_init(W,X,H)

	REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: W
	REAL(rp), DIMENSION(:), INTENT(IN) :: X
	REAL(rp), INTENT(INOUT) :: H
	REAL(rp) :: top !penser à mettre la topologie dans un tableau 
	INTEGER :: i
	
	
	CALL h_zer(W,X,H)
	CALL hu_zer(W,X,H)
	CALL hw_zer(W,X,H)
	CALL sig_zer(W,X,H)
	
	h_temp(:) = W(1,:)
	

END SUBROUTINE Cond_init		

	
SUBROUTINE h_zer(W,X,H)
	REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: W
	REAL(rp), DIMENSION(:), INTENT(IN) :: X
	REAL(rp), INTENT(INOUT) :: H
	INTEGER :: i
	REAL(rp) :: k,eps
	eps = 0.05_rp
!DO i = 1,Nx
!	W(1,i) = eps*(1._rp + exp(-(X(i)**2)*100._rp))
!END DO 

CALL h_theo(X,0._rp,W(1,:))
END SUBROUTINE h_zer 




SUBROUTINE hu_zer(W,X,H)
	REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: W
	REAL(rp), DIMENSION(:), INTENT(IN) :: X
	REAL(rp), INTENT(INOUT) :: H
	REAL(rp), DIMENSION(Nx) :: u_zer
	INTEGER :: i
	
CALL u_theo(x,0._rp,u_zer)

W(2,:) = W(1,:)*u_zer(:)

END SUBROUTINE hu_zer 




SUBROUTINE hw_zer(W,X,H)
	REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: W
	REAL(rp), DIMENSION(:), INTENT(IN) :: X
	REAL(rp), INTENT(INOUT) :: H
	REAL(rp), DIMENSION(Nx) :: w_zer
	
	REAL(rp) :: ul,ur
	INTEGER :: i
	
CALL w_theo(x,0._rp,w_zer)

W(3,:) = W(1,:)*w_zer(:)

END SUBROUTINE hw_zer




SUBROUTINE sig_zer(W,X,H)
	REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: W
	REAL(rp), DIMENSION(:), INTENT(IN) :: X
	REAL(rp), INTENT(INOUT) :: H
	REAL(rp), DIMENSION(Nx) :: sol
	
	REAL(rp) :: ul,ur 
	INTEGER :: i
CALL sig_theo(x,0._rp,sol)
W(4,:) = W(1,:)*sol(:)

END SUBROUTINE sig_zer


!=========================================================================
!
!	FONCTION DE CALCUL DES CONDITIONS AUX BORDS 
!
!=========================================================================

SUBROUTINE flux_godunov_WB_bord_Wall(W,X,F0, FNx,ns,g)

      IMPLICIT NONE

      INTEGER,      	        INTENT(IN)  :: ns
      REAL(rp), DIMENSION(:,:), INTENT(IN)  :: W
      REAL(rp), DIMENSION(4)  ,	INTENT(OUT) :: F0,FNx
      REAL(rp),                 INTENT(IN)  :: g
      REAL(rp), DIMENSION(:)  , INTENT(IN)  :: X

      INTEGER :: i
      REAL(rp)    :: hl, hr, pil, pir,tl,tr,  sl, sr 
      REAL(rp)    ::  ul, cl, ur, cr, wr, wl,sigl, sigr
      REAL(rp)    :: zl,zr, lambda 

    !=====================================================================
    !		Calcul du Flux au bord gauche
    !=====================================================================

	CALL topo(X(1),zl)
	CALL topo(X(1),zr)
	
         hl = W(1,1)
         hr = W(1,1)

         ul = -vitesse(W(1,1),W(2,1))
         ur = vitesse(W(1,1),W(2,1))      

         pil = hl*ul**2 + g*0.5*hl**2
         pir = hr*ur**2 + g*0.5*hr**2
         
         tl = hl*ul*wl
         tr = hr*ur*wr


         wl = omega(W(:,:),1)
         wr = omega(W(:,:),1)
         
         sigl = sigma(W(:,:),1)
         sigr = sigma(W(:,:),1)
         
         sl = hl*sigl*ul
         sr = hr*sigr*ur
         
         lambda = max(1.E-9_rp,max(0._rp,ul,ur) + sqrt(g*max(hr,hl,0._rp)))
         

         F0(1) = 0.5*(hl*ul+hr*ur) - 0.5*lambda*( (hr+zr) - (hl+zl))
         F0(2) = 0.5*(pil + pir) - 0.5*lambda*(hr*ur-hl*ul)
         F0(3) = 0.5*(tl + tr) - 0.5*lambda*(hr*wr - hl*wl)
         F0(4) = 0.5*(sl + sr) - 0.5*lambda*(hr*sigr - hl*sigl) 


    !=====================================================================
    !		Calcul du Flux au bord droit
    !=====================================================================


	CALL topo(X(ns),zl)
	CALL topo(X(ns),zr)
	
         hl = W(1,ns)
         hr = W(1,ns)

         ul = vitesse(W(1,ns),W(2,ns))
         ur = -vitesse(W(1,ns),W(2,ns))      

         pil = hl*ul**2 + g*0.5*hl**2
         pir = hr*ur**2 + g*0.5*hr**2
         
         tl = hl*ul*wl
         tr = hr*ur*wr


         wl = omega(W(:,:),ns)
         wr = omega(W(:,:),ns)
         
         sigl = sigma(W(:,:),ns)
         sigr = sigma(W(:,:),ns)
         
         sl = hl*sigl*ul
         sr = hr*sigr*ur
         
         lambda = max(1.E-9_rp,max(0._rp,ul,ur) + sqrt(g*max(hr,hl,0._rp)))
         

         FNx(1) = 0.5*(hl*ul+hr*ur) - 0.5*lambda*( (hr+zr) - (hl+zl))
         FNx(2) = 0.5*(pil + pir) - 0.5*lambda*(hr*ur-hl*ul)
         FNx(3) = 0.5*(tl + tr) - 0.5*lambda*(hr*wr - hl*wl)
         FNx(4) = 0.5*(sl + sr) - 0.5*lambda*(hr*sigr - hl*sigl) 
	

END SUBROUTINE flux_godunov_WB_bord_Wall



SUBROUTINE Cond_bord(W)
	REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: W
	INTEGER :: i
	REAL(rp), DIMENSION(Nx) :: h_sol, u_sol, w_sol, sig_sol
	REAL(rp), DIMENSION(4) :: F0,FNx
	
CALL h_theo(X,date,h_sol)
CALL u_theo(X,date,u_sol)

	!-- Conditions de Wall, ce qui implique que le flux en 0 est nul, aisni que le flux en Nx + 1
CALL flux_godunov_WB_Bord_Wall(W,X,F0,FNx,Nx,g)


	     W(1,1) = W(1,2) !W(1,1) - dt/dx*(Flux(1,1) - F0(1))
	     W(2,1) = W(2,2) !0._rp !W(2,1) - dt/dx*(Flux(2,1) - F0(2)) !-- La bathymétrie est identique aux indices 0 et 1
	     W(3,1) = W(3,2) !W(3,1) - dt/dx*(Flux(3,1) - F0(3))
	     W(4,1) = W(4,2) !- dt/dx*(Flux(4,1) - F0(4))
	     
	     W(1,Nx) = W(1,Nx-1) !W(1,Nx) - dt/dx*(FNx(1) - Flux(1,Nx))
	     W(2,Nx) = W(2,Nx-1) !0._rp !W(2,Nx) - dt/dx*(FNx(2) - Flux(2,Nx)) !-- La bathymétrie est identique aux indices Nx et Nx+1
	     W(3,Nx) = W(3,Nx-1) !W(3,Nx) - dt/dx*(FNx(3) - Flux(3,Nx))
	     W(4,Nx) = W(4,Nx-1) !W(4,Nx) - dt/dx*(FNx(4) - Flux(4,Nx))
	     

END SUBROUTINE Cond_bord



!=========================================================================
!
!	FONCTION DE CALCUL DES SOLUTIONS THEORIQUES
!
!=========================================================================


SUBROUTINE h_theo(x,t,sol)
	REAL(rp), DIMENSION(:) :: x
	REAL(rp)	       :: t
	REAL(rp), DIMENSION(:) :: sol
	INTEGER :: i	
	REAL(rp) :: k, c
	
   k = sqrt(3._rp*A/(4*(1+A)))/H
   c = sqrt(H*g*(1+A))

DO i = 1,Nx
	sol(i) =H*(1 + A/(cosh(k*(x(i) - x0 - c*t))**2))
END DO 

END SUBROUTINE h_theo

SUBROUTINE u_theo(x,t,sol)
	REAL(rp), DIMENSION(:), INTENT(IN) :: x
	REAL(rp)	      , INTENT(IN) :: t
	REAL(rp), DIMENSION(:), INTENT(OUT):: sol
	INTEGER :: i	
	REAL(rp) :: k, c
	REAL(rp), DIMENSION(Nx):: h_sol
	
CALL h_theo(x,t,h_sol)
	
DO i = 1,Nx
	sol(i) = (1 - H/h_sol(i))*sqrt((1+A)*g*H)
END DO 

END SUBROUTINE u_theo


SUBROUTINE w_theo(x,t,sol)
	REAL(rp), DIMENSION(:) :: x
	REAL(rp)	       :: t
	REAL(rp), DIMENSION(:) :: sol
	REAL(rp)	       :: ksi, k,c
	REAL(rp), DIMENSION(Nx):: h_sol
	
	INTEGER :: i 
	
	k = sqrt(3._rp*A/(4*(1+A)))/H
	c = sqrt(H*g*(1+A))
	
CALL h_theo(x,t,h_sol)

DO i = 1,Nx
	ksi = k*(x(i) - x0 - c*t)
	sol(i) = sqrt(0.75_rp*g)*(A*H)**(1.5_rp)*tanh(ksi)/(cosh(ksi)**2*h_sol(i))
END DO 


END SUBROUTINE w_theo



SUBROUTINE sig_theo(x,t,sol)
	REAL(rp), DIMENSION(:) :: x
	REAL(rp)	       :: t
	REAL(rp), DIMENSION(:) :: sol


CALL w_theo(x,t,sol)

sol(:) = sol(:)/sqrt(3._rp)


END SUBROUTINE sig_theo


END MODULE mod_cond 
