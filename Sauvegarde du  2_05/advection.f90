MODULE mod_advection


USE numerics
USE mod_vit_top


IMPLICIT NONE 

CONTAINS 

SUBROUTINE Advection(W)
	REAL(rp), DIMENSION(:,:) :: W
	INTEGER :: i
	
	CALL flux_advection()
	
	DO i = 1, Nx
		W(3,i) = W(3,i) - dt/dx*(Flux(3,i)-Flux(3,i-1)) 
        	W(4,i) = W(4,i) - dt/dx*(Flux(4,i)-Flux(4,i-1))
        END DO 

END SUBROUTINE Advection


SUBROUTINE flux_advection()

	

      INTEGER :: i
      REAL(rp)    :: hl,hr, tl,tr,  sl,sr 
      REAL(rp)    :: ul,ur, wl,wr, sigl,sigr
      REAL(rp) 	  :: Lambda


      DO i=1,Nx-1
      
         hl = W(1,i)
         hr = W(1,i+1)

         ul = vitesse(W(1,i),W(2,i))
         ur = vitesse(W(1,i+1),W(2,i+1))      
         
         tl = hl*ul*wl
         tr = hr*ur*wr

         wl = omega(W(:,:),i)
         wr = omega(W(:,:),i+1)
         
         sigl = sigma(W(:,:),i)
         sigr = sigma(W(:,:),i+1)
         
         sl = hl*sigl*ul
         sr = hr*sigr*ur
         
         lambda = max(1.E-9_rp,max(0._rp,ul,ur) + sqrt(g*max(hr,hl,0._rp)))
         
         Flux(3,i) = 0.5*(tl + tr) - 0.5*lambda*(hr*wr - hl*wl)
         Flux(4,i) = 0.5*(sl + sr) - 0.5*lambda*(hr*sigr - hl*sigl) 
	END DO 
END SUBROUTINE flux_advection

END MODULE mod_advection
