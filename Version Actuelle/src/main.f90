PROGRAM main


USE numerics
USE mod_cond 
USE mod_flux
USE mod_calcul
USE mod_dat 
USE mod_correction
USE mod_advection

IMPLICIT NONE

INTEGER :: i

REAL(rp) :: temp,t1,t2

INTEGER :: compt2

!====================================================================================
!
! Mettre tout les sigma et la bathymétrie  dans un tableau pour ne pas les recalculer
!
!====================================================================================

CALL cpu_time(t1)
CALL init_var(niter,date,dt,dx)

CALL init_tab()

CALL Cond_init(W(:,:),X(:),H)

OPEN(UNIT = 34, FILE = 'niv_eau.dat',STATUS= 'UNKNOWN')

compt2 = 0

DO WHILE (date < Tfin)
	h_temp(:) = W(1,:)

	CALL calcul_pas_de_temps(W,dt)

	date = date + dt
	niter = niter+1
	

IF(alp <= 1E-6) THEN !-- Distinction des cas avec et sans friction
	CALL Calcul_Shallow_Water()
ELSE
	CALL calcul_friction()	
END IF 

	CALL Cond_bord(W) !Calcul de W(:,1) et W(:,Nx) avant la correction

	CALL Correction(W,Nx)

	
write(*,*) 'Temps :',date, ' Temps final :', Tfin 

END DO 


CALL h_theo(x,Tfin,h_temp)
!======================================================================
!
!		SORTIES DANS LE TERMINAL ET AFFICHAGE DES GRAPHIQUES 
!
!
!======================================================================


CALL sol_dat()
CALL CPU_TIME(t2)

WRITE(*,*) 'Temps de la résolution : ', t2 - t1


END PROGRAM main 
