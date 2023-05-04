MODULE mod_dat

USE numerics
USE mod_cond

IMPLICIT NONE 


CONTAINS 


SUBROUTINE init_var(niter,date,dt,dx)

INTEGER, INTENT(INOUT) :: niter
REAL(rp), INTENT(INOUT) :: date
REAL(rp), INTENT(INOUT) :: dt 
REAL(rp), INTENT(INOUT) :: dx
 
OPEN(UNIT = 50, FILE = 'donnee', STATUS = 'OLD', ACTION = 'READ')

READ(50,*) xmin, xmax
READ(50,*) Nx
READ(50,*)Tfin
READ(50,*)Nobs
READ(50,*) Cas_test
READ(50,*) Cas_topo
READ(50,*) alp
READ(50,*) H
READ(50,*) A
READ(50,*) x0


CLOSE(50)

niter = 0
date = 0.0_rp
dt = 0.0_rp
dx = (xmax-xmin)/(Nx-1)

END SUBROUTINE init_var 


SUBROUTINE init_tab()
INTEGER :: i



ALLOCATE(W(4,Nx), W_temp(4,Nx), Flux(4,Nx), D(4,Nx), F(4,Nx))
ALLOCATE(X(Nx))
ALLOCATE(S(4,Nx))
ALLOCATE(h_temp(Nx))

W(:,:) = 0._rp
Flux(:,:) = 0._rp
D(:,:) = 0._rp
F(:,:) = 0._rp
h_temp = 0._rp

ALLOCATE(evo_nu(Nobs))

DO i = 1, Nx
	X(i) = xmin + (i-1)*(xmax-xmin)/(Nx-1) 
END DO

END SUBROUTINE init_tab




SUBROUTINE sol_dat()
REAL(rp) :: pos
INTEGER :: i
REAL(rp), DIMENSION(Nx) :: U
REAL(rp), DIMENSION(Nx) :: W_ex
REAL(rp), DIMENSION(Nx) :: sig_ex
REAL(rp) :: errmax, errL2
      write(6,*) 'fin de l''algorithme en',niter,' iterations'
      write(6,*) 'date de fin',date, Tfin


      CALL u_theo(X,Tfin,U)
      CALL w_theo(X,Tfin,W_ex)
      CALL sig_theo(X,Tfin,sig_ex)
      OPEN(unit=666, file = 'erreur',status = 'UNKNOWN')
      WRITE(666,*) "Nombez d'éléments :", Nx
      WRITE(666,*) 
      WRITE(666,*) "Variable h:"
      WRITE(666,*) "Erreur moyenne  :", sum(abs(h_temp(:)-W(1,:)))/Nx
      WRITE(666,*) "Erreur inf :", maxval(abs(h_temp(:)-W(1,:)))
      WRITE(666,*) "Variable U:"
      WRITE(666,*) "Erreur moyenne  :", sum(abs(U(:)-W(2,:)/W(1,:)))/Nx
      WRITE(666,*) "Erreur inf :", maxval(abs(U(:)-W(2,:)/W(1,:)))
      WRITE(666,*) "Variable W:"
      WRITE(666,*) "Erreur moyenne  :", sum(abs(W_ex(:)-W(3,:)/W(1,:)))/Nx
      WRITE(666,*) "Erreur inf :", maxval(abs(W_ex(:)-W(3,:)/W(1,:)))
      WRITE(666,*) "Variable sigma :"
      WRITE(666,*) "Erreur moyenne  :", sum(abs(sig_ex(:)-W(4,:)/W(1,:)))/Nx
      WRITE(666,*) "Erreur inf :", maxval(abs(sig_ex(:)-W(4,:)/W(1,:)))
      close(666)
      open(unit=21,file='solutions',status = 'UNKNOWN')
      do i=1,Nx
      CALL topo(X(i),pos)
        write(21,*) X(i),W(1,i)+pos,h_temp(i),W(2,i)/W(1,i),U(i), W(3,i)/W(1,i),W_ex(i),W(4,i)/W(1,i),sig_ex(i)
        
      enddo

	close(21)
	
END SUBROUTINE sol_dat 

END MODULE mod_dat 


