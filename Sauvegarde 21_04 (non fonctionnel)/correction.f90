MODULE mod_correction

USE numerics 
USE mod_algebre
USE mod_cond
IMPLICIT NONE 

CONTAINS 

SUBROUTINE Correction(W)
	REAL(rp), DIMENSION(:,:) :: W
	REAL(rp), DIMENSION(:,:),  ALLOCATABLE :: sol
	REAL(rp), DIMENSION(:,:), ALLOCATABLE :: A
	REAL(rp),DIMENSION(:),ALLOCATABLE :: U ! Vitesse horizontale 
	REAL(rp),DIMENSION(:),ALLOCATABLE :: omega ! Vitesse verticale
	REAL(rp),DIMENSION(:),ALLOCATABLE :: sigma ! Ecart-type 
	REAL(rp), DIMENSION(:), ALLOCATABLE :: B
	REAL(rp), DIMENSION(:), ALLOCATABLE :: h
	REAL(rp), DIMENSION(:), ALLOCATABLE :: um2
	REAL(rp), DIMENSION(:), ALLOCATABLE :: up2
	
	REAL(rp), DIMENSION(Nx) :: vec_debug

	INTEGER i,j

	!CALL h_theo(X,date,W(1,:))

	CALL Init_Vec_Cor(sol,A,B,U,h,um2,up2,omega,sigma)


	A = Diag(h,0) + Diag(up2,-2) + Diag(up2, 2)
	

	CALL Cond_Bord_Coeff(A,B) !--- Modification des coefficients liés aux conditions aux bords

	!WRITE(*,*) A
	!READ(*,*)
	CALL gauss_seidel (Nx, A, U, B) !Faire un Choleksy
	!WRITE(*,*) 'B :'
	
	!WRITE(*,*) B
	!CALL u_theo(X,date,U)	
	!CALL w_theo(X,date,omega)
	!CALL sig_theo(X,date,sigma)
	
	!vec_debug = matmul(A,U)
	!DO i = 1,Nx
	!	IF(vec_debug(i) /= B(i)) then 
	!		WRITE(*,*) i, B(i) - vec_debug(i), B(i), vec_debug(i)
	!	END IF 
	!END DO 
	
	
	CALL calcul_omega_sigma(U,omega,sigma)

	
	W(:,:) = assemble_correction(U,omega,sigma,Nx)
	
	

END SUBROUTINE Correction 



SUBROUTINE Init_Vec_Cor(sol,A,B,U,h,um2,up2,omega,sigma)
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT):: sol
	REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: A
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: U
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: B
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: h
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: um2
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: up2
	
	
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: omega
	REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: sigma
	INTEGER :: i 
	
	REAL :: cl, cr
	
! Faire une condition avec une variable choix qui vérifie si les variables sont allouées


	ALLOCATE(sol(4,Nx))
	ALLOCATE(A(Nx,Nx))
	ALLOCATE(B(Nx),h(Nx),U(Nx))
	ALLOCATE(um2(Nx),up2(nx))
	ALLOCATE(omega(Nx),sigma(Nx))

	sol(:,:) = 0._rp
	A(:,:) = 0._rp
	B(:) = 0._rp
	h(:) = 0._rp
	um2(:) = 0._rp
	up2(:) = 0._rp
	omega(:) = 0._rp
	sigma(:) = 0._rp
	U(:) = 0._rp
	
	DO i = 2,Nx-1
		h(i) = W(1,i) + 1._rp/(12._rp*dx**2)*(W(1,i+1)**3 + W(1,i-1)**3)
		!um2(i) =-1._rp/(12._rp*dx**2)*(W(1,i+1)**3) !-- i+1 car on regarde l'indice i-1 de deux lignes après (soit i+2)
		up2(i) =-1._rp/(12._rp*dx**2)*(W(1,i+1)**3)
		cl = W(1,i-1)*(sqrt(3._rp)*W(3,i-1) + W(4,i-1))
		cr = W(1,i+1)*(sqrt(3._rp)*W(3,i+1) + W(4,i+1))
		B(i) = W(2,i) + (cr-cl)/(8._rp*sqrt(3._rp)*dx)
	END DO 
	

END SUBROUTINE Init_Vec_Cor



SUBROUTINE Cond_bord_Coeff(A,B) 
	REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: A
	REAL(rp), DIMENSION(:), INTENT(INOUT) :: B
	REAL(rp) :: c1,c2,c3,cN2,cN1,cN
	INTEGER :: i,j
	
	!=========================================================================
	!
	!			CONDITIONS AU BORD GAUCHE 
	!
	!=========================================================================
	
	A(1,1) = W(1,1) + (W(1,2)**3 + W(1,1)**3)/(12._rp*dx**2)
	A(1,2) = 	       (W(1,1)**3)      /(12._rp*dx**2) !- retirer le "-"
	A(1,3) = 	      - (W(1,2)**3)	 /(12._rp*dx**2)
	A(2,2) = W(1,2) + (W(1,3)**3 + W(1,1)**3)/(12._rp*dx**2)
	A(2,4) = 	       -(W(1,3)**3)	 /(12._rp*dx**2)
	
	WRITE(*,*) A(1,1), A(1,2),A(1,3)
	DO i = 1,2
		DO j = i+1,i+2
			A(j,i) = A(i,j)
		END DO 
	END DO 
	

	c1 = W(1,1)*(sqrt(3._rp)*W(3,1) + W(4,1))
	c2 = W(1,2)*(sqrt(3._rp)*W(3,2) + W(4,2))
	c3 = W(1,3)*(sqrt(3._rp)*W(3,3) + W(4,3))
	
	B(1) = W(2,1) + (c2 - c1)/(4._rp*sqrt(3._rp)*dx)
	B(2) = W(2,2) + (c3 - c1)/(4._rp*sqrt(3._rp)*dx)

	!=========================================================================
	!
	!			CONDITIONS AU BORD GAUCHE 
	!
	!=========================================================================
	
	A(Nx,Nx)    = -(W(1,Nx)    + (W(1,Nx)**3 + W(1,Nx-1)**3) /(12._rp*dx**2))
	A(Nx,Nx-1)  = 		   (W(1,Nx)**3)	 	 /(12._rp*dx**2) !-- Retirer le "-" à cause des conditions aux bords 
	A(Nx,Nx-2)  = 		  - (W(1,Nx-1)**3)	         /(12._rp*dx**2)
	A(Nx-1,Nx-1) = W(1,Nx-1) + (W(1,Nx)**3 + W(1,Nx-2)**3)   /(12._rp*dx**2)
	A(Nx-1,Nx-3) = 		 - (W(1,Nx - 2)**3)	         /(12._rp*dx**2)	
	
	DO i = Nx-1,Nx
		DO j = i-2,i-1
			A(j,i) = A(i,j)
		END DO 
	END DO 

	cN =  W(1,Nx)  *(sqrt(3._rp)*W(3,Nx  ) + W(4,Nx  ))
	cN1 = W(1,Nx-1)*(sqrt(3._rp)*W(3,Nx-1) + W(4,Nx-1))
	cN2 = W(1,Nx-2)*(sqrt(3._rp)*W(3,Nx-2) + W(4,Nx-2))


	
	B(Nx)   = W(2,Nx)     + (cN - cN1)/(4._rp*sqrt(3._rp)*dx)
	B(Nx-1) = W(2,Nx - 1) + (cN - cN2)/(4._rp*sqrt(3._rp)*dx)

END SUBROUTINE Cond_bord_Coeff


SUBROUTINE calcul_omega_sigma(U,omega,sigma)
REAL(rp), DIMENSION(Nx), INTENT(IN) :: U
REAL(rp), DIMENSION(Nx), INTENT(INOUT) :: omega
REAL(rp), DIMENSION(Nx), INTENT(INOUT) :: sigma

INTEGER :: i

DO i = 2,Nx-1
	omega(i) = -W(1,i)*(U(i+1) - U(i-1))/(4._rp*dx) !-- Gradiant de bathymétrie ici, négligé car on est sur fond plat
	sigma(i) = -W(1,i)*(U(i+1) - U(i-1))/(4._rp*dx*sqrt(3._rp))
END DO 

omega(1) = -W(1,1)*(U(2) + U(1))/(4._rp*dx)
sigma(1) = -W(1,1)*(U(2) + U(1))/(4._rp*dx*sqrt(3._rp))

omega(Nx) = W(1,Nx)*(U(Nx) + U(Nx-1))/(2._rp*dx)
sigma(Nx) = W(1,Nx)*(U(Nx) + U(Nx-1))/(2._rp*dx*sqrt(3._rp))


END SUBROUTINE calcul_omega_sigma


FUNCTION Diag(vec,k) 
	REAL(rp), DIMENSION(Nx) :: vec
	INTEGER :: k
	REAL(rp), DIMENSION(Nx,Nx) :: Diag
	INTEGER :: i
	
	Diag(:,:) = 0._rp
	IF(k <= 0) THEN
		DO i = 1, Nx+k
			Diag(i,i-k) = vec(i)
		END DO 
	ELSE IF(k > 0) THEN
		DO i = 1, Nx - k
			Diag(i+k,i) = vec(i)
		END DO 
	END IF 

END FUNCTION Diag
	
	
FUNCTION assemble_correction(U,omega,sigma,Nx)
	INTEGER :: Nx
	REAL(rp),DIMENSION(Nx) :: U
	REAL(rp),DIMENSION(Nx) :: omega
	REAL(rp),DIMENSION(Nx) :: sigma
	
	REAL(rp), DIMENSION(4,Nx) :: assemble_correction
	
	REAL(rp), DIMENSION(Nx) :: h
	
	h(:) = W(1,:)
	
	
	assemble_correction(1,:) = h(:)
	assemble_correction(2,:) = h(:)*U(:)
	assemble_correction(3,:) = h(:)*omega(:)
	assemble_correction(4,:) = h(:)*sigma(:)
	
END FUNCTION assemble_correction
	




	
		
	
	
END MODULE mod_correction
