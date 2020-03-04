module parametros
implicit none

character(100) :: solution, convergence
integer :: i, j, k !contadores
double precision, parameter :: pi = 4.0d0*datan(1.0d0)
integer, parameter :: dimensao = 3				 !dimensão do domínio omega, que normalmente é uma bola de raio 1 centrada no zero
integer, parameter :: N = 2001					 !número de pontos da malha 1D
double precision, parameter :: a = 0.0d0, b = 2.0d0, p = 3.0d0 	 !da equação
double precision, parameter :: c = 0.0d0, d = 1.0d0	 	 !domínio
double precision, parameter :: tol = 1e-8			 !tolerância para convergência da norma H0_1
double precision, parameter :: w_sor = 1.8
double precision :: dr, max_sor
double precision :: alpha_new, beta_old, beta_new
double precision :: aux_1, aux_2, aux_3, aux_4
double precision :: norma_v, norma_solution, energy_solution
double precision, allocatable, dimension(:) :: r, u_old, u_new, grad_u_new, v_old, v_new, v_comp, w_old, w_new

end module parametros
