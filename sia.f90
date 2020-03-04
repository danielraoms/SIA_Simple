!Scaling Iterative Algorithm
!implementation by Daniel Raom Santiago - 02/03/2020
!problems of type -Δu + a*u - b*(u^p) = 0
!solved on a unit ball of R^3 in this case
!radial symmetry of the solution: u = u(r)

include 'parametros.f90' 	!declarations
include 'subs.f90' 		!subroutines necessárias a este main

program sia
use parametros
use subs
implicit none

write(solution, "(a,i0,a)") "solution.dat"
write(convergence, "(a,i0,a)") "convergence.dat"

allocate(r(N), u_old(N), u_new(N), grad_u_new(N), v_old(N), v_new(N), v_comp(N), w_old(N), w_new(N))

!calculando refinamento da malha
dr = (d-c)/N

!definindo pontos da malha
do j = 1, N
	r(j) = j*dr
end do

!definindo iteração inicial (initial guess, so to speak)
do j = 1, N
	u_old(j) = u_o(r(j))
	w_old(j) = u_old(j)
end do

	!inicializando contador das iterações
	k = 0	


!LOOP das iterações
DO
	k = k + 1	 
	write(*,*) "Iteração", k

	!PASSO 1
	!************************************************************************************
	!definindo beta_old
	beta_old = maxval(u_old)

	!para a primeira iteração, defina v_old
	if (k .eq. 1) then
		!definindo v_old
		v_old(:) = (u_old(:))/beta_old
	end if

	!PASSO 2
	!************************************************************************************	
	!C.C. de Dirichlet na fronteira: w(0) = 0
	w_new(N) = 0.0d0
	w_old(N) = 0.0d0


	!SOR step - resolvendo a EDP para achar w_new 
	do i = 1, 10000000


		!C.C. de Neumann na origem: w'(0) = 0 (forward finite difference, second order)
		w_new(1) = (4.0d0*w_old(2) - w_old(3))/3.0d0

		!pontos centrais
		do j = 2, N-1

			!termo auxiliar em comum ao problema p/ N=2 e N=3
			aux_1 = 1.0d0/(2.0d0 + a*(dr**2.0d0))

			if (dimensao .eq. 2) then
				aux_2 = -aux_1 + (dr/(4.0d0*r(j) + (2.0d0*a*r(j)*(dr**2.0d0))))
				aux_3 = -(aux_1 + (dr/(4.0d0*r(j) + (2.0d0*a*r(j)*(dr**2.0d0)))))
				aux_4 = -(dr**2)/(2.0d0 + a*(dr**2))
			else if (dimensao .eq. 3) then
				!calculando termos auxiliares 
				aux_2 = aux_1*(-1.0d0 + (dr/r(j)))
				aux_3 = -aux_1*(1.0d0 + (dr/r(j)))
				aux_4 = -(dr**2)/(2.0d0 + a*(dr**2))
			end if
			
			!calculando w_new
			w_new(j) = (-b*(v_old(j)**p)*aux_4 - aux_2*w_new(j-1) - aux_3*w_old(j+1))
		end do

		!calculando diferença máxima entre iterações
		if (mod(i,100) .eq. 0) then
			max_sor = abs(w_new(1) - w_old(1))
			do j = 1, N
				if ((abs(w_new(j) - w_old(j))) .GT. max_sor) then
					max_sor = abs(w_new(j) - w_old(j))
				end if
			end do

			write(*,*) "maxsor", i, max_sor
			!read(*,*)

			!critério de convergência
			if (max_sor .GT. 1.0e-8) then
				CONTINUE 
			else				
				EXIT
			end if
		end if


		if (mod(i,10000000) .eq. 0) then
			write(*,*) "SOR ainda não convergiu!"
			read(*,*)
		end if

		do j = 1, N-1
			w_old(j) = w_new(j)
		end do

	end do


	!PASSO 3
	!************************************************************************************
	!calculando alpha_new
	alpha_new = 1.0d0/maxval(w_new)

	!calculando beta_new
	beta_new = (beta_old**p)/alpha_new

	!calculando v_new - condição v_new(x_0) = 1 é satisfeita (o primeiro ponto é normalizado)
	v_new(:) = (w_new(:))/maxval(w_new)
	write(*,*) maxval(w_new), w_new(1), alpha_new, beta_old, beta_new
	write(*,*) 'v', v_new(1)
	!read(*,*)
	
	!PASSO 4
	!************************************************************************************
	!calculando norma H0_1 de v_new
	v_comp(:) = v_new(:) - v_old(:)

	call H0_1_norm(r, dimensao, dr, v_comp, norma_v)

	!critério de parada: convergência da norma H0_1 de (v_new - v_old)
	if (norma_v .lt. tol) then
		!calculando u_new, a SOLUÇÃO
		!u_new(:) = beta_new*(v_new(:))
		u_new(:) = (alpha_new**(1.0d0/(p-1)))*(v_new(:))

		!calculando a norma H0_1 de u_new
		call H0_1_norm(r, dimensao, dr, u_new, norma_solution)
		
		!calculando a energia do funcional associado para u_new
		call J_functional(r, dimensao, dr, u_new, energy_solution)

		!escreve num arquivo .dat os dados da convergência da solução
		open(unit = 101, file = trim(solution), status = "unknown")
	
		!escreve na tela os resultados 
		write(*,*) "The algorithm has converged with", tol, "tolerance."
		write(*,*) "H0_1 norm of solution:", norma_solution
		write(*,*) "L^infty norm of solution:", maxval(u_new), alpha_new
		write(*,*) "Energy of the functional J(v) = ", energy_solution

		write(101,*) "The algorithm has converged with", tol, "tolerance."
		write(101,*) "H0_1 norm of solution:", norma_solution
		write(101,*) "L^infty norm of solution:", maxval(u_new), alpha_new
		write(101,*) "Energy of the functional J(v) = ", energy_solution

		close(unit = 101)
	
		!escreve num arquivo .dat a solução
		open(unit = 102, file=trim(solution), status = "unknown")
		
		do j = 1, N
			write(102,*) j, r(j), u_new(j)
		end do

		close(unit = 102)

		EXIT

	else
		continue

	end if

	write(*,*) "Norma H0_1 na iteração:", norma_v, alpha_new

	!resetting arrays para próxima iteração
	v_old(:) = v_new(:)
	w_old(:) = w_new(:)
	!beta_old = beta_new

END DO

write(*,*) "Simulation has ended successfully."


end program sia
