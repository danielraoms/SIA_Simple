module subs
	use parametros
	implicit none

contains


!FUNCTIONS*************************************************************************************************
!integral de uma função (array N-dimensional, onde N é o número de pontos da malha 1D) de uma variável pela regra do trapézio (terceira ordem)
function integral(f)  
implicit none
	integer :: j
	double precision, dimension(N) :: f
	double precision :: p1 			!média do valor da função no primeiro e último pontos da malha
	double precision :: p2 			!soma dos valores da função nos pontos internos
	double precision :: integral 		!output	

	p1 = (f(1) + f(N))/2.0d0

	p2 = 0.0d0
	do j = 2, (N-1)
		p2 = p2 + f(j)
	end do

	integral = (p1 + p2)*dr
	
	return	
end function

!função do initial guess, u_old
function u_o(x) 
	double precision :: u_o, p_dummy
	double precision, intent(in) :: x

	p_dummy = x
	u_o = 1.0d0 

	return 

end function 



!SUBROUTINES************************************************************************************************
!gradiente de uma função f num ponto dado por diferenças finitas (segunda ordem) 
!disclaimer: simetria radial: grad(f(r, theta, phi)) = (df/dr, 0, 0)
subroutine gradiente(f, j, d_r, grad_dummy)
implicit none
	integer, intent(in) :: j
	double precision, intent(in) :: f(j)
	double precision, intent(in) :: d_r
	double precision, intent (out) :: grad_dummy

	!diferença central (pontos internos da malha)
	if ((j .GT. 1) .AND. (j .LT. N)) then 
		grad_dummy =  (f(j+1) - f(j-1))/(2.0d0*d_r) 
	!simetria radial: condição de Neumann na origem (diferença adiantada) (primeiro ponto da malha)
	else if (j .EQ. 1) then
		grad_dummy = (-3.0d0*f(j) + 4.0d0*f(j+1) - f(j+2))/(2.0d0*dr)
	!diferença atrasada (último ponto da malha)
	else if (j .EQ. N) then 
		grad_dummy = (f(j-2) - 4.0d0*f(j-1) + 3.0d0*f(j))/(2.0d0*d_r) 
	end if

	return
end subroutine gradiente


!calcula a norma de Sobolev H0_1
subroutine H0_1_norm(r_dummy, dimensao_dummy, d_r, f_dummy, norm_dummy)
implicit none
double precision, dimension(N), intent(in) :: r_dummy, f_dummy
double precision, intent(in) :: d_r
integer, intent(in) :: dimensao_dummy
double precision, intent(out) :: norm_dummy
integer :: j
double precision, dimension(N) :: f_quadrado, grad_f, grad_quadrado
double precision :: p1, p2


	if (dimensao_dummy .eq. 2) then
		do j = 1, N
			call gradiente(f_dummy, j, d_r, grad_f(j))
			f_quadrado(j) = (f_dummy(j)**2.0d0)*r_dummy(j)	
			grad_quadrado(j) = (grad_f(j)**2.0d0)*r_dummy(j)
		end do

	 	!norma de f em L2
		p1 = (2.0d0*pi)*(integral(f_quadrado))

		!norma de f' em L2
		p2 = (2.0d0*pi)*(integral(grad_quadrado))

	else if (dimensao_dummy .eq. 3) then

		do j = 1, N
			call gradiente(f_dummy, j, d_r, grad_f(j))
			f_quadrado(j) = (f_dummy(j)**2.0d0)*(r_dummy(j)**2.0d0)
			grad_quadrado(j) = (grad_f(j)**2.0d0)*(r_dummy(j)**2.0d0)
		end do

	 	!norma de f em L2
		p1 = (4.0d0*pi)*(integral(f_quadrado))

		!norma de f' em L2
		p2 = (4.0d0*pi)*(integral(grad_quadrado))
	end if

	
	norm_dummy = sqrt(p2)
	

	!PS.: o espaço H0_1(omega) é o fecho do espaço das funções C^infty de suporte compacto com respeito à norma H1(omega).
	!Pela desigualdade de Poincaré para H0_1(omega), a norma de H0_1(omega) é equivalente à norma L2 do gradiente da f.
	!Por isso, suprimos o termo da f da norma de H0_1(omega).

	return
end subroutine H0_1_norm


!calcula a energia do funcional associado, J(g) - dimensao define a dimensão da bola de omega
subroutine J_functional(r_dummy, dimensao_dummy, d_r, g_dummy, energy_dummy)
implicit none
integer :: i
double precision, dimension(N), intent(in) :: r_dummy, g_dummy
double precision, intent(in) :: d_r
integer, intent(in) :: dimensao_dummy
double precision, dimension(N) :: grad_dummy
double precision, dimension(N) :: absol !=|grad|^2
double precision, dimension(N) :: integrando_A !=(|grad|^2)*r^2
double precision, dimension(N) :: integrando_B_parcial !=lambda*(w^2) - F(w)
double precision, dimension(N) :: integrando_B !=(lambda*(w^2) - F(w))*r^2
double precision :: termo_A, termo_B
double precision, intent(out) :: energy_dummy

	!calculando o gradiente da g_dummy
	do i = 1, N
		call gradiente(g_dummy, i, d_r, grad_dummy(i))
	end do

	if (dimensao_dummy .eq. 2) then
		do i = 1, N
			absol(i) = ((grad_dummy(i))**2.0d0) + (a*(g_dummy(i)**2.0d0))
			integrando_A(i) = r_dummy(i)*(absol(i)/2.0d0)
		end do

		do i = 1, N
			integrando_B_parcial(i) = -(b/(p+1.0d0))*(g_dummy(i)**(p+1.0d0))
			integrando_B(i) = r_dummy(i)*integrando_B_parcial(i)
		end do

		termo_A = (2.0d0*pi)*integral(integrando_A)
		
		termo_B = (2.0d0*pi)*integral(integrando_B)


	else if (dimensao_dummy .eq. 3) then
		do i = 1, N
			absol(i) = (grad_dummy(i))**2.0d0
			integrando_A(i) = ((r_dummy(i))**2.0d0)*(absol(i)/2.0d0)	
		end do
			
		do i = 1, N
			integrando_B_parcial(i) = ((a*(g_dummy(i)**2.0d0))/2.0d0) -  (b/(p+1.0d0))*(g_dummy(i)**(p+1.0d0))
			integrando_B(i) = (r_dummy(i)**2.0d0)*integrando_B_parcial(i)
		end do
		
		termo_A = (4.0d0*pi)*integral(integrando_A)

		termo_B = (4.0d0*pi)*integral(integrando_B)

	end if

	energy_dummy = (termo_A + termo_B)

	return
end subroutine J_functional


end module subs
