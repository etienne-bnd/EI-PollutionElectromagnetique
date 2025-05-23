# -*- coding: utf-8 -*-


# Python packages
import matplotlib.pyplot
import numpy
import os


# MRG packages
import _env
import preprocessing
import processing
import postprocessing
#import solutions

# not use
def your_optimization_procedure(domain_omega, spacestep, omega, f, f_dir, f_neu, f_rob,
                           beta_pde, alpha_pde, alpha_dir, beta_neu, beta_rob, alpha_rob,
                           Alpha, mu, chi, V_obj):
    """This function return the optimized density.

    Parameter:
        cf solvehelmholtz's remarks
        Alpha: complex, it corresponds to the absorbtion coefficient;
        mu: float, it is the initial step of the gradient's descent;
        V_obj: float, it characterizes the volume constraint on the density chi.
    """

    k = 0
    (M, N) = numpy.shape(domain_omega)
    numb_iter = 100
    energy = numpy.zeros((numb_iter+1, 1), dtype=numpy.float64)
    while k < numb_iter and mu > 10**(-5):
        print('---- iteration number = ', k)
        print('1. computing solution of Helmholtz problem, i.e., u')
        print('2. computing solution of adjoint problem, i.e., p')
        print('3. computing objective function, i.e., energy')
        print('4. computing parametric gradient')
        while ene >= energy[k] and mu > 10 ** -5:
            print('    a. computing gradient descent')
            print('    b. computing projected gradient')
            print('    c. computing solution of Helmholtz problem, i.e., u')
            print('    d. computing objective function, i.e., energy (E)')
            ene = compute_objective_function(domain_omega, u, spacestep)
            if bool_a:
                # The step is increased if the energy decreased
                mu = mu * 1.1
            else:
                # The step is decreased is the energy increased
                mu = mu / 2
        k += 1

    print('end. computing solution of Helmholtz problem, i.e., u')

    return chi, energy, u, grad

# not used
def your_compute_objective_function(domain_omega, u, spacestep):
    """
    This function compute the objective function:
    J(u,domain_omega)= \int_{domain_omega}||u||^2 

    Parameter:
        domain_omega: Matrix (NxP), it defines the domain and the shape of the
        Robin frontier;
        u: Matrix (NxP), it is the solution of the Helmholtz problem, we are
        computing its energy;
        spacestep: float, it corresponds to the step used to solve the Helmholtz
        equation.
    """

    energy = 0.0

    return energy

# used

def compute_projected(chi, domain, V_obj):
    """This function performs the projection of $\chi^n - mu*grad

    To perform the optimization, we use a projected gradient algorithm. This
    function caracterizes the projection of chi onto the admissible space
    (the space of $L^{infty}$ function which volume is equal to $V_{obj}$ and whose
    values are located between 0 and 1).

    :param chi: density matrix
    :param domain: domain of definition of the equations
    :param V_obj: characterizes the volume constraint
    :type chi: numpy.array((M,N), dtype=float64)
    :type domain: numpy.array((M,N), dtype=complex128)
    :type float: float
    :return:
    :rtype:
    """

    (M, N) = numpy.shape(domain)
    S = 0
    for i in range(M):
        for j in range(N):
            if domain[i, j] == _env.NODE_ROBIN:
                S = S + 1

    B = chi.copy()
    l = 0
    chi =preprocessing.set2zero(chi, domain)

    V = numpy.sum(numpy.sum(chi)) / S #c'est le beta
    debut = -numpy.max(chi)
    fin = numpy.max(chi)
    ecart = fin - debut
    # We use dichotomy to find a constant such that chi^{n+1}=max(0,min(chi^{n}+l,1)) is an element of the admissible space
    while ecart > 10 ** -4: # on cherche à obtenir un beta qui soit égal au beta objectif
        # calcul du milieu
        l = (debut + fin) / 2
        for i in range(M):
            for j in range(N):
                chi[i, j] = numpy.maximum(0, numpy.minimum(B[i, j] + l, 1))
        chi = preprocessing.set2zero(chi, domain)
        V = sum(sum(chi)) / S
        if V > V_obj:
            fin = l
        else:
            debut = l
        ecart = fin - debut
        # print('le volume est', V, 'le volume objectif est', V_obj)

    return chi

def solve_robin(alpha, k, N, chi, solve_type="adjoint", u=None): # N nombre d'élément de la discrétisation / on peut préciser le type pour savoir quel type de problème on veut résoudre le normal one et le adjoint one
    if solve_type == "adjoint" and u is None:
        raise ValueError("Parameter 'u' is required for solve_type='adjoint'")
# number of points along x-axis
    M = 2 * N  # number of points along y-axis
    level = 0 # level of the fractal
    spacestep = 1.0 / N  # mesh size
    wavenumber = k
    beta_pde, alpha_pde, alpha_dir, beta_neu, alpha_rob, beta_rob = preprocessing._set_coefficients_of_pde(M, N)
    
    
    
    
    # -- set right hand sides of the partial differential equation
    f, f_dir, f_neu, f_rob = preprocessing._set_rhs_of_pde(M, N)
    f_dir[:, :] = 0.0 # c'est les murs réfléchissant 
    f_dir[0, 0:N] = 1.0 # c'est la source
    ## f_neu, f_rob are 0 matrix
    
    if solve_type=="adjoint":
        row_s=round(3*N/8) # on défini le domaine d'intéret dans ces deux lignes
        row_e=round(5*N/8)
        f=-2*numpy.conjugate(u)
        f[row_s:row_e,row_s:row_e]=-4*numpy.conjugate(u[row_s:row_e,row_s:row_e])
    # -- set geometry of domain
    domain_omega, x, y, _, _ = preprocessing._set_geometry_of_domain(M, N, level)



    # -- initialize
    alpha_rob[:, :] = - wavenumber * 1j

    # -- define absorbing material
    Alpha = alpha
    alpha_rob = Alpha * chi

    # -- set parameters for optimization
    S = 0  # surface of the fractal
    for i in range(0, M):
        for j in range(0, N):
            if domain_omega[i, j] == _env.NODE_ROBIN:
                S += 1
    # ----------------------------------------------------------------------
    # -- Do not modify this cell, these are the values that you will be assessed against.
    # ----------------------------------------------------------------------
    # -- compute finite difference solution
    u = processing.solve_helmholtz(domain_omega, spacestep, wavenumber, f, f_dir, f_neu, f_rob,
                        beta_pde, alpha_pde, alpha_dir, beta_neu, beta_rob, alpha_rob)

    return u
def cal_energy(sol,area):
    sol=sol.reshape(-1,1)
    power2 = numpy.real(sol)**2
    nnodes=len(power2)
    energy=numpy.sum(power2)*area/nnodes
    return energy
def BelongsInteriorDomain(node):
	if (node < 0):
		return 1
	if node == 3:
# 		print("Robin")
		return 2
	else:
		return 0


def compute_gradient_descent(chi, grad, domain, mu):
	"""This function makes the gradient descent.
	This function has to be used before the 'Projected' function that will project
	the new element onto the admissible space.
	:param chi: density of absorption define everywhere in the domain
	:param grad: parametric gradient associated to the problem
	:param domain: domain of definition of the equations
	:param mu: step of the descent
	:type chi: numpy.array((M,N), dtype=float64
	:type grad: numpy.array((M,N), dtype=float64)
	:type domain: numpy.array((M,N), dtype=int64)
	:type mu: float
	:return chi:
	:rtype chi: numpy.array((M,N), dtype=float64

	.. warnings also: It is important that the conditions be expressed with an "if",
			not with an "elif", as some points are neighbours to multiple points
			of the Robin frontier.
	"""

	(M, N) = numpy.shape(domain)
	# for i in range(0, M):
	# 	for j in range(0, N):
	# 		if domain_omega[i, j] != _env.NODE_ROBIN:
	# 			chi[i, j] = chi[i, j] - mu * grad[i, j]
	# # for i in range(0, M):
	# 	for j in range(0, N):
	# 		if preprocessing.is_on_boundary(domain[i , j]) == 'BOUNDARY':
	# 			chi[i,j] = chi[i,j] - mu*grad[i,j]
	# print(domain,'jesuisla')
	#chi[50,:] = chi[50,:] - mu*grad[50,:]
	for i in range(1, M - 1):
		for j in range(1, N - 1):
			#print(i,j)
			#chi[i,j] = chi[i,j] - mu * grad[i,j]
			a = BelongsInteriorDomain(domain[i + 1, j])
			b = BelongsInteriorDomain(domain[i - 1, j])
			c = BelongsInteriorDomain(domain[i, j + 1])
			d = BelongsInteriorDomain(domain[i, j - 1])
			if a == 2:
# 				print(i+1,j, "-----", "i+1,j")
				chi[i + 1, j] = chi[i + 1, j] - mu * grad[i, j]
			if b == 2:
# 				print(i - 1, j, "-----", "i - 1, j")
				chi[i - 1, j] = chi[i - 1, j] - mu * grad[i, j]
			if c == 2:
# 				print(i, j + 1, "-----", "i , j + 1")
				chi[i, j + 1] = chi[i, j + 1] - mu * grad[i, j]
			if d == 2:
# 				print(i, j - 1, "-----", "i , j - 1")
				chi[i, j - 1] = chi[i, j - 1] - mu * grad[i, j]

	return chi
if __name__ == '__main__':

    N = 50  # number of points along x-axis
    M = 2 * N  # number of points along y-axis
    level = 0 # level of the fractal
    spacestep = 1.0 / N  # mesh size la taille d'un cube de la grille


    # on défini les différents paramètres
    wavenumber = 23*numpy.pi 
    alpha=17.23-94.01j
    domain_omega, x, y, _, _ = preprocessing._set_geometry_of_domain(M, N, level) # on défini le domaine en fonction des paramètres que l'on a initialement rentré
    V_obj=0.4 # le coefficient beta qui défini combien on a en proportion de matériau 
    lr=0.0001 # l'erreur qu'on se permet

    chi = preprocessing._set_chi(M, N, x, y) # définir un chi de départ qui n'est pas encore optimal
    chi = preprocessing.set2zero(chi, domain_omega)
    for i in range(1000):
        u1=solve_robin(alpha, wavenumber, N, chi, solve_type="normal", u=None) # on résout le système normal
        energy=cal_energy(u1,4)
        print(energy)
        q1=solve_robin(alpha, wavenumber, N, chi, solve_type="adjoint", u=u1) # on résout le système adjoint
        dJ=numpy.real(alpha*u1*q1) # on trouve ainsi la dérivé de J
        chi=compute_gradient_descent(chi,dJ, domain_omega, lr) # on récupère une nouvelle valeur de chi
        # chi[N,:]=chi[N,:]+lr*dJ[N-1,:]
        chi=compute_projected(chi, domain_omega, V_obj) # on la projette sur Uab* et on recommence la boucle pour améliorer tout ça
        
