This folder contains files for solving numerically the stochastic heat equation.

Files included in the folder:
- StocHeatEq_implicitEuler.m: this program solves numerically the stochastic heat 
			      equation with an implicit Euler method and produces a 
		              trajectory in space-time in the interval [0,1], given some 
			      initial conditions

- LU_tridiag.m: this function calculates the vectors defining the matrices L and U in 
                which a tridiagonal matrix can be factorized 

- solve_Aud.m: this function solves the equation (LU)*v = d, where L and U are LU factors
	       of a tridiagonal matrix

