# Matlab code

C_main2D implements the SUPG method through a parameter in the C_dati struct.

In order to use CG (Continuous Galerkin):
either set CG in the field corresponding to the method or set tau = 0.

It is possible to use the theoretical estimate for the stabilization parameter,
i.e. 

tau_k = 1 / sqrt(1/(Dati.dt)^2 + 30*Dati.mu^2 /femregion.h^4 + norm(beta)^2 / femregion.h^2);


using C_main2D_theory, which calls the same C_dati function
but a different C_matrix function that assembles the matrices by exploiting the theoretical formula for tau

main_ANN implements the creation and the training of the Neural Network,
calling Output1RegressionLayer, which implements the final layer of 
the Network, in which the FE Library is called in order to compute the loss

The analysis is carried out in comparison.m