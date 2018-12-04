function[X1,X2,X3] = PCE_Lorenz(X1,X2,X3,Dt,sigma,beta,rho)
% PCE expansion for the Stochastic Galerkin to be applied to the PCE-EnKF
% For the Lorenz 63' system
%% Initializes ODE solvers
[loc,wei] = GaussHermite(3);
loc = loc*sqrt(2);
wei = wei/sqrt(pi);
%% Creates the initial condition of the coefficients
X = [X1(:);X2(:);X3(:)];
%% Creates the unintegrated rhs of the coefficient odes
V_vec = ode45(@(t,x) lorenz_63_coeff(t,x,beta,sigma,rho,wei,loc),[0,Dt],X);
coeff = V_vec.y(:,end);
X1 = reshape((coeff(1:27)),3,3,3);
X2 = reshape((coeff(28:54)),3,3,3);
X3 = reshape((coeff(55:81)),3,3,3);
end