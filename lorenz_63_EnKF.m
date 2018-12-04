%% Driver Script for the PCE-EnKF Lorenz 63' system
% Based on the paper "A generalized polynomial chaos based ensemble Kalman filter
% with high accuracy" from Xiu et al. 
% Name: Andrew Pensoneault
% Date: 2/3/2018
close all;
clc;
%% Relevant parameters for DA 
rng(1)
Q = 50;                                    % Variance of the model error
R = 10;                                    % Variance of the data error
Q_init = 100;                             % Variance of the initial condition
v0 = [1.508870, -1.531271, 25.46091]';    % True initial condition of the Lorenz 63 model
nbv = 1000;                                % Number of ensemble members
n_cycle = 40;                            % Number of data assimilation cycles
H = eye(1,3);                             % Measurement operator
[meas_dim , state_dim] = size(H);         % Measurement dim and state dim
Dt = .05;                                 % Frequency of measurement    
sigma = 10;
beta = 8/3;
rho = 28;
gen_func = @(X1,X2,X3,x,y,z) [sum(sum(sum(hermite_matrix(x,y,z).*X1)));sum(sum(sum(hermite_matrix(x,y,z).*X2)));sum(sum(sum(hermite_matrix(x,y,z).*X3)))];
%% Initialize all containers
y = zeros(meas_dim,n_cycle);
true = zeros(state_dim,n_cycle+1);
true2 = zeros(state_dim,n_cycle+1);
ma = zeros(state_dim,n_cycle+1);
ma2 = zeros(state_dim,n_cycle+1);
C = zeros(state_dim,state_dim,n_cycle+1);
C2 = zeros(state_dim,state_dim,n_cycle+1);
true(:,1) = v0;
true2(:,1) = v0;
ma(:,1) = v0;
ma2(:,1) = v0;
C(:,:,1) = Q_init*eye(state_dim);
C2(:,:,1) = Q_init*eye(state_dim);
Xprior = repmat(v0,1,nbv) + sqrt(Q_init)*randn(state_dim,nbv);
X = Xprior;
%% Perform standard SREnKF
tic();
for i=1:n_cycle
    true_tmp = ode45(@(t,x) lorenz_63(t,x,10,8/3,28),[0,Dt],true(:,i));
    true(:,i+1) = true_tmp.y(:,end);
    for j=1:nbv
        samp_tmp = ode45(@(t,x) lorenz_63(t,x,10,8/3,28),[0,Dt],X(:,j));
        Xprior(:,j) = samp_tmp.y(:,end)+sqrt(Q_init)*randn(state_dim,1);
    end
    y(:,i) = H*true(:,i+1) + sqrt(R)*randn(meas_dim,1);
    X = SREnKF(Xprior,H,y(:,i),R,nbv);
    ma(:,i+1) = mean(X,2);
    C(:,:,i+1) = (X-ma(:,i+1))*(X-ma(:,i+1))'/(nbv-1);
end
toc()

%% Perform PCE-SREnKF
rng(1)
Xprior = repmat(v0,1,nbv) + sqrt(Q_init)*randn(state_dim,nbv);
X = Xprior;
tic();
X1 = zeros(3,3,3);
X2 = zeros(3,3,3);
X3 = zeros(3,3,3);
X1(1,1,1) = v0(1);
X2(1,1,1) = v0(2);
X3(1,1,1) = v0(3);
X1(2,1,1) = sqrt(Q_init);
X2(1,2,1) = sqrt(Q_init);
X3(1,1,2) = sqrt(Q_init);
for i=1:n_cycle
    true_tmp = ode45(@(t,x) lorenz_63(t,x,10,8/3,28),[0,Dt],true2(:,i));
    true2(:,i+1) = true_tmp.y(:,end);
    [X1f,X2f,X3f] = PCE_Lorenz(X1,X2,X3,Dt,sigma,beta,rho);
    y2(:,i) = H*true2(:,i+1) + sqrt(R)*randn(meas_dim,1);
    [ma2(:,i+1),C2(:,:,i+1),X1,X2,X3] = SREnKF2(X1f,X2f,X3f,gen_func,H,y2(:,i),R,nbv);
end
toc()
%% Comparison plots SREnKF
figure()
hold on
plot3(ma(1,:),ma(2,:),ma(3,:))
plot3(true(1,:), true(2,:), true(3,:), 'o')
title({['Plot of the SEnKF mean of the Lorenz 63 model with '] ['nbv=' num2str(nbv) ' and Q=' num2str(Q), ' and R=' num2str(R)]})
legend('Assimilated solution','Truth')

figure()
hold on
errorbar(0:Dt:Dt*n_cycle,ma(1,:),squeeze(sqrt(C(1,1,:))))
plot(0:Dt:Dt*n_cycle,true(1,:),'o')
title({['Plot of the first dimension of the SEnKF mean of the Lorenz 63 model with '] ['nbv=' num2str(nbv) ' and Q=' num2str(Q), ' and R=' num2str(R)]})
legend('Assimilated solution','Truth')

figure()
hold on
errorbar(0:Dt:Dt*n_cycle,ma(2,:),squeeze(sqrt(C(2,2,:))))
plot(0:Dt:Dt*n_cycle,true(2,:),'o')
title({['Plot of the second dimension of the SEnKF mean of the Lorenz 63 model with '] ['nbv=' num2str(nbv) ' and Q=' num2str(Q), ' and R=' num2str(R)]})
legend('Assimilated solution','Truth')

figure()
hold on
errorbar(0:Dt:Dt*n_cycle,ma(3,:),squeeze(sqrt(C(3,3,:))))
plot(0:Dt:Dt*n_cycle,true(3,:),'o')
title({['Plot of the third dimension of the SEnKF mean of the Lorenz 63 model with '] ['nbv=' num2str(nbv) ' and Q=' num2str(Q), ' and R=' num2str(R)]})
legend('Assimilated solution','Truth')

%% Comparison plots PCE-SREnKF
figure()
hold on
plot3(ma2(1,:),ma2(2,:),ma2(3,:))
plot3(true2(1,:), true2(2,:), true2(3,:), 'o')
title({['Plot of the PCE-SEnKF mean of the Lorenz 63 model with '] ['nbv=' num2str(nbv) ' and Q=' num2str(Q), ' and R=' num2str(R)]})
legend('Assimilated solution','Truth')

figure()
hold on
errorbar(0:Dt:Dt*n_cycle,ma2(1,:),squeeze(sqrt(C2(1,1,:))))
plot(0:Dt:Dt*n_cycle,true2(1,:),'o')
title({['Plot of the first dimension of the PCE-SEnKF mean of the Lorenz 63 model with '] ['nbv=' num2str(nbv) ' and Q=' num2str(Q), ' and R=' num2str(R)]})
legend('Assimilated solution','Truth')

figure()
hold on
errorbar(0:Dt:Dt*n_cycle,ma2(2,:),squeeze(sqrt(C2(2,2,:))))
plot(0:Dt:Dt*n_cycle,true2(2,:),'o')
title({['Plot of the second dimension of the PCE-SEnKF mean of the Lorenz 63 model with '] ['nbv=' num2str(nbv) ' and Q=' num2str(Q), ' and R=' num2str(R)]})
legend('Assimilated solution','Truth')

figure()
hold on
errorbar(0:Dt:Dt*n_cycle,ma2(3,:),squeeze(sqrt(C2(3,3,:))))
plot(0:Dt:Dt*n_cycle,true2(3,:),'o')
title({['Plot of the third dimension of the PCE-SEnKF mean of the Lorenz 63 model with '] ['nbv=' num2str(nbv) ' and Q=' num2str(Q), ' and R=' num2str(R)]})
legend('Assimilated solution','Truth')