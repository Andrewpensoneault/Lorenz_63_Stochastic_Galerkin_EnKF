%% Driver Script for the PCE-EnKF Lorenz 63' system
% Based on the paper "A generalized polynomial chaos based ensemble Kalman filter
% with high accuracy" from Xiu et al.
% Name: Andrew Pensoneault
% Date: 2/3/2018
close all;
clc;
%% Relevant Model Parameters
v0 = [1.508870, -1.531271, 25.46091]';    % True initial condition of the Lorenz 63 model
sigma = 10;
beta = 8/3;
rho = 28;
H = eye(1,3);                             % Measurement operator

%% Relevant parameters for DA
Q = 10;                                   % Variance of the model error
R = 1;                                   % Variance of the data error
Q_init = 100;                               % Variance of the initial condition
n_cycle = 100;                            % Number of data assimilation cycles
[meas_dim , state_dim] = size(H);         % Measurement dim and state dim
Dt = .1;                                  % Frequency of measurement
nbv = 1000;                              % Ensemble number
test_which = 1;                           % 0 - EnKF / 1 - PCE / 2 - Both / 3 - Both only plot comparison and leave it open
folder = './Figures';
%% PCE Function
gen_func = @(X1,X2,X3,x,y,z) [sum(sum(sum(hermite_matrix(x,y,z).*X1)));sum(sum(sum(hermite_matrix(x,y,z).*X2)));sum(sum(sum(hermite_matrix(x,y,z).*X3)))];
%% Initialize all containers
y = zeros(meas_dim,n_cycle);
y2 = zeros(meas_dim,n_cycle);
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

for i = 1:n_cycle     
    true_tmp = ode45(@(t,x) lorenz_63(t,x,10,8/3,28),[0,Dt],true(:,i));
    true(:,i+1) = true_tmp.y(:,end);
    true2(:,i+1) = true(:,i+1);
end

if (test_which == 0) || (test_which == 2)
    %% Perform standard SREnKF
    rng(1)
    tic();
    for i=1:n_cycle
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
end
if (test_which == 1) || (test_which == 2)
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
        [X1f,X2f,X3f] = PCE_Lorenz(X1,X2,X3,Dt,sigma,beta,rho);
        y2(:,i) = H*true2(:,i+1) + sqrt(R)*randn(meas_dim,1);
        [ma2(:,i+1),C2(:,:,i+1),X1,X2,X3] = PCEEnKF(X1f,X2f,X3f,gen_func,H,y2(:,i),R,Q,nbv);
    end
    toc()
end
%% Create Plots
plot_EnKF(folder,test_which,ma,C,true,ma2,C2,true2,Dt,n_cycle,nbv,Q,R,Q_init)
