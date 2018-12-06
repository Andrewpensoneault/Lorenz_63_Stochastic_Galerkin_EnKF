function[X] = SREnKF(Xprior,H,y,R,nbv)
[meas_dim, state_dim] = size(H);
m = mean(Xprior,2);
C = ((Xprior-m)*(Xprior-m)')/(nbv-1);
CHT =C*H';
Xpert = Xprior-mean(Xprior,2);
S = H*CHT + R*eye(meas_dim);
K = CHT/S;
Ktilde = CHT/((sqrtm(H*CHT+R*eye(meas_dim)))')/(sqrtm(H*CHT+R*eye(meas_dim))+sqrt(R)*eye(meas_dim));
mp = m + K*(y-H*m);
Xperta = Xpert - Ktilde*H*Xpert;
X = mp + Xperta;
end