function[X] = SREnKF(Xprior,H,y,R,nbv)

[meas_dim, state_dim] = size(H);
m = mean(Xprior,2);
Xpert = Xprior - m;

C = Xpert*Xpert'/(nbv-1);
CHT =C*H';
S = H*CHT + R*eye(meas_dim);
K = CHT*inv(S);
KH = K*H;
KY = K*y;
mp = (eye(state_dim)- KH)*m+KY;

HX = H*Xpert/sqrt(nbv - 1);
T = R*inv(R*eye(nbv)+HX'*HX);

Tsqrt = sqrtm(T);
Zeta = Xpert*Tsqrt;
X = mp + Zeta;
end