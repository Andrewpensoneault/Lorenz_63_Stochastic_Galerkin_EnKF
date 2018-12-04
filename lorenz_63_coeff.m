function[y] = lorenz_63_coeff(t,x,beta,sigma,rho,wei,loc)
%RHS for ode45 solver for Lorenz 63 PCE coefficients
%% Reshape into matricies
maty1 = zeros(3,3,3);
maty2 = zeros(3,3,3);
maty3 = zeros(3,3,3);
X1 = reshape(x(1:27),3,3,3);
X2 = reshape(x(28:54),3,3,3);
X3 = reshape(x(55:81),3,3,3);
%% Evaluate integral
for i=1:3
    for j=1:3
        for k=1:3
            H = hermite_matrix(loc(i),loc(j),loc(k));
            x = repmat(sum(sum(sum(H.*X1))),3,3,3);
            y = repmat(sum(sum(sum(H.*X2))),3,3,3);
            z = repmat(sum(sum(sum(H.*X3))),3,3,3);
            maty1 = maty1 + wei(i)*wei(j)*wei(k)*sigma*H.*(y-x);
            maty2 = maty2 + wei(i)*wei(j)*wei(k)*H.*(x.*(rho-z)-y);
            maty3 = maty3 + wei(i)*wei(j)*wei(k)*H.*(x.*y-beta*z);
        end
    end
end
%% Reshape to vector
y = [maty1(:);maty2(:);maty3(:)];
end