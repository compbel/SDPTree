function obj = calcObjMatrix(T,D,unrooted)

[D_T,rho] = leafDist(T);
Sigma = (rho + rho' - D_T)/2;
rho_norm = 2.^(-rho);
U = 4.^(Sigma);
if unrooted
    U(Sigma==0) = 2;
end
R = D.*U;
obj = rho_norm'*R*rho_norm;

