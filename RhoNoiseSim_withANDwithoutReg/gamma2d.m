%% Gamma Matrices in 2d case
% nu = 1...16
function G=gamma2d(nu)
sigma=cell(4,1); % Pauli Matrices
sigma{1}=[1,0;0,1];
sigma{2}=[0,1;1,0];
sigma{3}=[0, -1i;1i,0];
sigma{4}=[1,0;0,-1];

i = ceil(nu/4);
j = nu-(i-1)*4;

G = 0.5*kron(sigma{i},sigma{j}); % Normalized with 0.5!
end

