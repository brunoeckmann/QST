%% Quantum State Tomography
% Calculate reconstruced density matrix
function [rho_reconstr, M] = qst_linearinversion(measurementBasis,countSignal)
    % Calculate B Matrix and Inverse
    B_Matrix=zeros(16);
    for nu = 1:16
        for mu=1:16
            B_Matrix(nu,mu)=measurementBasis{nu}'*gamma2d(mu)*measurementBasis{nu};
        end
    end
    
    B_Matrix_Inv=inv(B_Matrix);
    s = size(measurementBasis);

    % Calculate M Matrices
    M=cell(s(1),1);
    for nu=1:s(1)
        M{nu}=zeros(4);
        for mu=1:s(1)
            M{nu}=M{nu}+B_Matrix_Inv(mu,nu)*gamma2d(mu);
        end
    end
    


    r=zeros(s(1),1);
    for nu=1:s(1)
        r(nu)=dot(B_Matrix_Inv(nu,:),countSignal);
    end
    
    %%Method 1 to calculate rho
    %rho=zeros(4);
    %for nu=1:16
    %   rho=rho+gamma2d(nu)*r(nu); 
    %end

    %% Method 2 to calculate rho
    rho_reconstr=zeros(4);
    for nu=1:s(1)
        rho_reconstr=rho_reconstr+M{nu}*countSignal(nu);
    end

end