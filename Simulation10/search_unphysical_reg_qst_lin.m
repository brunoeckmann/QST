%% Artificial signalling behaviour
p_sig=zeros(2,2,3,3);

p_sig(:,:,1,1) = [0.5 0;0 0.5];
p_sig(:,:,1,2) = [0.5 0;0 0.5];
p_sig(:,:,1,3) = [1 0;0 0];
p_sig(:,:,2,1) = [0.5 0;0 0.5];
p_sig(:,:,2,2) = [0.5 0;0 0.5];
p_sig(:,:,2,3) = [0.5 0;0 0.5];
p_sig(:,:,3,1) = [0.5 0;0 0.5];
p_sig(:,:,3,2) = [0.5 0;0 0.5];
p_sig(:,:,3,3) = [0.5 0;0 0.5];

checkNorm=checkNormalization(p_sig);
checkNS=checkNonSignaling(p_sig);

fprintf('Normalization: %i\n',checkNorm);
fprintf('Non-Signalling: %i\n',checkNS);

%%

counter=0;
for i=1:size(rho_Noise,1)
   e = eigs(rho_Noise{i,1});
   if sum(e<0)
       %fprintf('%i\n',i);
       counter=counter+1;
   end
end

fprintf('Total unphysical density matrices: %i\n',counter);

