%% Selects the relevant behaviours for qst
% Input: P-Matrix from projection to polytop
% Output: shortened P vector for QST
function P=P_QST_Selection(Behaviour)
    s=size(Behaviour);
    %P=zeros(16,1);
    
    P=[];
    
    Behaviour(2,:,2,:)=nan;
    Behaviour(:,2,:,2)=nan;
    Behaviour(2,:,3,:)=nan;
    Behaviour(:,2,:,3)=nan;
    
    
    for ii=1:3
        for jj=1:3
            P=[P,Behaviour(1,1,ii,jj),Behaviour(1,2,ii,jj)];
        end
        for jj=1:3
            P=[P,Behaviour(2,1,ii,jj),Behaviour(2,2,ii,jj)];
        end
    end
    
    
    P = reshape(P,2*2*3*3,1);
    P(isnan(P(:)))=[];
end
