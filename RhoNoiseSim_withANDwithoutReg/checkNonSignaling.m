function checkNS=checkNonSignaling(Behaviour)
% CHECKNONSIGNALING Checks if a Behaviour is part of the Non-Signaling
% regime given by eq (7), Brunner et al.
%
%   CHECKNS=CHECKNONSIGNALING(BEHAVIOUR) returns wheather a given Behaviour P(ab|xy) is part of NS. Returns TRUE if Behaviour does not show signaling.
%   The input Behaviour is a matrix of shape Mat(A,B,X,Y) where a=1,...,A,
%   b=1,...,B, x=1,...,X and y=1,...,Y

    checkNS=true; % true if Behaviour fullfilles NS conditions
    [Atot,Btot,Xtot,Ytot]=size(Behaviour);
    tol = eps(1e+10);
    
    for a=1:Atot
        for x=1:Xtot
            for y=1:Ytot
                Sum1=sum(Behaviour(a,:,x,y));
                for yy = 1:Ytot
                    Sum2=sum(Behaviour(a,:,x,yy));
                    % compare if the two sums are equal
                    if abs(Sum1-Sum2) > tol
                        checkNS=false;
                        %fprintf('NS Condition not fullfilled for Sum(P(%i :|%i %i)) ~= Sum(P(%i :|%i %i))',a-1,x-1,y-1,a-1,x-1,yy-1);
                        %fprintf(' (%f ~= %f)\n',Sum1,Sum2);
                        return
                    end
                end
            end
        end
    end
    
    for b=1:Btot
        for y=1:Ytot
            for x=1:Xtot
                Sum1=sum(Behaviour(:,b,x,y));
                for xx = 1:Xtot
                    Sum2=sum(Behaviour(:,b,xx,y));
                    % compare if the two sums are equal
                    if abs(Sum1-Sum2) > tol
                        checkNS=false;
                        %fprintf('NS Condition not fullfilled for Sum(P(: %i|%i %i)) ~= Sum(P(: %i|%i %i))',b-1,x-1,y-1,b-1,xx-1,y-1)
                        %fprintf(' (%f ~= %f)\n',Sum1,Sum2);
                        return
                    end
                end
            end
        end
    end
end
