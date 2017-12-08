function check=checkNormalization(Behaviour)
% CHECKNORMALIZATION Checks if a Behaviour P(a,b|x,y) is normalized
% Sum(:,:,x,y) != 1 for all x,y
%  
%   Behaviour is a (A x B x Y x Y) - Matrix
%   a = 1...A
%   b = 1...B
%   x = 1...X
%   y = 1...Y
%
%   Returns true if Behaviour is normalized


    check=true;
    tol=eps(1e+9);

    [A,B,X,Y]=size(Behaviour);
    for x=1:X
        for y=1:Y
            if abs(sum(reshape(Behaviour(:,:,x,y),A*B,1))-1) > tol
                check = false;
                return;
            end
        end
    end


end
