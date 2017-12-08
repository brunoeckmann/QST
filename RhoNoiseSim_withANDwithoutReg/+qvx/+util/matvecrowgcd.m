function [A b] = matvecrowgcd(A, b)
% simplify the rows of the linear system (A*x, b) that are not relatively prime
    for r = 1:size(A, 1)
        g = b(r);
        for c = 1:size(A, 2)
            g = gcd(g, A(r,c));
        end
        if g ~= 1
            A(r, :) = A(r, :)/g;
            b(r) = b(r)/g;
        end
    end            
end
