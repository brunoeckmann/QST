function M = matrowgcd(M)
% simplify the rows that are not relatively prime
    for r = 1:size(M, 1)
        g = M(r, 1);
        for c = 2:size(M, 2)
            g = gcd(g, M(r,c));
        end
        if g ~= 1
            M(r, :) = M(r, :)/g;
        end
    end            
end
