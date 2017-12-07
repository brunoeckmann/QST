function pd=checkPositiveSemiDefinite(A)
    pd=false;
    eigValues = eig((A+A')/2);
    tol=eps(1e6);
    
    if(eigValues>=0)
        pd=true;
    end
        
end