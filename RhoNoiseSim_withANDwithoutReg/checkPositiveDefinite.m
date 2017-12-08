function pd=checkPositiveSemiDefinite(A)
    pd=false;
    eigValues = eig((A+A')/2);
    tol=eps(1e6);
    
    if( abs((eigValues>0)-ones(size(eigValues,1),1))<= tol)
        pd=true;
    end
        
end