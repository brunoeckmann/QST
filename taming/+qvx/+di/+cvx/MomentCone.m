function x = MomentCone(momentMatrix)
    f = momentMatrix.nFreeMoments;
    d = size(momentMatrix.indexMatrix, 1);
    n = momentMatrix.outputDim;
    cvx_begin set sdp
    variable x(n)
    variable freeMoments(f)
    chi = zeros(d, d);
    for i = 1:n
        chi = chi + momentMatrix.matrixForMoment(i) * x(i);
    end
    for i = 1:f
        chi = chi + momentMatrix.matrixForMoment(i + n) * freeMoments(i);
    end
    if ~isequal(momentMatrix.iA, [])
        momentMatrix.iA * [x;freeMoments] <= momentMatrix.ib
    end
    chi == semidefinite(d)
    cvx_end
end
