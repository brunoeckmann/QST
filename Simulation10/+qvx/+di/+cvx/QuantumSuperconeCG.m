function CG = QuantumSuperconeCG(scenario, hierarchyLevel)
% QuantumSuperConeCG Represents a conic extension of QuantumSupersetCG
    MM = qvx.di.MomentMatrix.quantumCollinsGisin(scenario, hierarchyLevel);
    d = qvx.di.Correlations.dimension(scenario, 'NG');
    cvx_begin set sdp
    variable CG(d)
    CG == qvx.di.cvx.MomentCone(MM);
    cvx_end
end
