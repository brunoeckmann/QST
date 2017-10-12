function CG = QuantumSupersetCG(scenario, hierarchyLevel)
% QuantumSetCG Represents an approximation quantum set using the NPA hierarchy
%
% scenario       - the scenario to consider
% hierarchyLevel - an instance of the HierarchyLevel class 
    MM = qvx.di.MomentMatrix.quantumCollinsGisin(scenario, hierarchyLevel);
    d = qvx.di.Correlations.dimension(scenario, 'NG');
    cvx_begin set sdp
    variable CG(d)
    CG == qvx.di.cvx.MomentCone(MM);
    CG(1) == 1;
    cvx_end
end
