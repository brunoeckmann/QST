function [PNA flag prob] = FindNA_KL_PENLAB(PObs, hierarchyLevel, options)
% FindNA_KL_PENLAB Find the nearest (approximate quantum)/(nonsignaling) distribution
%
% This variant parameterizes the candidate distribution in the Collins-Gisin basis,
% and computes the distance in the Pabxy space using the Kullback-Leibler divergence.
%
% INPUTS
% PObs           Observed frequencies, of size [A B.. X Y..], with indices (a,b..,x,y..)
% hierarchyLevel Hierarchy level to use (of type qvx.di.HierarchyLevel)
%                Pass [] to use the nonsignaling polytope as output set
% options        Options to pass to PENLAB, defaults to []
%
% Note: non-negativity constraints are *always* added, because they are required for the 
%       KL-divergence to be properly defined
%
% OUTPUTS
% PNA            Nearest approximation obtained
% flag           Returns 1 if optimal, 2 if suboptimal, > 2 if error (see prob.solvermsg{flag})
% prob           PENLAB problem object
    
    if nargin < 3
        options = [];
    end
    
    S = qvx.di.Scenario.fromPabxy(PObs);
    % converts in the internal vector format
    PObs = qvx.di.Correlations.convert(S, 'P', 'Pabxy', PObs);

    CG = S.conversionMatrix('P', 'NG');
    c = CG(:,1);
    a = CG(:,2:end);
    penm = [];
    penm.probname = 'KL minimization';
    penm.comment = 'SDP constraints with nonlinear objective';
    dP = length(PObs); % dim of the full P
    dCG = size(CG, 2); % dim of the CG param.
    dCG1 = dCG - 1; % removing normalization
    % variables 1..dP = the full PReg(abxy)
    % variables dP+1..dP+dCG-1 = the Collins-Gisin param. without the constant 1
    % variables dP+dCG.. = additional moments
    if ~isequal(hierarchyLevel, [])
        
        % for the basis
        % 1 = constant
        % 2..dCG+1 = Collins-Gisin moments
        % dCG+2.. = additional hidden moments
        MM = qvx.di.MomentMatrix.quantumCollinsGisin(S, hierarchyLevel);
        CGBasis = cell(1, dCG);
        for i = 1:dCG
            CGBasis{i} = MM.matrixForMoment(i);
        end
        HigherBasis = cell(1, MM.nFreeMoments);
        for i = 1:MM.nFreeMoments
            HigherBasis{i} = MM.matrixForMoment(i + dCG);
        end
        dim = size(CGBasis{1}, 1);
        penm.NALIN = 1;
        penm.lbA = zeros(1,1); % lower-bound on eigenvalues is 0
        penm.mconfun = @mymconfun;
        penm.mcongrad = @mymcongrad;
        Nh = length(HigherBasis);
    else
        Nh = 0;
    end
    penm.Nx = dP + dCG1 + Nh;
    penn.NY = 0;
    % lower bounds for the P(ab|xy)
    penm.lbx = [zeros(1,dP) ones(1,dCG1+Nh)*(-Inf)];
    penm.ubx = [ones(1,dP+dCG1+Nh)*Inf];
    penm.lbxbar = 1:dP; % the bounds on the variables P(ab|xy) are strict > 0
    penm.ubxbar = [];
    penm.objfun = @myobjfun;
    penm.objgrad = @myobjgrad;
    penm.objhess = @myobjhess;
    % equality constraints between the P(ab|xy) and the CG formulation
    penm.NgLIN = dP;
    penm.confun = @myconfun;
    penm.congrad = @mycongrad;
    penm.lbg = zeros(dP, 1);
    penm.ubg = zeros(dP, 1);
    % defines the problem
    prob = penlab(penm)
    xinit = zeros(dP + dCG1 + Nh, 1);
    % initialize with the uniformly random distribution
    Pur_P = qvx.di.Correlations.uniformlyRandom(S, 'P');
    Pur_CG = qvx.di.Correlations.uniformlyRandom(S, 'NG');
    xinit(1:dP) = Pur_P;
    xinit(dP+1:dP+dCG1) = Pur_CG(2:end);
    if Nh > 0
        xinit(dP+dCG1+1:end) = min(Pur_CG(2:end));
    end
    prob.xinit = xinit;
    prob.opts = options;
    flag = prob.solve();
    x = prob.x;
    PNA = x(1:dP);
    PNA = qvx.di.Correlations.toPabxy(S, PNA);
    
    function [g userdata] = mymconfun(x, Y, k, userdata)
    % Evaluation of SDP matrix for the NPA hierarchy
        g = CGBasis{1};
        for i = 1:dCG1
            g = g + CGBasis{i+1} * x(i+dP);
        end
        for i = 1:Nh
            g = g + HigherBasis{i} * x(i+dP+dCG1);
        end
    end
    function [dg userdata] = mymcongrad(x, Y, k, i, userdata)
    % Gradient of SDP matrix for the NPA hierarchy (is linear, so a constant)
        if i > dP+dCG1
            dg = sparse(HigherBasis{i-dP-dCG1});
        elseif i > dP
                dg = sparse(CGBasis{i-dP+1}); % (for CHSH scenario) 17..24 maps to 2..9
        else
            dg = sparse(dim, dim); % or zero
        end
    end
    
    function [g userdata] = myconfun(x, Y, userdata)
    % Equality constraints
        g = x(1:dP) - c - a * x(dP+1:dP+dCG1);
        assert(isequal(size(g), [dP 1]));
    end
    function [dg userdata] = mycongrad(x, Y, userdata)
        dg = sparse(dP+dCG1+Nh, dP);
        dg(1:dP,1:dP) = eye(dP);
        dg(dP+1:dP+dCG1,1:dP) = -a';
        assert(isequal(size(dg), [dP+dCG1+Nh dP]));
    end
    
    function [kl userdata] = myobjfun(x, Y, userdata)
    % Computation of KL divergence
        kl = 0;
        for i = 1:dP
            f = PObs(i);
            p = x(i);
            if f > 0
                kl = kl + f * log( f/p );
            end
        end
    end
    function [dkl userdata] = myobjgrad(x, Y, userdata)
        dkl = sparse(dP+dCG1+Nh, 1);
        for i = 1:dP
            f = PObs(i);
            p = x(i);
            dkl(i) = -f/p;
        end
    end
    function [dkl2 userdata] = myobjhess(x, Y, userdata)
        dkl2 = sparse(dP+dCG1+Nh, dP+dCG1+Nh);
        for i = 1:dP
            f = PObs(i);
            p = x(i);
            dkl2(i,i) = f/p/p;
        end
    end
                       
end
