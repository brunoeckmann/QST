function PProj = FindProjection(Pobs)
    S = qvx.di.Scenario.fromPabxy(Pobs);
    PProj = qvx.di.Correlations.nonsignalingProjection(S, 'Pabxy', Pobs);
end
