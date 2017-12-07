%% Measurement
% Signal: create count signal
function n_count = countSignal(measurementBasis, densityMatrix)
    s = size(measurementBasis);
    n_count=zeros(s(1),1);
    for nu=1:s(1)
       n_count(nu)=measurementBasis{nu}'*densityMatrix*measurementBasis{nu};
    end
end

