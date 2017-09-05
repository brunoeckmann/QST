function rho_reconstr = qst_maximumlikelihood(rhoLinearInversion,measurementBasis,countSignal)
    

    % Estimate parameters for minimization from QST_LinearInversion
    T=zeros(4,4);
    D=det(rhoLinearInversion);
    M1_11=det(rhoLinearInversion(2:4,2:4));
    M1_12=det(rhoLinearInversion(2:4,[1 3:4]));
    M2_1122=det(rhoLinearInversion(3:4,3:4));
    M2_1223=det(rhoLinearInversion(3:4,[1 4]));
    M2_1123=det(rhoLinearInversion(3:4,[2 4]));
    
    T(1,1) = sqrt(D/M1_11);
    T(2,1) = M1_12 / sqrt(M1_11*M2_1122);
    T(2,2) = sqrt(M1_11/M2_1122);
    T(3,1) = M2_1223 / sqrt(rhoLinearInversion(4,4) * M2_1122);
    T(3,2) = M2_1123 / sqrt(rhoLinearInversion(4,4) * M2_1122);
    T(3,3) = sqrt(M2_1122/rhoLinearInversion(4,4));
    T(4,1) = rhoLinearInversion(4,1)/sqrt(rhoLinearInversion(4,4));
    T(4,2) = rhoLinearInversion(4,2)/sqrt(rhoLinearInversion(4,4));
    T(4,3) = rhoLinearInversion(4,3)/sqrt(rhoLinearInversion(4,4));
    T(4,4) = sqrt(rhoLinearInversion(4,4));
    
    T(isnan(T)) = 0; 
    
    B=T'*T;
    trace(B);
    g=B/trace(B);
    
    % get start values
    start(1) = T(1,1);
    start(2) = T(2,2);
    start(3) = T(3,3);
    start(4) = T(4,4);
    start(5) = real(T(2,1));
    start(6) = imag(T(2,1));
    start(7) = real(T(3,2));
    start(8) = imag(T(3,2));
    start(9) = real(T(4,3));
    start(10) = imag(T(4,3));
    start(11) = real(T(3,1));
    start(12) = imag(T(3,1));
    start(13) = real(T(4,2));
    start(14) = imag(T(4,2));
    start(15) = real(T(4,1));
    start(16) = imag(T(4,1));
    
%     mB = cat(3,measurementBasis{:});
%     cS = reshape(countSignal,1,1,16);
%     
%     
    TT = @(t) [t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)];
    rho = @(t) TT(t)'*TT(t) / trace(TT(t)'*TT(t));
%     
%     custpdf = @(data,t) (measurementBasis{1}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{1}-data(1))^2 / 2*(measurementBasis{1}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{1})...
%          + (measurementBasis{2}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{2}-data(2))^2 / 2*(measurementBasis{2}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{2})...
%          + (measurementBasis{3}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{3}-data(3))^2 / 2*(measurementBasis{3}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{3})...
%          + (measurementBasis{4}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{4}-data(4))^2 / 2*(measurementBasis{4}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{4})...
%          + (measurementBasis{5}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{5}-data(5))^2 / 2*(measurementBasis{5}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{5})...
%          + (measurementBasis{6}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{6}-data(6))^2 / 2*(measurementBasis{6}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{6})...
%          + (measurementBasis{7}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{7}-data(7))^2 / 2*(measurementBasis{7}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{7})...
%          + (measurementBasis{8}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{8}-data(8))^2 / 2*(measurementBasis{8}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{8})...
%          + (measurementBasis{9}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{9}-data(9))^2 / 2*(measurementBasis{9}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{9})...
%          + (measurementBasis{10}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{10}-data(10))^2 / 2*(measurementBasis{10}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{10})...
%          + (measurementBasis{11}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{11}-data(11))^2 / 2*(measurementBasis{11}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{11})...
%          + (measurementBasis{12}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{12}-data(12))^2 / 2*(measurementBasis{12}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{12})...
%          + (measurementBasis{13}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{13}-data(13))^2 / 2*(measurementBasis{13}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{13})...
%          + (measurementBasis{14}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{14}-data(14))^2 / 2*(measurementBasis{14}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{14})...
%          + (measurementBasis{15}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{15}-data(15))^2 / 2*(measurementBasis{15}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{15})...
%          + (measurementBasis{16}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{16}-data(16))^2 / 2*(measurementBasis{16}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{16});
%      
%      L = @(t) (measurementBasis{1}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{1}-countSignal(1))^2 / 2*(measurementBasis{1}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{1})...
%          + (measurementBasis{2}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{2}-countSignal(2))^2 / 2*(measurementBasis{2}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{2})...
%          + (measurementBasis{3}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{3}-countSignal(3))^2 / 2*(measurementBasis{3}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{3})...
%          + (measurementBasis{4}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{4}-countSignal(4))^2 / 2*(measurementBasis{4}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{4})...
%          + (measurementBasis{5}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{5}-countSignal(5))^2 / 2*(measurementBasis{5}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{5})...
%          + (measurementBasis{6}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{6}-countSignal(6))^2 / 2*(measurementBasis{6}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{6})...
%          + (measurementBasis{7}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{7}-countSignal(7))^2 / 2*(measurementBasis{7}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{7})...
%          + (measurementBasis{8}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{8}-countSignal(8))^2 / 2*(measurementBasis{8}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{8})...
%          + (measurementBasis{9}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{9}-countSignal(9))^2 / 2*(measurementBasis{9}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{9})...
%          + (measurementBasis{10}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{10}-countSignal(10))^2 / 2*(measurementBasis{10}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{10})...
%          + (measurementBasis{11}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{11}-countSignal(11))^2 / 2*(measurementBasis{11}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{11})...
%          + (measurementBasis{12}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{12}-countSignal(12))^2 / 2*(measurementBasis{12}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{12})...
%          + (measurementBasis{13}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{13}-countSignal(13))^2 / 2*(measurementBasis{13}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{13})...
%          + (measurementBasis{14}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{14}-countSignal(14))^2 / 2*(measurementBasis{14}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{14})...
%          + (measurementBasis{15}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{15}-countSignal(15))^2 / 2*(measurementBasis{15}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{15})...
%          + (measurementBasis{16}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{16}-countSignal(16))^2 / 2*(measurementBasis{16}'*[t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)]*measurementBasis{16});

L = @(t) (measurementBasis{1}'*rho(t)*measurementBasis{1}-countSignal(1))^2 / 2*(measurementBasis{1}'*rho(t)*measurementBasis{1})...
         + (measurementBasis{2}'*rho(t)*measurementBasis{2}-countSignal(2))^2 / 2*(measurementBasis{2}'*rho(t)*measurementBasis{2})...
         + (measurementBasis{3}'*rho(t)*measurementBasis{3}-countSignal(3))^2 / 2*(measurementBasis{3}'*rho(t)*measurementBasis{3})...
         + (measurementBasis{4}'*rho(t)*measurementBasis{4}-countSignal(4))^2 / 2*(measurementBasis{4}'*rho(t)*measurementBasis{4})...
         + (measurementBasis{5}'*rho(t)*measurementBasis{5}-countSignal(5))^2 / 2*(measurementBasis{5}'*rho(t)*measurementBasis{5})...
         + (measurementBasis{6}'*rho(t)*measurementBasis{6}-countSignal(6))^2 / 2*(measurementBasis{6}'*rho(t)*measurementBasis{6})...
         + (measurementBasis{7}'*rho(t)*measurementBasis{7}-countSignal(7))^2 / 2*(measurementBasis{7}'*rho(t)*measurementBasis{7})...
         + (measurementBasis{8}'*rho(t)*measurementBasis{8}-countSignal(8))^2 / 2*(measurementBasis{8}'*rho(t)*measurementBasis{8})...
         + (measurementBasis{9}'*rho(t)*measurementBasis{9}-countSignal(9))^2 / 2*(measurementBasis{9}'*rho(t)*measurementBasis{9})...
         + (measurementBasis{10}'*rho(t)*measurementBasis{10}-countSignal(10))^2 / 2*(measurementBasis{10}'*rho(t)*measurementBasis{10})...
         + (measurementBasis{11}'*rho(t)*measurementBasis{11}-countSignal(11))^2 / 2*(measurementBasis{11}'*rho(t)*measurementBasis{11})...
         + (measurementBasis{12}'*rho(t)*measurementBasis{12}-countSignal(12))^2 / 2*(measurementBasis{12}'*rho(t)*measurementBasis{12})...
         + (measurementBasis{13}'*rho(t)*measurementBasis{13}-countSignal(13))^2 / 2*(measurementBasis{13}'*rho(t)*measurementBasis{13})...
         + (measurementBasis{14}'*rho(t)*measurementBasis{14}-countSignal(14))^2 / 2*(measurementBasis{14}'*rho(t)*measurementBasis{14})...
         + (measurementBasis{15}'*rho(t)*measurementBasis{15}-countSignal(15))^2 / 2*(measurementBasis{15}'*rho(t)*measurementBasis{15})...
         + (measurementBasis{16}'*rho(t)*measurementBasis{16}-countSignal(16))^2 / 2*(measurementBasis{16}'*rho(t)*measurementBasis{16});
         
     
    function [F,J] = LL(t,mB,cS)
        TTT = @(t) [t(1), 0, 0, 0; t(5) + 1i*t(6), t(2), 0, 0; t(11) + 1i*t(12), t(7) + 1i*t(8), t(3), 0; t(15) + 1i*t(16), t(13) + 1i* t(14), t(9) + 1i*t(10), t(4)];
        rho2 = @(t) TTT(t)'*TTT(t) / trace(TTT(t)'*TTT(t));
        F = (mB{1}'*rho2(t)*mB{1}-cS(1))^2 / 2*(mB{1}'*rho2(t)*mB{1})...
         + (mB{2}'*rho2(t)*mB{2}-cS(2))^2 / 2*(mB{2}'*rho2(t)*mB{2})...
         + (mB{3}'*rho2(t)*mB{3}-cS(3))^2 / 2*(mB{3}'*rho2(t)*mB{3})...
         + (mB{4}'*rho2(t)*mB{4}-cS(4))^2 / 2*(mB{4}'*rho2(t)*mB{4})...
         + (mB{5}'*rho2(t)*mB{5}-cS(5))^2 / 2*(mB{5}'*rho2(t)*mB{5})...
         + (mB{6}'*rho2(t)*mB{6}-cS(6))^2 / 2*(mB{6}'*rho2(t)*mB{6})...
         + (mB{7}'*rho2(t)*mB{7}-cS(7))^2 / 2*(mB{7}'*rho2(t)*mB{7})...
         + (mB{8}'*rho2(t)*mB{8}-cS(8))^2 / 2*(mB{8}'*rho2(t)*mB{8})...
         + (mB{9}'*rho2(t)*mB{9}-cS(9))^2 / 2*(mB{9}'*rho2(t)*mB{9})...
         + (mB{10}'*rho2(t)*mB{10}-cS(10))^2 / 2*(mB{10}'*rho2(t)*mB{10})...
         + (mB{11}'*rho2(t)*mB{11}-cS(11))^2 / 2*(mB{11}'*rho2(t)*mB{11})...
         + (mB{12}'*rho2(t)*mB{12}-cS(12))^2 / 2*(mB{12}'*rho2(t)*mB{12})...
         + (mB{13}'*rho2(t)*mB{13}-cS(13))^2 / 2*(mB{13}'*rho2(t)*mB{13})...
         + (mB{14}'*rho2(t)*mB{14}-cS(14))^2 / 2*(mB{14}'*rho2(t)*mB{14})...
         + (mB{15}'*rho2(t)*mB{15}-cS(15))^2 / 2*(mB{15}'*rho2(t)*mB{15})...
         + (mB{16}'*rho2(t)*mB{16}-cS(16))^2 / 2*(mB{16}'*rho2(t)*mB{16});
        
        if nargout > 1   % Two output arguments
            
        end
        
    end

    %phat = mle(countSignal,'pdf',custpdf,'start',start)
    %options = optimset('MaxIter',4000);
    options = optimset('MaxFunEvals',80000,'MaxIter',80000);

    opt = optimoptions('lsqnonlin','MaxFunctionEvaluations', 80000, 'MaxIterations',80000,'UseParallel',1, 'Algorithm', 'levenberg-marquardt');

    %x= fminsearch(L,ones(16,1),options);
    %x=fminsearch(L,start,options);
    %x = lsqnonlin(L,start,[],[],opt);
    
    x = ga(L,16);
    
    T_QST(1,1) = x(1);
    T_QST(2,2) = x(2);
    T_QST(3,3) = x(3);
    T_QST(4,4) = x(4);
    T_QST(2,1) = x(5) + 1i*x(6);
    T_QST(3,1) = x(11) + 1i*x(12);
    T_QST(3,2) = x(7) + 1i*x(8);
    T_QST(4,1) = x(15) + 1i*x(16);
    T_QST(4,2) = x(13) + 1i*x(14);
    T_QST(4,3) = x(9) + 1i*x(10);
    
    rho_reconstr = T_QST'*T_QST / trace(T_QST'*T_QST);
    
end