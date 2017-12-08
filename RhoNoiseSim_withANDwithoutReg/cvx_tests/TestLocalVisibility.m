tol = 1e-7;
import qvx.di.*;
import qvx.di.gallery.*;

S = CHSH.scenario;
P = CHSH.PRBox('P');
I = CHSH.expression('P');
P0 = Correlations.uniformlyRandom(S, 'P');

% local visibility of PR box is 1/2
cvx_begin
variable v
v*P + (1-v)*P0 == qvx.di.cvx.LocalSet(S)
maximize v
cvx_end
assert(abs(v-1/2) < tol);
