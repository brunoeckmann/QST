%% Calculate Behaviour P(ab|xy)
% P = |<Chi_n * Chi_m | Psi>|^2 = <Chi_n * Chi_m | rho | Chi_n * Chi_m >
function P=behaviour(MM,rho)

%P=cell(36,1);
P=zeros(2,2,3,3);
    for x = 1:3
       for y = 1:3
          for a = 1:2 
            for b = 1:2
                %P{12*(x-1)+4*(y-1)+2*(a-1)+b}=real(kron(MM{x,a},MM{y,b})'*rho*kron(MM{x,a},MM{y,b}));
                P(a,b,x,y)=real(kron(MM{x,a},MM{y,b})'*rho*kron(MM{x,a},MM{y,b}));
            end
          end
       end
    end
    %P=cell2mat(P);
end

