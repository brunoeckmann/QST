%% Calculate Behaviour P(ab|xy)
% P = |<Chi_n * Chi_m | Psi>|^2 = <Chi_n * Chi_m | rho | Chi_n * Chi_m >
function P=behaviourString(MM,rho)

MString=[1,2;3,5;4,6];

%P=cell(36,1);
P=cell(2,2,3,3);
    for x = 1:3
       for y = 1:3
          for a = 1:2 
            for b = 1:2
                P{a,b,x,y}=strcat(num2str(MString(x,a)),'*',num2str(MString(y,b))); % real(kron(MM{x,a},MM{y,b})'*rho*kron(MM{x,a},MM{y,b}));
            end
          end
       end
    end
    %P=cell2mat(P);
end