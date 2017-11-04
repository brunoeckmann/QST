function D=kullback_leibler_divergence(behaviourReg, behaviourNoise)
% Calculate the Kullback-Leibler Divergence between two Behaviours P
% Behaviour must be P(a b | x y)
% a Matrix A x B x X x Y 
% a=1..A, b=1..B, x=1...X, y=1...Y
[A B X Y]=size(behaviourReg);
D=0;
for x=1:X
   for y=1:Y 
       for a=1:A
            for b=1:B
                if behaviourNoise(a,b,x,y)~=0
                    D=D+behaviourReg(a,b,x,y)*log(behaviourReg(a,b,x,y) / behaviourNoise(a,b,x,y));
                end
            end
       end
   end
end


end