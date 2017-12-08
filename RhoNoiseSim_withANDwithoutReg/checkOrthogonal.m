%% Check if all vectors inside extended measurment basis are orthogonal
% MeasurementBasis is a cell(i,j)
function r=checkOrthogonal(MeasurementBasis)
    [I,J]=size(MeasurementBasis);
    tol = eps(1);

    r=true;
    for ii=1:I
       for j=1:J
          for jj=1:J
              dotProduct = dot(MeasurementBasis{ii,j},MeasurementBasis{ii,jj});
              if dotProduct > tol && j~=jj
                  fprintf('MM{%f,%f}xMM{%f,%f}\n',ii,j,ii,jj);
                  r=false;
              end
          end
       end
    end
end
