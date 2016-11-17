function [C] = cauchyCoordinates(cage, z)
% Compute regular Cauchy coordiantes for given cage at strictly inside points z (z are not allowed to be on the boundary)
%
% Input parameters:
% cage - cage
% z - points inside the cage

if size(z,2) == 1; z = z.'; end


Aj = cage - cage([end 1:end-1], :);
Ajp = Aj([2:end 1], :);

Bj = bsxfun(@minus, cage.', z.');
Bjp = Bj(:, [2:end 1]);
Bjm = Bj(:, [end 1:end-1]);


oneOver2pi_i = 1/(2*pi*1i);

%C = oneOver2pi_i*((Bjp./Ajp).*log(Bjp./Bj) - (Bjm./Aj).*log(Bj./Bjm));

C = oneOver2pi_i*(Bjp.*bsxfun(@rdivide, log(Bjp./Bj), Ajp.') - Bjm.*bsxfun(@rdivide, log(Bj./Bjm), Aj.'));

end