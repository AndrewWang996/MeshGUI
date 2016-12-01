function [D, E] = derivativesOfCauchyCoord(cage, z)
% Compute first and second derivatives of regular Cauchy coordiantes
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

D = oneOver2pi_i*(bsxfun(@rdivide, log(Bj./Bjm), Aj.') - bsxfun(@rdivide, log(Bjp./Bj), Ajp.'));

if(nargout > 1)
    E = oneOver2pi_i*(Bjp-Bjm)./(Bjm.*Bj.*Bjp);
end

end
