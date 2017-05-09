function argumentPrincipleIntegral = integralOfDerivativeOfLog_fz_onPath(cage, phi, pathVertices, prodFunc)
%insidePolygon must be fully contained in the cage (which is also a polygon).
%for example insidePolygon can be an inward offset of the cage polygon.


argumentPrincipleIntegral = integral(@(z)...
    prodFunc(z) * derivativeOfLog_fz(cage, z, phi), ...
    complex(pathVertices(1,:)), ...
    complex(pathVertices(end,:)), ...
    'Waypoints', complex(pathVertices(2:end-1,1), pathVertices(2:end-1,2))...
);

end


function r = derivativeOfLog_fz(cage, z, phi)
%z can be a scalar or vector of complex points

[D, E] = derivativesOfCauchyCoord(cage, z);
%size(z)
r = (E*phi)./(D*phi);
r = gather(r.');

end