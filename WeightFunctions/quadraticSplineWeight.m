
%{
This is assuming that the number of keyframes (n) is 3.

%}
function w = quadraticSplineWeight(numKeyframes, t)

if numKeyframes ~= 3
    error('must have 3 keyframes when using quadraticSplineWeight interpolation');
end

w = zeros(3, 1);
w(1) = (1-t)*(1-t);
w(2) = 2*(t*(1-t));
w(3) = t*t;

end