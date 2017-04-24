
function w = hermiteWeight(n, t)

t2 = t*(n-1)+1;
i = floor(t2);
w = zeros(2*n, 1);

t3 = t2 - i;

w(i) = 2*t3^3 - 3*t3^2 + 1;
w(n+i) = t3^3 - 2*t3^2 + t3;
if i<n
    w(i+1) = -2*t3^3 + 3*t3^2;
    w(n+i+1) = t3^3 - t3^2;
end
    

end