function w = graphComplex(complexValues, loop)

if nargin < 2
   loop = true;
end

display(complexValues);
n = length(complexValues);

real_z = real(complexValues);
imag_z = imag(complexValues);

plot(real_z, imag_z, 'g*');
hold on;
quiver(real_z(1:n-1), imag_z(1:n-1), real_z(2:n) - real_z(1:n-1), imag_z(2:n) - imag_z(1:n-1), 0)
if loop
    quiver(real_z(n), imag_z(n), real_z(1) - real_z(n), imag_z(1) - imag_z(n))
end
hold off;

end