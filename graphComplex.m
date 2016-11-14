function w = graphComplex(complexValues)

display(complexValues);
n = length(complexValues);

real_z = real(complexValues);
imag_z = imag(complexValues);

plot(real_z, imag_z, 'g*');
hold on;
quiver(real_z(1:n-1), imag_z(1:n-1), real_z(2:n) - real_z(1:n-1), imag_z(2:n) - imag_z(1:n-1), 0)
hold off;

end