%{
This function takes an integer n and a float t between [0,1].
In practice n is commonly given as the number of keyframes while 
    t is used as time.

%}
function w = linearWeight(n, t)

    t2 = t*(n-1)+1;
    indices = floor(t2);
    w = zeros(n, length(t));
   
    for i = 1 : length(indices)
        w( indices(i), i ) = 1 + indices(i) - t2(i);
        if indices(i) < n
            w( indices(i) + 1, i ) = t2(i) - indices(i);
        end
    end
end
