function coordinates = toComples(coords)
    coordinates = coords(1,:) + 1i * coords(2,:);
end

function coordinates = toNormal(coords)
    coordinates = [real(coords); imag(coords)];
end