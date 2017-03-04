function indices = getIndex(pts, vertices)
    ptsX = real(pts);
    ptsY = imag(pts);
    indices = knnsearch(vertices, [ptsX, ptsY]);
end