function G = meshToGraph(vertices, faces)
    edges = zeros(size(faces,1) * 3, 2);
    for i = 1:size(faces, 1)
        a = faces(i,1);
        b = faces(i,2);
        c = faces(i,3);
        
        edges(3*i - 2, 1) = a;
        edges(3*i - 2, 2) = b;
        
        edges(3*i - 1, 1) = a;
        edges(3*i - 1, 2) = c;
        
        edges(3*i - 0, 1) = b;
        edges(3*i - 0, 2) = c;
    end
    
    edges = unique(sort(edges,2), 'rows');
    edgeA = round(edges(:,1));
    edgeB = round(edges(:,2));
    
    edgeLengths = getDistance( vertices(edgeA), vertices(edgeB) );
    G = graph(edgeA, edgeB, edgeLengths);
end

function dist = getDistance(p1, p2)
    dist = abs(p1 - p2);
end

