function plotTree(vertices, edges)
    x = vertices(edges(:,1), 1);
    y = vertices(edges(:,1), 2);
    u = vertices(edges(:,2), 1) - x;
    v = vertices(edges(:,2), 2) - y;
    
    quiver(x,y,u,v,1);
end