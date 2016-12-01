function boundaryLoop = getBoundaryEdges(X, T)
    nv = size(X,1);
    nf = size(T,1);

    [edges,~,trianglesToEdges] = unique(sort([T(:,1) T(:,2) ; T(:,2) T(:,3) ; T(:,3) T(:,1)],2),'rows');
    trianglesToEdges = reshape(trianglesToEdges,[],3);

    neighboringTriangles = accumarray([trianglesToEdges(:) ones(3*nf,1)],ones(3*nf,1));
    boundaryEdges = find(neighboringTriangles==1);

    be = edges(boundaryEdges,:);
    adj = sparse(be(:,1),be(:,2),ones(size(be,1),1),nv,nv);
    adj = adj+adj';
    boundaryLoop = dfsearch(graph(adj),be(1));

    % Make sure it's counter-clockwise -- magic formula
    xb = real(X(boundaryLoop));
    yb = imag(X(boundaryLoop));
    dx = circshift(xb,-1)-xb;
    yy = circshift(yb,-1)+yb;
    ind = sum(dx.*yy);

    if ind > 0 % clockwise!
        boundaryLoop = boundaryLoop(end:-1:1);
    end
end

