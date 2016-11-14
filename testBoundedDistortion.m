%% Problem setup

clear

% We'll randomly generate a deformation with these parameters, just to
% demonstrate what the algorithm is doing
nhandles = 6;
displacementAmount = 0.05;

% Deformation bounds, seem to be reasonable from the paper (Fig. 4)
k = .2; % seem to be reasonable from Figure 4 in the paper
sigma2 = .5;
sigma1 = 1.5;
lambda = 100;

% Load meshj
[X,T] = readOff('Meshes/llr2.off');
X = X / max(abs(X(:))); % rescale to unit box
X = X(:,1)+1i*X(:,2); % convert vertices to complex numbers

% Number of vertices/triangles
nv = size(X,1);
nf = size(T,1);

% Randomly choose handles and their displacements to demo the deformation
handles = randperm(nv,nhandles);
r = X(handles); % rest positions
q = r; % fix all handles
q(1) = q(1) + displacementAmount*(randn + 1i*randn); % but perturb one randomly

% Display input mesh with handles
trimesh(T,real(X),imag(X),'color','k');
axis equal; axis off; hold on;
for i=1:nhandles, plot([r(i) q(i)],'b-'); end % plot destination of each handle
plot(X(handles),'r.','markersize',50,'linewidth',10);
title('Desired deformation');

%% Find boundary of the triangle mesh

[edges,~,trianglesToEdges] = unique(sort([T(:,1) T(:,2) ; T(:,2) T(:,3) ; T(:,3) T(:,1)],2),'rows');
trianglesToEdges = reshape(trianglesToEdges,[],3);

neighboringTriangles = accumarray([trianglesToEdges(:) ones(3*nf,1)],ones(3*nf,1));
boundaryEdges = find(neighboringTriangles==1);

be = edges(boundaryEdges,:);
adj = sparse(be(:,1),be(:,2),ones(size(be,1),1),nv,nv);
adj = adj+adj';
boundaryLoop = dfsearch(graph(adj),be(1));

nb = length(boundaryLoop);

% Make sure it's counter-clockwise -- magic formula
xb = real(X(boundaryLoop));
yb = imag(X(boundaryLoop));
dx = circshift(xb,-1)-xb;
yy = circshift(yb,-1)+yb;
ind = sum(dx.*yy);

if ind > 0 % clockwise!
    boundaryLoop = boundaryLoop(end:-1:1);
end
                         
%% Do Cauchy coordinates

boundaryPositions = X(boundaryLoop); % nb x 1

% Compute vertex normals on boundary
edgeTangent = circshift(boundaryPositions,-1)-boundaryPositions;
edgeNormal = imag(edgeTangent)-1i*real(edgeTangent);
vtxNormal = .5*(edgeNormal + circshift(edgeNormal,1));
vtxNormal = vtxNormal ./ abs(vtxNormal);

% Dispalce boundary
boundaryPositions = boundaryPositions + vtxNormal*.001*sum(abs(edgeTangent));

% Use notation of the paper to define assorted useful variables
% See figure 2 and appendix A
B = repmat(boundaryPositions',nv,1) - repmat(X,1,nb); % nv x nb
A = diff(boundaryPositions)'; 
A(end+1) = boundaryPositions(1)-boundaryPositions(end); % 1 x nb
A = circshift(A,1); % A is difference z_j - z_{j-1}!

nextA = circshift(A,-1); % 1 x nb

nextB = circshift(B,-1,2); % nv x nb
prevB = circshift(B,1,2); % nv x nb

C = ( bsxfun(@rdivide,nextB,nextA).*log(nextB./B) ...
     -bsxfun(@rdivide,prevB,A).*log(B./prevB) )/(2*pi*1i);

Cprime = ( bsxfun(@rdivide,log(B./nextB),nextA) ...
          +bsxfun(@rdivide,log(B./prevB),A) )/(2*pi*1i);

%% Compute map

theta = zeros(nb,1);

M = boundaryLoop; % not bothering with their active set method
nm = length(M);

clear Phi Psi 

for i=1:3 % paper says 1-3 iterations are enough...
    cvx_begin % see equation (19)
        cvx_solver mosek
        cvx_precision low % foresake quality for speed

        variable Phi(nb,1) complex
        variable Psi(nb,1) complex

        variable f(nv,1) complex
        variable fz(nv,1) complex
        variable fzbar(nv,1) complex

        minimize sum_square_abs(fz(M).*exp(1i*theta) - 1)/nm...
                   +sum_square_abs(fzbar(M))/nm...
                   +lambda*sum_square_abs(f(handles) - q)
        subject to
            % Set up the basis
            f == C*Phi + conj(C*Psi); % slow -- only really need on boundary (should fix this)
            fz == Cprime*Phi;
            fzbar == conj(Cprime*Psi);

            % Constraints from the paper
            Psi(1) == 0;
            abs(fzbar(boundaryLoop)) <= k*real(fz(boundaryLoop).*exp(1i*theta))
            abs(fz(boundaryLoop))+abs(fzbar(boundaryLoop)) <= sigma1
            abs(fzbar(boundaryLoop)) <= real(fz(boundaryLoop).*exp(1i*theta)) - sigma2
    cvx_end

    theta = -angle(fz(boundaryLoop)); % Update angles for next iteration
end

% Display deformed mesh
figure
trimesh(T,real(f),imag(f),'color','k');
axis equal; axis off; hold on;
title(sprintf('Iteration %d',i));
plot(X(handles),'g.','markersize',50,'linewidth',10);
plot(f(handles),'r.','markersize',50,'linewidth',10);
plot(q,'b.','markersize',50,'linewidth',10);
drawnow;