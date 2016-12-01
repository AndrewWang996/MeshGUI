%% Problem setup

clear

% We'll randomly generate a deformation with these parameters, just to
% demonstrate what the algorithm is doing
nhandles = 3;
displacementAmount = 0.05;
boundaryNormalLength = 0.0001;

% Deformation bounds, seem to be reasonable from the paper (Fig. 4)
k = .2; % seem to be reasonable from Figure 4 in the paper
sigma2 = .5;
sigma1 = 1.5;
lambda = 100;

% Load meshj
filename = 'Meshes/elephant.off';
[X,T] = readOff(filename);
X = X(:,1)+1i*X(:,2); % convert vertices to complex numbers

%% choose boundary positions
boundaryPositions = boundaryPicker(filename);
nb = length(boundaryPositions);

%% Number of vertices/triangles
nv = size(X,1);
nf = size(T,1);

%% Randomly choose handles and their displacements to demo the deformation
handles = randperm(nv,nhandles);
r = X(handles); % rest positions
q = r; % fix all handles
q(1) = q(1) + displacementAmount*(randn + 1i*randn); % but perturb one randomly

%% Display input mesh with handles
trimesh(T,real(X),imag(X),'color','k');
axis equal; axis off; hold on;
for i=1:nhandles, plot([r(i) q(i)],'b-'); end % plot destination of each handle
plot(X(handles),'r.','markersize',50,'linewidth',10);
title('Desired deformation');





%% Do Cauchy coordinates

C = cauchyCoordinates(boundaryPositions, X);
Cprime = derivativesOfCauchyCoord(boundaryPositions, X);

%% Find boundary of the triangle mesh

boundaryLoop = getBoundaryEdges(X,T);

%% Compute map

M = boundaryLoop; % not bothering with their active set method
nm = length(M);

theta = zeros(nm,1);

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






