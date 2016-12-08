function killme = interpolationPreparation(X,T)


    fprintf('\n Running bdhInterpPrep');
    assert( all(signedAreas(X, T)>0), 'source triangulation not in correct order!');

    %%
    nBEMSample = 20;
    BEMSampleOffset = cage_offset; 

    nVirtual = 1;
    mxsub = subdivPolyMat(offsetCage, nVirtual*size(offsetCage,1));

    nx = size(mxsub, 1);
    S = subdivPolyMat(mxsub*offsetCage, ceil(nBEMSample*nx));
    w = S*mxsub*polygonOffset(offsetCage, BEMSampleOffset, false);

    %% init, precomputations
    % [C0, D0] = regularCauchyCoord(mxsub*offsetCage, w);
    C0 = cauchyCoordinates(mxsub*offsetCage, w.');
    D0 = derivativesOfCauchyCoord(mxsub*offsetCage, w.');

    C = cauchyCoordinates(mxsub*offsetCage, X);
    D = derivativesOfCauchyCoord(mxsub*offsetCage, X);

    C2 = cauchyCoordinates(offsetCage, X);

    % assert( norm(sum(C0,2)-1)<1e-8, 'samples outside cage?' );
    % precomputed data
    invC0 = pinv( [real(C0) -imag(C0)] );
    invC0 = complex(invC0(1:nx,:), invC0(nx+1:2*nx,:));
    invD0 = pinv(D0);


    if ~exist('PhiPsyKF', 'var')
        PhiPsyKF = [offsetCage*[1 0] Phi Psy];
    end


    C1 = cauchyCoordinates(offsetCage, w.');
    % invC1 = pinv( [real(C1) -imag(C1)] );
    % invC1 = complex(invC1(1:end/2,:), invC1(end/2+1:end,:));


    phipsy = invC0*real(C1*PhiPsyKF);
    phipsy = phipsy + repmat( imag(mean(PhiPsyKF)-mean(phipsy))*1i, size(phipsy,1), 1); % fix global translation in y-axis


    numKeyFrame = size(phipsy,2)/2;
    ikeyframe = min(numKeyFrame-1, ikeyframe);


    keyframes = (1:numKeyFrame)*2-1;


    if ~var(keyframes), warning('interpolating same shape!'); end

    fprintf('### interpolating %d keyframes\n', numKeyFrame);

    fz = D0*phipsy(:, keyframes);
    fzbar = D0*phipsy(:, keyframes+1);  % fzbarbar actaully

    fzX = D*phipsy(:, keyframes);
    fzbarX = D*phipsy(:, keyframes+1);

    fWtFun = @(w) linearWeight(numKeyFrame, w);

    fComposeHarmonicMap = @(phipsy) phipsy(:,1:2:end)+conj(phipsy(:,2:2:end));



    %% extract g
    gFromFz = @(fz) complex( log(abs(fz)), angle(fz(end)) + cumsum(angle(fz./circshift(fz, 1))) );
    g = gFromFz(fz);


    TR = triangulation(T, real(X), imag(X));
    edges = TR.edges;
    nv = numel(X);

    gdifFromFz = @(fz) complex(log(abs(fz(edges(:,2),:)./fz(edges(:,1),:))), angle( fz(edges(:,2),:)./fz(edges(:,1),:) ));
    MV2Evec = sparse( repmat(1:size(edges,1), 2, 1)', edges, repmat( [-1 1], size(edges,1), 1), size(edges,1), nv );
    gX = [MV2Evec; sparse(1, 1, 1, 1, numel(X))]\[gdifFromFz(fzX); log(fzX(1,:))];


    %% for numerical integration
    if exist('interpAnchID', 'var') && ~isempty(interpAnchID)
        anchorVertexIndex = interpAnchID(1);   
    else
        anchorVertexIndex = 1; 
    end

    adjacencyGraph = sparse(edges(:, 2), edges(:, 1), 1, nv, nv);
    assert(nnz(triu(adjacencyGraph)) == 0);
    [disc, pred, closed] = graphtraverse(adjacencyGraph, anchorVertexIndex, 'Directed', false, 'Method', 'BFS');
    endIndices = uint32(disc(2:end))';
    startIndices = uint32(pred(endIndices))';
    eVec = X(endIndices) - X(startIndices);
    e4treecumsum = [startIndices endIndices];





    if exist('bdhiMethod', 'var')~=1, bdhiMethod = 'metric'; end

    fprintf('### prepare to do BDH interpolate by %s\n', bdhiMethod);

    if strcmp( bdhiMethod, 'nu' )
    %% interp nu, does not maintain/interp stretch direction
        fInterpFz = @(wt) exp(g*fWtFun(wt)); % here g is just the definition of log(f_z)
        nu = fzbar./fz;
        fInterpNu = @(wt) nu*fWtFun(wt);
        fInterpFzbar = @(wt) fInterpFz(wt).*fInterpNu(wt);

        nuX = fzbarX./fzX;
        fInterpNuX = @(wt) nuX*fWtFun(wt);
        fInterpFzX = @(wt) exp(gX*fWtFun(wt));
        fInterpFzbarX = @(wt) fInterpFzX(wt).*fInterpNuX(wt);
        fInterpEtaX = @(wt) fInterpFzX(wt).^2.*fInterpNuX(wt);
    elseif strcmp( bdhiMethod, 'eta' )
        scaleEta = 1;
        scaleGlobal = 1;
        boundGlobalk = 0;

        if boundGlobalk
            kmax = max(abs(fzbar(:)./fz(:)));
        else
            kmax = max(abs(fzbar./fz), [], 2);
        end

        %% interp with BEM
        fInterpEta0  = @(wt) (fz.*fzbar)*fWtFun(wt);
        fInterpFz0   = @(wt) exp(g*fWtFun(wt));
        fInterpk0    = @(wt) abs(fInterpEta0(wt)./fInterpFz0(wt).^2);
        fInterpFzbar0 = @(wt) fInterpEta0(wt)./fInterpFz0(wt);

        fMyHilbert = @(s) C0*invC0*s - 1i*mean(imag(C0*invC0*s));
        if scaleGlobal
            maxsiga = max( abs(fz)+abs(fzbar), [], 2 );
            minsigb = min( abs(fz)-abs(fzbar), [], 2 );

            fRhoa = @(wt) (maxsiga-abs(fInterpFz0(wt)))./abs(fInterpFzbar0(wt));
            fRhob = @(wt) (abs(fInterpFz0(wt))-minsigb)./abs(fInterpFzbar0(wt));

            fSScale = @(wt) min(1, min( min([kmax./fInterpk0(wt), fRhoa(wt), fRhob(wt)], [], 2) ) );
            fSScale = @(wt) 1; % for bounded distortion 
        else
            fSScale = @(wt) exp( fMyHilbert( min(0, log(kmax./fInterpk0(wt))) ) );
        end

        fInterpEta   = fInte,kkrpEta0;
        fInterpFz    = fInterpFz0;

        fInterpEtaX0  = @(wt) fzX.*fzbarX*fWtFun(wt);
        fInterpFzX0   = @(wt) exp(gX*fWtFun(wt));
        fInterpEtaX  = fInterpEtaX0;
        fInterpFzX   = fInterpFzX0;

        if ~scaleEta
            fInterpFz   = @(wt) fSScale(wt).^-0.5.*fInterpFz0(wt);
            fInterpFzX  = @(wt) fSScale(wt).^-0.5.*fInterpFzX0(wt);
        else
            fInterpEta  = @(wt) fSScale(wt).*fInterpEta0(wt);
            fInterpEtaX = @(wt) fSScale(wt).*fInterpEtaX0(wt);
        end

        fInterpFzbar = @(wt) fInterpEta(wt)./fInterpFz(wt);
    elseif strcmp( bdhiMethod, 'metric' )
        % can be linear, hermite, quadratic
        etaInterpType = 'hermite';

        %% abs2(fz) interp
        fMyHilbert = @(s) C0*(invC0*s);


        eta = fz.*fzbar;
        dfnorm2 = abs(fz).^2 + abs(fzbar).^2;

        if strcmp(etaInterpType, 'hermite')
            fA = @(wt) (hermiteWeight(numKeyFrame, wt, abs(eta))).^2;
            fB = @(wt) hermiteWeight(numKeyFrame, wt, dfnorm2);
            fAB2LogFz2 = @(A, B) log( (B+sqrt(B.^2-4*A))/2 );
            fInterpFz = @(wt) exp( 0.5*fMyHilbert( fAB2LogFz2( fA(wt), fB(wt) ) ) );
            fInterpEta  = @(wt) hermiteWeight(numKeyFrame, wt, eta, ones(size(eta)) );
            fInterpFzbar = @(wt) fInterpEta(wt)./fInterpFz(wt);


            fMyHilbertX = @(s) C*(invC0*s);
            fInterpFzX = @(wt) exp( 0.5*fMyHilbertX( fAB2LogFz2( fA(wt), fB(wt) ) ) );
            etaX = fzX.*fzbarX;
            fInterpEtaX = @(wt) hermiteWeight(numKeyFrame, wt, etaX, ones(size(etaX)) );
        elseif strcmp(etaInterpType, 'linear')

            fA = @(wt) (abs(eta)*fWtFun(wt)).^2;
            fB = @(wt) dfnorm2*fWtFun(wt);
            fAB2LogFz2 = @(A, B) log( (B+sqrt(B.^2-4*A))/2 );
            fInterpFz = @(wt) exp( 0.5*fMyHilbert( fAB2LogFz2( fA(wt), fB(wt) ) ) );

            % conformal only
            % fInterpFz = @(wt) exp( 0.5*fMyHilbert( log( abs(fz(:,1:2)).^2*fWtFun(wt) ) ) );  % interp fz(t)^2
            % fInterpFzbar = @(wt) fzbar(:,1)*0;

            fInterpEta  = @(wt) eta*fWtFun(wt);
            fInterpFzbar = @(wt) fInterpEta(wt)./fInterpFz(wt);

            fMyHilbertX = @(s) C*(invC0*s);
            fInterpFzX = @(wt) exp( 0.5*fMyHilbertX( fAB2LogFz2( fA(wt), fB(wt) ) ) );
            etaX = fzX.*fzbarX;
            fInterpEtaX = @(wt) etaX*fWtFun(wt);
        end

    end


    % global pose
    if ~exist('interpAnchID', 'var'), interpAnchID = []; end
    if isempty(interpAnchID), warning('anchors for interpolation not set!'); end

    anchSrc = fComposeHarmonicMap( C(interpAnchID, :)*phipsy(:, reshape([keyframes; keyframes+1], [], 1)) );



    fBdhInterpX2 = @(wt) XFromFzEta(fInterpFzX(wt), fInterpEtaX(wt), eVec, e4treecumsum, uint32(anchorVertexIndex));
    fBdhInterpX = @(wt) alignPoseByAnchors( fBdhInterpX2(wt), interpAnchID, anchSrc, fWtFun(wt) );

end



