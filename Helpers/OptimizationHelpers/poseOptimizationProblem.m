function [H,Aeq,beq] = poseOptimizationProblem(meshname, indices, velocities, whichKeyframe, anchorIndex)
    [vertices, faces] = getMesh(meshname);
    nV = size(vertices,1);
    L = - cotmatrix(vertices, faces);
    H = horzcat( vertcat(L, zeros(nV)), vertcat(zeros(nV), L) );

    Aeq_comp = getAeq(meshname, indices, whichKeyframe, anchorIndex);
    Aeq_real = real(Aeq_comp);
    Aeq_imag = imag(Aeq_comp);
    Aeq = horzcat( vertcat(Aeq_real, Aeq_imag), vertcat(-Aeq_imag, Aeq_real) );

    beq_comp = getBeq(meshname, indices, velocities, whichKeyframe, anchorIndex);
    beq_real = real(beq_comp);
    beq_imag = imag(beq_comp);
    beq = vertcat(beq_real, beq_imag);
end