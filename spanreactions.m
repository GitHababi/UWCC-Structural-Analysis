function reactionForces = spanreactions(X,load,supportPositions,Idistribution,E)
%SPANREACTIONS Solves for reactions at positions specified using beam FEM
% X: row vector of node positions
% load: row vector of distributed load q(x) (negative = downward)
% supportPositions: row vector of positions where v = 0
% Idistribution: row vector of I(x)
% E: modulus

n = length(X);
ndof = 2*n; % [v1 θ1 v2 θ2 ...]

K = zeros(ndof);
F = zeros(ndof,1);

%%--- Assemble stiffness and load ---
for i = 1:n-1
    Le = X(i+1) - X(i);
    
    % average EI over element
    EI = E * (Idistribution(i) + Idistribution(i+1))/2;
    
    % element stiffness matrix
    ke = EI/Le^3 * ...
        [12     6*Le   -12     6*Le;
         6*Le   4*Le^2 -6*Le   2*Le^2;
        -12    -6*Le    12    -6*Le;
         6*Le   2*Le^2 -6*Le   4*Le^2];
    
    % consistent load vector (uniform q over element)
    q = (load(i) + load(i+1))/2; % average load
    fe = q*Le/2 * [1; Le/6; 1; -Le/6];
    
    % DOF indices
    dof = [2*i-1, 2*i, 2*(i+1)-1, 2*(i+1)];
    
    % assemble
    K(dof,dof) = K(dof,dof) + ke;
    F(dof) = F(dof) + fe;
end

%%--- Apply boundary conditions (supports: v = 0) ---
fixedDOF = [];
for i = 1:length(supportPositions)
    pos = supportPositions(i);
    
    % find closest node
    [~, idx] = min(abs(X - pos));
    
    fixedDOF = [fixedDOF, 2*idx-1]; % displacement DOF only
end

allDOF = 1:ndof;
freeDOF = setdiff(allDOF, fixedDOF);

%%--- Solve system ---
U = zeros(ndof,1);
U(freeDOF) = K(freeDOF,freeDOF) \ F(freeDOF);

%%--- Compute reactions ---
R = K*U - F;

%%--- Extract reactions at supports (vertical forces only) ---
reactionForces = zeros(size(supportPositions));
for i = 1:length(supportPositions)
    pos = supportPositions(i);
    [~, idx] = min(abs(X - pos));
    
    reactionForces(i) = R(2*idx-1);
end

end