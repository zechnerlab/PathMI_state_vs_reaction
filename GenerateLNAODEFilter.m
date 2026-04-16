%% GenerateLNAODEFilter(S, h, params, condL)
%  This function is used to generate the conditional moment equations based
%  on the linear noise approximation. Note that the function is adapted to
%  this particular system and does not apply to other systems in this form.
%  
%  Parameters
%  S: stoichiomtery matrix
%  h: rate functions (defined symbolically)
%  params: list of parameter names (strings)
%  condL: this is used to specificy of the ligand is also conditioned on
%  (condL = 1) or not (condL = 0)

function [X, h, A, BBT, F, BBT_Y] = GenerateLNAODEFilter(S, h, params, condL)
numStates = size(S, 2);

% Define symbolic system state:
% X1: ligand
% X2: downstream species (X in the paper)
% X3: inactive state
% X4,X5,...: active states
X = sym('X', [numStates, 1]);
h = sym(h(:));

% Calculate jacobian
A = simplify(jacobian(S'*h, X));

% Eliminate off state (index 3) based on conservation. This leads to additional fluxes in the matrix A
% because the off state is substitutated in the first on state (index 4) by the sum of the active states. 
% This elimination is numerically advantagous. 
AR = A;
AR(4, 4:end) = A(4, 4:end) - A(4, 3); %the elimination adds fluxes coming from other active states.

% Calculate noise tensor
BBT = S'*diag(h)*S;

% Calculate rate equation fluxes
F = S'*h;

% Here we select the components that are part of the "hidden" network. If
% we condition on the ligand (condL=1), then those are only the active
% states, so we have N dimensions. If we don't condition on the ligand
% (condL=0), then we have the ligand included as a hidden state, leading to
% N+1 dimensions.
stateIdx = 1:numStates;
hiddenStates = setdiff(stateIdx, [2, 3, condL]); %take out X, inactive state and ligand (if condL=1)
numHiddenStates = length(hiddenStates);
A_X = AR(hiddenStates, hiddenStates); %Select relevant part of A

% Vector H such that H*X gives the drift of the observation process Y
H = AR(2, hiddenStates);
BBT_X = S(:, hiddenStates)'*diag(h)*S(:, hiddenStates);

% Noise of the observation process
BBT_Y = S(:, 2)'*diag(h)*S(:, 2);
BBT_Y = subs(BBT_Y, 'beta', 0);


% Write expressions to file
fileName = sprintf('LNAMomentsODE_%d.m', condL);
fid = fopen(fileName, 'w');

fprintf(fid, 'function dN = LNAMomentsODE_%d(t, N', condL);

for k=1:length(params)
   fprintf(fid, ', %s', params{k});
end
fprintf(fid, ')\n');

for k=1:numStates
    fprintf(fid, 'X%d=N(%d);\n', k, k);
end

fprintf(fid, 'C = reshape(N(%d:end), %d, %d);\n\n', numStates+1, numHiddenStates, numHiddenStates);

for k=1:numStates
fprintf(fid, 'dX(%d) = %s;\n', k, char(F(k, :)));
end

fprintf(fid, '\n');

for k=1:numHiddenStates
   for j=1:numHiddenStates
    fprintf(fid, 'A_X(%d, %d) = %s;\n', k, j, char(A_X(k, j)));
   end
end

fprintf(fid, '\n');

for k=1:numHiddenStates
   for j=1:numHiddenStates
    fprintf(fid, 'BBT_X(%d, %d) = %s;\n', k, j, char(BBT_X(k, j)));
   end
end

fprintf(fid, '\n');


for k=1:size(BBT_Y, 1)
   for j=1:size(BBT_Y, 2)
    fprintf(fid, 'BBT_Y(%d, %d) = %s;\n', k, j, char(BBT_Y(k, j)));
   end
end

for k=1:size(H, 1)
    for j=1:size(H, 2)
    fprintf(fid, 'H(%d, %d) = %s;\n', k, j, char(H(k, j)));
    end
end
   
fprintf(fid, '\n');

% This is the formula for the covariance equation
fprintf(fid, 'dC = A_X*C + C*A_X'' + BBT_X - C*H''*1/BBT_Y*H*C;');
fprintf(fid, '\n');
fprintf(fid, 'dN = [dX(:); reshape(dC, %d, 1)];\n', numHiddenStates^2);

fclose(fid);