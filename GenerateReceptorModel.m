% GenerateReceptorModel(N, c1, c2, kOn, mu, gamma, beta)
% This function defines the stoichiometry of the multi-state receptor with
% N active states.
function [Pre, Post, c, X0, h] = GenerateReceptorModel(N, c1, c2, kOn, mu, gamma, beta)

X0 = [100   1  ]'; %ligand and downstream species X
X0 = [X0; zeros(N, 1); 10]; %Add receptor states

% set parameters
c = [
    c1
    c2
    beta
    kOn
    mu*ones(N, 1);
    gamma*ones(N, 1);
    ];

% define rate functions
h = {str2sym('c1')
    str2sym('c2*X1')
    str2sym('beta*X2')
    str2sym('X3*X1*kOn')
    };
% add rates for the switching between receptor states with rate mu
for i=1:N
    xstr = sprintf('mu*X%d', i+3);
    h{end+1} = str2sym(xstr);
end

% add production reactions of downstream species X for each active state
for i=1:N
    xstr = sprintf('gamma*X%d', i+3);
    h{end+1} = str2sym(xstr);
end


% define reactant and product matrices
Pre0 = [0 0 0 0 
        1 0 0 0 
        0 1 0 0
        1 0 1 0 
    ];

Post0 = [1 0 0 0
         0 0 0 0
         0 0 0 0
         1 0 0 1
    ];

% extend by parts corresponding to the N active states
Pre = zeros(4 + 2*N, 3 + N);
Post = zeros(4 + 2*N, 3 + N);

% switching reactions
for k=1:N-1
    Pre(4+k, 3+k) = 1;
    Post(4+k, 4+k) = 1;
end

Pre(4+N, end) = 1;
Post(4+N, 3) = 1;

% production reactions of downstream species X
for k=1:N
    Pre(4+N+k, 3+k) = 1;
    Post(4+N+k, 3+k) = 1;
    Post(4+N+k, 2) = 1;
end

% put everything together
Pre(1:4, 1:4) = Pre0;
Post(1:4, 1:4) = Post0;

end