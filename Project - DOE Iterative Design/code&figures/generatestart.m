function X0 = generatestart(imagetype,Index)
%% generate the starting point
i         = sqrt(-1);
N11       = Index(2,1);
N22       = Index(2,2);
S1start   = Index(3,1);
S1end     = Index(3,2);
S2start   = Index(4,1);
S2end     = Index(4,2);

X0        = zeros(N11,N22);
if imagetype == 1
    X0(S1start:S1end,S2start:S2end)  ...
    = rand(length(S1start:S1end),length(S2start:S2end));
elseif imagetype == 2
    X0(S1start:S1end,S2start:S2end)  ...
    = exp(i*pi/2*rand(length(S1start:S1end),length(S2start:S2end)));
elseif imagetype == -2
    X0(S1start:S1end,S2start:S2end)  ...
    = exp(i*2*pi*rand(length(S1start:S1end),length(S2start:S2end)));
end
