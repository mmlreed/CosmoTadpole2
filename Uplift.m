function g = Uplift(p,g)

% Uplift/subsidence relative to boundaries. Note that uplift occurs only where C == 1
M = g.C;
Uminus1 = g.U;
Sminus1 = g.S;
g.U = g.U + p.dt.*p.E.*g.C;
g.S = g.U + g.H;
g.S(M==0) = Sminus1(M==0);


% if lake(s) are turned on
if p.doLake
    lake_level = datasample(g.S(g.C2 == 1),1);
    g.C2(g.S <= lake_level) = 1; % this could be changed to g.U so that the lake overtakes the land when its level is the same
    g.S(g.C2 == 1) = lake_level;
    g.U(g.C2 == 1) = lake_level;
    g.H(g.C2 == 1) = 1e-20;
    g.Channels(g.C2 == 1) = 0;
end
end


