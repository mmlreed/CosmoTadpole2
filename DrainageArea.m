function g = DrainageArea(p,g)

% compute upslope contributing area using the requested flow routing method
[g.A g.minima g.fl g.W] = eval(['mex' p.routing '(g.S,p.dy/p.dx,1*(~g.C + g.C2),p.bvec,p.flood)']);

% multiply by cell area
g.A = p.dx*p.dy*g.A;
