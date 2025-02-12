function [p g] = FT(p,g)

% calculate RHS of the NKWE, recording the Courant number and which points 
% exceed the channel incision threshold
[RHS g Courant] = NKWE(g,p);

% adjust the time step based on the specified maximum Courant number and
% maximum time step
if p.doAdaptiveTimeStep
    p.dt = min( [p.dtmax p.Courant/Courant*p.dt]);
end

% Here we would have use the something like g.S where g.S is g.U + g.H (a
% soil thickness grid) -- then update the g.H based on integrated erosion 
C = g.Channels;
M = g.C;

% eroding soil
Sminus1 = g.S;
Hminus1 = g.H;
g.S = g.S + p.dt*RHS;
g.H = g.S - g.U;

if sum(isnan(g.H(:)))>0
    fprintf('FT.m: NaN(s) detected in H\n')
end

g.fluv_erosion = RHS;


Uminus1 = g.U;

% assures that eroding into rock is harder than into soil when g.H goes
% negative (e.g., both soil and bedrock were eroded)
neg_H1 = g.H(g.H < 0);

g.H(g.H < 0) = g.H(g.H < 0)./(p.Kf_s/p.Kf_r);

neg_H2 = g.H(g.H < 0);

neg_H_diff = abs(neg_H1 - neg_H2);

g.S(g.H < 0) = g.S(g.H <0) - neg_H_diff;

g.U(g.H<=0) = g.S(g.H<=0);

% These can go back to not being 'channels', so soil and transport can
% occur -- or g.Channels can represent anything beyond the threshold
% with no soil production and transport possible

S1 = mexSlope(g.S,p.dx,p.dy,p.bvec);

g.Channels(g.H > 0.002) = 0;
g.Channels(g.A >= (p.dx*p.dy)*p.Acr2) = 1;
g.Channels(S1 >= p.Sc) = 0;

C = g.Channels;

% maintain the infinitescimal soil thickness in bedrock channels 

g.H(C==1 & g.H<=0) = 1e-20;

g.H(g.C2==1) = 1e-20;
g.H(M==0) = 1e-20;
g.H(p.F==1) = 1e-20;


g.S = g.U + g.H;

if sum(isnan(g.S(:)))>0
    fprintf('FT.m: NaN(s) detected in S at end of calc\n')
end

if sum(isinf(g.S(:)))>0
    fprintf('FT.m: Infs(s) detected in S at end of calc\n')
end

% for calculating basin-averaged/lake-averaged Dinf later
% volumetric flux of sediment (easier to work with) -- no soil flux from
% bedrock channels and vice versa



g.Ns(C==1 & g.H <= 0.002) = g.Nzb(C==1 & g.H <= 0.002);
