function [RHS g Courant] = NKWE(g,p)

% evaluates the RHS of the nonlinear kinematic wave term


% calculate drainage area

g = DrainageArea(p,g);

% calculate slope
S1 = mexSlope(g.S,p.dx,p.dy,p.bvec);

[g S] = UpwindSlope(p,g);
S(S<0) = 0; % This will only happen where we've routed flow uphill to flood depressions


% calculate Keff, which accounts for channel width < grid spacing
if p.wexp == 0
    Keff_s = p.Kf_s * p.Kw / p.dx;
    Keff_r = p.Kf_r * p.Kw / p.dx;
elseif p.wexp == 1
    Keff_r = p.Kf_r * p.Kw * g.A / p.dx;
    Keff_s = p.Kf_s * p.Kw * g.A / p.dx;
else
    Keff_r = p.Kf_r * p.Kw * g.A.^p.wexp ./ p.dx;
    Keff_s = p.Kf_s * p.Kw * g.A.^p.wexp ./ p.dx;
end

Am = g.A.^p.m;

if p.w ~= 1 % ditto for slope   
    Sw = S.^p.w;
else
    Sw = S;
end

AmSw = Am.*Sw;


% record locations of channels (points that exceed threshold)
Channels = AmSw > p.thetac;

% No channels below 2 grid cells -- prevents artifactual 'spires' that
% cannot erode
Channels(g.A < (p.dx*p.dy)*2) = 0;

if p.Acr == p.Acr2
    Channels(g.A < (p.dx*p.dy)*p.Acr) = 0;
end

% Channels cannot be greater than the threshold slope
Channels(S1 >= p.Sc) = 0;


g.Channels = Channels;

% calculate right-hand side of NKWE 

Keff = zeros(size(g.Channels));
Keff(g.H > 1e-20) = Keff_s; 
Keff(g.H == 1e-20) = Keff_r; 

if p.thetac == 0
        
    RHS = -Keff .* AmSw;    
    
else
    
    % Case 1: if erosion goes as excess shear stress (conventional):
    RHS = -Keff .* (AmSw - p.thetac);
    
%     % Case 2: if erosion goes as shear stress:
%     RHS = -Keff .* AmSw;

    
    % In either case, this term is only nonzero where AmSw > thetac (bedrock can only erode)
    RHS = RHS.*(RHS < 0); 
    
end



% no-erosion rules 

% no erosion in:
%
% 1. areas that are flooded (or local minima, if the user did not request
% to flood depressions when doing flow routing). We assume these are areas of
% temporary deposition.
%
% 2. Fixed boundaries and other locations requested by user, as recorded in the matrix C.
%
% 3. Non-channels. Note that if we use a midpoint method like RK2, we 
% should use the channel grid determined at the beginning of
% the time step, not the grid identified locally, since this call to NKWE
% could be one of the sub-steps in the RK2 method.

% if p.flood==1
%     RHS = RHS.*(~fl & g.C==1 & g.Channels); 
% else
%     RHS = RHS.*(~min & g.C==1 & g.Channels);
% end

if p.flood==1
    RHS = RHS.*(~g.fl & g.C==1 & Channels); 
else
    RHS = RHS.*(~g.minima & g.C==1 & Channels);
end

% No fluvial incision in lake cells
RHS(g.C2==1)=0;

% maximum Courant #, neglecting threshold:
% if dz/dt = -Keff*A^m*S^n
% then C = Keff*A^m*S^(n-1)*dt/dx

Courant = max(Keff(:).*AmSw(:)./S(:))*sqrt(2)*p.dt/p.dx;
