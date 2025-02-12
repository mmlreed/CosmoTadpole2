function [g,p] = SoilProd_FCI(g,p)
%
% Update soil thickness H, bedrock elevation B, mineral abundances X, and
% cosmogenic radionuclide concentrations in bedrock and soil according to 
% the perturbations to these quantities by soil production.
% This is part of a splitting method, whereby this function computes a
% portion of the temporal changes in H, B, and X, and other functions
% compute the other portions of those changes.


%% Update soil thickness H
% Use exact solution for change in soil thickness over an interval dt.
% This solution can be derived via separation of variables.

% Parameters

H1 = g.H;
H1(H1 > 6) = 6; % interpolant only valid to 6 m

g.FCI = interp2(p.T1, p.Hc, p.Ci, p.T, H1); % frost-cracking intensity (FCI) in terms of porosity increase per unit area? (mm3/mm2)



p.P0 = p.spe.*g.FCI; % conversion of FCI to maximum soil production rate -- p.spe (1/T)




alpha = p.alpha;
dt = p.dt;
rhos = p.rhos;
rhor = p.rhor;

% Grids
B = g.U;
H = g.H;
C = g.Channels;
C2 = g.C2;
M = g.C;
X = g.X;

% below, log() does like zeros, which can occur with the intro of temp dependence 


% p.P0(p.P0<=0.350)=0.350; 



p.P0(p.P0<=p.P0_orig)=p.P0_orig;

g.P0 = p.P0;

% soil production
Hnplus1 = (1/alpha).*log(alpha.*p.P0.*dt./rhos + exp(alpha.*H));

if sum(isnan(Hnplus1(:)))>0
    fprintf('SoilProd_FCI.m: NaN(s) detected in H after soil production\n')
end


% This is to limit extremely thin soils
Hnplus1(Hnplus1<1e-20)=1e-20;

% boundary condition: no soil is produced, and therefore no bedrock
% lowering due to soil production, at boundaries.

% Hnplus1(C==1) = 1e-20;
Hnplus1(C==1) = H(C==1);

Hnplus1(C2==1) = H(C2==1);
Hnplus1(M==0) = 1e-20;

%% Update bedrock elevation B
% Conservation of mass gives new bedrock surface elevation

B_minus1 = B;
B = B + (rhos/rhor)*(H - Hnplus1);

if sum(isnan(B(:)))>0
    fprintf('SoilProd_FCI.m: NaN(s) detected in U after soil production\n')
end

g.U = B;
g.U(C==1) = B_minus1(C==1); % process isolated to hillslopes
g.U(C2==1) = B_minus1(C2==1);
g.U(M==0) = B_minus1(M==0);
% g.U(g.fl == 1) = B_minus1(g.fl == 1);

% Maintain extremely low levels of soil in channels and boundaries


% bedrock lowering velocity for semi-Langrangian advection -- will implement as needed
BLV = (p.P0./rhor).*exp(-alpha.*Hnplus1);
%BLV(C==1)=E;
%BLV(C==2)=0;
%BLV(M==0)=E;

% if sum(isnan(BLV(:)))>0
    % fprintf('SoilProd.m: NaN(s) detected in BLV\n')
% end

g.BLV = BLV;

%% Update soil mineral abundances X
for i=1:p.nX
    
    Xtemp = X(:,:,i).*H./Hnplus1 + p.xr(i)*(Hnplus1-H)./Hnplus1;
    
    % Apply boundary conditions.
    Xtemp2 = X(:,:,i);
    Xtemp(C==1) = Xtemp2(C==1); % enforce boundary condition: no change
    Xtemp(C2==1) = Xtemp2(C2==1); % lake cells
    Xtemp(M==0) = Xtemp2(M==0); % outside modeling domain

    Xtemp(p.F==1)=Xtemp2(p.F==1); 
    % in X at boundaries
    Xtemp(Xtemp<0)=0;
    g.X(:,:,i) = Xtemp;
end

for i = 1:p.nX
    Xtemp = g.X(:,:,i);
    Xtemp(isnan(Xtemp(:))) = p.xr(i);
    g.X(:,:,i) = Xtemp;
end

X_ratio = sum(g.X,3);
g.X = g.X./X_ratio;

H = Hnplus1;

H(C2==1) = 1e-20;
H(M==0) = 1e-20;

g.H = H;

if sum(isnan(X(:)))>0
    fprintf('SoilProd_FCI.m: NaN(s) detected in X\n')
end