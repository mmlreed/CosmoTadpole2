function g = SoilProd(g,p)
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
P = p.P0; 
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


% soil production
Hnplus1 = (1/alpha)*log(alpha*P*dt/rhos + exp(alpha*H));

if sum(isnan(Hnplus1(:)))>0
    fprintf('SoilProd.m: NaN(s) detected in Hnplus1 after soil production\n')
end

% sometimes needed during spin-up with high dt
Hnplus1(isnan(Hnplus1))=1e-20;

% This is to limit extremely thin soils
Hnplus1(Hnplus1<1e-20)=1e-20;

% boundary condition: no soil is produced, and therefore no bedrock
% lowering due to soil production, at boundaries.

% Hnplus1(C==1) = 1e-20;
Hnplus1(C==1) = H(C==1);

Hnplus1(C2==1) = H(C2==1);
Hnplus1(M==0) = 1e-20;

Hnplus1(p.catch==1) = 0.5;



%% Update bedrock elevation B
% Conservation of mass gives new bedrock surface elevation

B_minus1 = B;
B = B + (rhos/rhor)*(H - Hnplus1);

if sum(isnan(B(:)))>0
    fprintf('SoilProd.m: NaN(s) detected in U after soil production\n')
end

g.U = B;
g.U(C==1) = B_minus1(C==1); % process isolated to hillslopes
g.U(C2==1) = B_minus1(C2==1);
g.U(M==0) = B_minus1(M==0);

g.U(p.catch==1)=B_minus1(p.catch==1);

% Maintain extremely low levels of soil in channels and boundaries


% bedrock lowering velocity for semi-Langrangian advection -- will implement as needed
% BLV = (P/rhor)*exp(-alpha.*Hnplus1);
%BLV(C==1)=E;
%BLV(C==2)=0;
%BLV(M==0)=E;
% 
% % if sum(isnan(BLV(:)))>0
%     % fprintf('SoilProd.m: NaN(s) detected in BLV\n')
% % end
% 
% g.BLV = BLV;

%% Update soil mineral abundances X
for i=1:p.nX
    Xtemp = X(:,:,i).*H./Hnplus1 + p.xr(i)*(Hnplus1-H)./Hnplus1;
    % Apply boundary conditions.
    Xtemp2 = X(:,:,i);
    Xtemp(C==1) = Xtemp2(C==1); % enforce boundary condition: no change
    Xtemp(C2==1) = Xtemp2(C2==1);
    Xtemp(M==0) = Xtemp2(M==0);
    Xtemp(p.F==1)=Xtemp2(p.F==1);
    Xtemp(p.catch==1)=p.xr(i);
    g.X(:,:,i) = Xtemp;
    g.X(g.X<=0)=p.xr(i);
end

% During spin-up this can keep things turning over

for i = 1:p.nX
    Xtemp = g.X(:,:,i);
    Xtemp(isnan(Xtemp(:))) = p.xr(i);
    g.X(:,:,i) = Xtemp;
end

X_ratio = sum(g.X,3);
g.X = g.X./X_ratio;

H = Hnplus1;

H(C==1) = 1e-20;
H(C2==1) = 1e-20;
H(M==0) = 1e-20;



g.H = H;

g.H(g.H>10)=10;
g.H(isinf(g.H))=10;

g.P0=p.P0.*ones(size(g.H));

if sum(isnan(g.X(:)))>0
    fprintf('SoilProd.m: NaN(s) detected in g.X\n')
end

if sum(isnan(g.U(:)))>0
    fprintf('SoilProd.m: NaN(s) detected in g.U\n')
end
if sum(isnan(g.S(:)))>0
    fprintf('SoilProd.m: NaN(s) detected in g.S\n')
end
if sum(isnan(g.H(:)))>0
    fprintf('SoilProd.m: NaN(s) detected in g.H\n')
end
