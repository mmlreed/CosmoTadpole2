function g = NL8(g,p) 

% Setting up parameters
del_x = p.dx;
del_y = p.dy; 
dx = p.dx; % grid resolution X
dy = p.dy; % grid resolution Y
dt = p.dt; % timestep (yr); ignore Courant for now
S_c = p.Sc; % critical hillslope gradient
K = p.K; % number of cells in X dir
J = p.J; % number of cells in Y dir
K3 = p.D; % non-linear diffusion-like constant; Roering (2008) notation

% Grids
S = g.S; % surface elevation
C = g.Channels; % marked channels
C2 = g.C2; % lake cells
M = g.C; % boundaries; this is confusing; fix later
F = p.F; % areas of no change/ghost boundary
H = g.H; % soil thickness
X = g.X; % mineral concentrations
Ns = g.Ns; % nuclide concentration
rhos = p.rhos;
rhor = p.rhor;
host_min = p.host_min;
nX = p.nX;

% Derived parameters
diag_k = sqrt(dx^2 + dy^2);
diag_delk = (2*dx*dy)/diag_k;


if sum(isnan(H(:)))>0
    fprintf('NL8.m: NaN(s) detected in H before diffusion!\n')
end

if sum(isnan(S(:)))>0
    fprintf('NL8.m: NaN(s) detected in S before diffusion!\n')
end

% save previous 
Snminus_1 = S;

% Let's call these right way from now on... (x,y) - (row,col)
% mirror y direction
% top side
S_y_L = [S(:,2) S];
% bottom side
S_y_R = [S S(:,end-1)];

% mirror x direction
% left side
S_x_T = [S(2,:); S]; % S(end-1,:)]; 
% right side
S_x_B = [S;  S(end-1,:)];


% mirrored for easier diagonals
S_diag = [S(:,2) S S(:,end-1)];
S_diag = [S_diag(2,:); S_diag; S_diag(end-1,:)]; 

S_xy_TL = S_diag(1:end-1, 1:end-1);
S_xy_BL = S_diag(2:end, 1:end-1);
S_xy_BR = S_diag(2:end, 2:end);
S_xy_TR = S_diag(1:end-1, 2:end);


% 
% Need to compute spatial 1st derivatives for all eight directions with one-way
% differencing (z_ij - z_k)/dk going in counter-clockwise direction starting 
% from directly above z_ij. This choice is abritrary

S_k1 = (S_x_T(2:end, 1:end) - S_x_T(1:end-1, 1:end))/dx;

S_k2 = (S_xy_TL(2:end, 2:end) - S_xy_TL(1:end-1, 1:end-1))/diag_k; 

S_k3 = (S_y_L(1:end, 2:end) - S_y_L(1:end, 1:end-1))/dy;

S_k4 = (S_xy_BL(1:end-1,2:end) - S_xy_BL(2:end, 1:end-1))/diag_k;

S_k5 = (S_x_B(1:end-1, 1:end) - S_x_B(2:end, 1:end))/dx;

S_k6 = (S_xy_BR(1:end-1,1:end-1) - S_xy_BR(2:end, 2:end))/diag_k;

S_k7 = (S_y_R(1:end, 1:end-1) - S_y_R(1:end, 2:end))/dy;

S_k8 = (S_xy_TR(2:end, 1:end-1) - S_xy_TR(1:end-1,2:end))/diag_k;



% Slopes at Sc or above will produce infinite fluxes leading to NaNs 

S_k1(S_k1>=S_c)=S_c - 0.01;
S_k1(S_k1<=-S_c)=-S_c + 0.01;

S_k2(S_k2>=S_c)=S_c - 0.01;
S_k2(S_k2<=-S_c)=-S_c + 0.01;

S_k3(S_k3>=S_c)=S_c - 0.01;
S_k3(S_k3<=-S_c)=-S_c + 0.01;

S_k4(S_k4>=S_c)=S_c - 0.01;
S_k4(S_k4<=-S_c)=-S_c + 0.01;

S_k5(S_k5>=S_c)=S_c - 0.01;
S_k5(S_k5<=-S_c)=-S_c + 0.01;

S_k6(S_k6>=S_c)=S_c - 0.01;
S_k6(S_k6<=-S_c)=-S_c + 0.01;

S_k7(S_k7>=S_c)=S_c - 0.01;
S_k7(S_k7<=-S_c)=-S_c + 0.01;

S_k8(S_k8>=S_c)=S_c - 0.01;
S_k8(S_k8<=-S_c)=-S_c + 0.01;

% The following algorithm allows for the setting of lake-side facing cell
% slopes to be a certain value (rough simulation of erosive waves) -- Slow
% model down

% if p.t >= p.lake_time
% 
%     [m, n] = size(g.C2);
%     
%     % Initialize a cell array to store neighbors for each cell
%     allNeighbors = cell(m, n);
%     
%     % Iterate over each cell
%     for i = 1:m
%         for j = 1:n
%             % Check neighbors for the current cell
%             neighbors = g.C2(max(1, i-1):min(m, i+1), max(1, j-1):min(n, j+1));
%             
%             % Store neighbors in the cell array
%             allNeighbors{i, j} = neighbors;
%         end
%     end
%     
%     
%     S_k1_abs = abs(S_k1);
%     S_k2_abs = abs(S_k2);
%     S_k3_abs = abs(S_k3);
%     S_k4_abs = abs(S_k4);
%     S_k5_abs = abs(S_k5);
%     S_k6_abs = abs(S_k6);
%     S_k7_abs = abs(S_k7);
%     S_k8_abs = abs(S_k8);
%     
%     S_k1_sign = zeros(size(S_k1));
%     S_k1_sign(S_k1 > 0) = 1;
%     S_k1_sign(S_k1 < 0) = -1;
%     
%     S_k2_sign = zeros(size(S_k2));
%     S_k2_sign(S_k2 > 0) = 1;
%     S_k2_sign(S_k2 < 0) = -1;
%     
%     S_k3_sign = zeros(size(S_k3));
%     S_k3_sign(S_k3 > 0) = 1;
%     S_k3_sign(S_k3 < 0) = -1;
%     
%     S_k4_sign = zeros(size(S_k4));
%     S_k4_sign(S_k4 > 0) = 1;
%     S_k4_sign(S_k4 < 0) = -1;
%     
%     S_k5_sign = zeros(size(S_k5));
%     S_k5_sign(S_k5 > 0) = 1;
%     S_k5_sign(S_k5 < 0) = -1;
%     
%     S_k6_sign = zeros(size(S_k6));
%     S_k6_sign(S_k6 > 0) = 1;
%     S_k6_sign(S_k6 < 0) = -1;
%     
%     S_k7_sign = zeros(size(S_k7));
%     S_k7_sign(S_k7 > 0) = 1;
%     S_k7_sign(S_k7 < 0) = -1;
%     
%     S_k8_sign = zeros(size(S_k8));
%     S_k8_sign(S_k8 > 0) = 1;
%     S_k8_sign(S_k8 < 0) = -1;
%     
%     
%     
%     for i = 2:m-1
%         for j = 2:n-1
%             neighbs_lake = cell2mat(allNeighbors(i,j));
%             if neighbs_lake(1,1) == 1 & S_k8_abs(i,j) < p.lake_slope
%                 S_k8(i,j) = p.lake_slope*S_k8_sign(i,j);
%             end
%             if neighbs_lake(1,2) == 1 & S_k1_abs(i,j) < p.lake_slope
%                 S_k1(i,j) = p.lake_slope*S_k1_sign(i,j);
%             end
%             if neighbs_lake(1,3) == 1 & S_k2_abs(i,j) < p.lake_slope
%                 S_k2(i,j) = p.lake_slope*S_k2_sign(i,j);
%             end
%             if neighbs_lake(2,1) == 1 & S_k7_abs(i,j) < p.lake_slope
%                 S_k7(i,j) = p.lake_slope*S_k7_sign(i,j);
%             end
%             if neighbs_lake(2,3) == 1 & S_k3_abs(i,j) < p.lake_slope
%                 S_k3(i,j) = p.lake_slope*S_k3_sign(i,j);
%             end
%             if neighbs_lake(3,1) == 1 & S_k6_abs(i,j) < p.lake_slope
%                 S_k6(i,j) = p.lake_slope*S_k6_sign(i,j);
%             end
%             if neighbs_lake(3,2) == 1 & S_k5_abs(i,j) < p.lake_slope
%                 S_k5(i,j) = p.lake_slope*S_k5_sign(i,j);
%             end
%             if neighbs_lake(3,3) == 1 & S_k4_abs(i,j) < p.lake_slope
%                 S_k4(i,j) = p.lake_slope*S_k4_sign(i,j);
%             end
%         end
%     end
% end


% calculate vector flux qk for each direction using nonlinear diffusion
q_k1 = (-K3.*S_k1)./(1-(S_k1./S_c).^2);
q_k2 = (-K3.*S_k2)./(1-(S_k2./S_c).^2);
q_k3 = (-K3.*S_k3)./(1-(S_k3./S_c).^2);
q_k4 = (-K3.*S_k4)./(1-(S_k4./S_c).^2);
q_k5 = (-K3.*S_k5)./(1-(S_k5./S_c).^2);
q_k6 = (-K3.*S_k6)./(1-(S_k6./S_c).^2);
q_k7 = (-K3.*S_k7)./(1-(S_k7./S_c).^2);
q_k8 = (-K3.*S_k8)./(1-(S_k8./S_c).^2);


% uncomment for linear gradient relationship as per Andersen et al. (2015)
% Does not change much in terms of Dinf response
% if p.t >= p.frost_time 
%     q_k1 = (-K3.*S_k1);
%     q_k2 = (-K3.*S_k2);
%     q_k3 = (-K3.*S_k3);
%     q_k4 = (-K3.*S_k4);
%     q_k5 = (-K3.*S_k5);
%     q_k6 = (-K3.*S_k6);
%     q_k7 = (-K3.*S_k7);
%     q_k8 = (-K3.*S_k8);
% end
    
% No flux within streams; flux to streams 'transported' away
% q_k1(C==1)=0;
% q_k2(C==1)=0;
% q_k3(C==1)=0;
% q_k4(C==1)=0;
% q_k5(C==1)=0;
% q_k6(C==1)=0;
% q_k7(C==1)=0;
% q_k8(C==1)=0;

% flux into lake -- turned on for calculating lake flux but can be turned
% off here
% q_k1(C2==1)=0;
% q_k2(C2==1)=0;
% q_k3(C2==1)=0;
% q_k4(C2==1)=0;
% q_k5(C2==1)=0;
% q_k6(C2==1)=0;
% q_k7(C2==1)=0;
% q_k8(C2==1)=0;

% No flux outside of watershed
q_k1(M==0)=0;
q_k2(M==0)=0;
q_k3(M==0)=0;
q_k4(M==0)=0;
q_k5(M==0)=0;
q_k6(M==0)=0;
q_k7(M==0)=0;
q_k8(M==0)=0;

% No change in masked areas/buffer/ghost boundary
q_k1(F==1)=0;
q_k2(F==1)=0;
q_k3(F==1)=0;
q_k4(F==1)=0;
q_k5(F==1)=0;
q_k6(F==1)=0;
q_k7(F==1)=0;
q_k8(F==1)=0;

% make shifted grids for p.F (mask) for no-flow boundary to speed up computation
pF_k1 = circshift(F, [1 0]);
pF_k2 = circshift(F, [1 1]);
pF_k3 = circshift(F, [0 1]);
pF_k4 = circshift(F, [-1 1]);
pF_k5 = circshift(F, [-1 0]);
pF_k6 = circshift(F, [-1 -1]);
pF_k7 = circshift(F, [0 -1]);
pF_k8 = circshift(F, [1 -1]);

% pF matrix sizes do not match q_k
q_k1(pF_k1 == 1) = 0;
q_k2(pF_k2 == 1) = 0;
q_k3(pF_k3 == 1) = 0;
q_k4(pF_k4 == 1) = 0;
q_k5(pF_k5 == 1) = 0;
q_k6(pF_k6 == 1) = 0;
q_k7(pF_k7 == 1) = 0;
q_k8(pF_k8 == 1) = 0;

% turn these into Qk's; Qk = del_k*qk
Q_k1 = q_k1*dy;
Q_k2 = q_k2*diag_delk;
Q_k3 = q_k3*dx;
Q_k4 = q_k4*diag_delk;
Q_k5 = q_k5*dy;
Q_k6 = q_k6*diag_delk;
Q_k7 = q_k7*dx;
Q_k8 = q_k8*diag_delk;

% Change in elevation is the sum of Qk's divided by area
dz_dt = (1/(dx*dy))*(Q_k1 + Q_k2 + Q_k3 + Q_k4 + Q_k5 + Q_k6 + Q_k7 + Q_k8)*dt;

% just ensuring that zeroing out works
dz_dt(C==1)=0;
dz_dt(C2==1)=0;
dz_dt(M==0)=0;
dz_dt(F==1)=0;

dz_dt_orig = dz_dt;


if sum(isnan(dz_dt(:)))>0
    fprintf('NL8.m: NaN(s) detected in dz_dt\n')
end

% Correct for high outgoing fluxes -- This does "lock in" very low soil
% depths in areas with slopes near the critical slope (a feature not a
% bug?)..

qs_ratio=ones(size(S));


neg_H = -1*H;

dz_dt1 = dz_dt;

dz_dt1(dz_dt <= neg_H) = neg_H(dz_dt <= neg_H);

dz_dt = dz_dt1;
% 
qs_ratio = dz_dt_orig./dz_dt;

qs_ratio(isnan(qs_ratio))=1;

% % correct fluxes for high outgoing fluxes -- qs_ratio is 1 in cells that do not need corrected, eliminating the loops below
q_k1(q_k1 < 0) = q_k1(q_k1 < 0)./qs_ratio(q_k1 < 0);
q_k2(q_k2 < 0) = q_k2(q_k2 < 0)./qs_ratio(q_k2 < 0);
q_k3(q_k3 < 0) = q_k3(q_k3 < 0)./qs_ratio(q_k3 < 0);
q_k4(q_k4 < 0) = q_k4(q_k4 < 0)./qs_ratio(q_k4 < 0);
q_k5(q_k5 < 0) = q_k5(q_k5 < 0)./qs_ratio(q_k5 < 0);
q_k6(q_k6 < 0) = q_k6(q_k6 < 0)./qs_ratio(q_k6 < 0);
q_k7(q_k7 < 0) = q_k7(q_k7 < 0)./qs_ratio(q_k7 < 0);
q_k8(q_k8 < 0) = q_k8(q_k8 < 0)./qs_ratio(q_k8 < 0);

 
% Now deposition fluxes must be recalculated
% Loop through everything again to correct for high incoming fluxes

q_k1_p = circshift(-q_k5, [1 0]);
q_k2_p = circshift(-q_k6, [1 1]);
q_k3_p = circshift(-q_k7, [0 1]);
q_k4_p = circshift(-q_k8, [-1 1]);
q_k5_p = circshift(-q_k1, [-1 0]);
q_k6_p = circshift(-q_k2, [-1 -1]);
q_k7_p = circshift(-q_k3, [0 -1]);
q_k8_p = circshift(-q_k4, [1 -1]);

q_k1(q_k1 >= 0) = q_k1_p(q_k1 >= 0);
q_k2(q_k2 >= 0) = q_k2_p(q_k2 >= 0);
q_k3(q_k3 >= 0) = q_k3_p(q_k3 >= 0);
q_k4(q_k4 >= 0) = q_k4_p(q_k4 >= 0);
q_k5(q_k5 >= 0) = q_k5_p(q_k5 >= 0);
q_k6(q_k6 >= 0) = q_k6_p(q_k6 >= 0);
q_k7(q_k7 >= 0) = q_k7_p(q_k7 >= 0);
q_k8(q_k8 >= 0) = q_k8_p(q_k8 >= 0);


Q_k1 = q_k1*dy;
Q_k2 = q_k2*diag_delk;
Q_k3 = q_k3*dx;
Q_k4 = q_k4*diag_delk;
Q_k5 = q_k5*dy;
Q_k6 = q_k6*diag_delk;
Q_k7 = q_k7*dx;
Q_k8 = q_k8*diag_delk;

% total volumetric flux -- useful later for catchment/lake-averaged
% denudation rates
Q_k_tot = (Q_k1 + Q_k2 + Q_k3 + Q_k4 + Q_k5 + Q_k6 + Q_k7 + Q_k8);

g.Q_k_tot = Q_k_tot;

dz_dt = (1/(dx*dy))*(Q_k1 + Q_k2 + Q_k3 + Q_k4 + Q_k5 + Q_k6 + Q_k7 + Q_k8)*dt;


% No change zones
dz_dt(C==1) = 0;
dz_dt(M==0) = 0;
dz_dt(C2==1) = 0;
dz_dt(F==1) = 0;


% disp(max(dz_dt(:)))


S = Snminus_1 + dz_dt;

S(isnan(S)) = Snminus_1(isnan(S));

g.S = S;

if sum(isnan(S(:)))>0
    fprintf('NL_NPS.m: NaN(s) detected in S\n')
end

% volumetric flux to sinks for different post-processing calculations


% Calc changes in X due to transport in terms of mass (+ -> gain; - ->
% loss)


qx_k1 = q_k1.*dx.*rhos.*X;
qx_k2 = q_k2.*diag_k.*rhos.*X;
qx_k3 = q_k3.*dy.*rhos.*X;
qx_k4 = q_k4.*diag_k.*rhos.*X;
qx_k5 = q_k5.*dx.*rhos.*X;
qx_k6 = q_k6.*diag_k.*rhos.*X;
qx_k7 = q_k7.*dy.*rhos.*X;
qx_k8 = q_k8.*diag_k.*rhos.*X;


% positive fluxes require usage of data from point of origin instead of local; allocate these at there are threeeee lewps


qx_k1_p = q_k1.*dx.*rhos.*circshift(X, [1 0]);
qx_k2_p = q_k2.*diag_k.*rhos.*circshift(X, [1 1]);
qx_k3_p = q_k3.*dy.*rhos.*circshift(X, [0 1]);
qx_k4_p = q_k4.*diag_k.*rhos.*circshift(X, [-1 1]);
qx_k5_p = q_k5.*dx.*rhos.*circshift(X, [-1 0]);
qx_k6_p = q_k6.*diag_k.*rhos.*circshift(X, [-1 -1]);
qx_k7_p = q_k7.*dy.*rhos.*circshift(X, [0 -1]);
qx_k8_p = q_k8.*diag_k.*rhos.*circshift(X, [1 -1]);



q_k1_r = repmat(q_k1, [1,1,p.nX]);
q_k2_r = repmat(q_k2, [1,1,p.nX]);
q_k3_r = repmat(q_k3, [1,1,p.nX]);
q_k4_r = repmat(q_k4, [1,1,p.nX]);
q_k5_r = repmat(q_k5, [1,1,p.nX]);
q_k6_r = repmat(q_k6, [1,1,p.nX]);
q_k7_r = repmat(q_k7, [1,1,p.nX]);
q_k8_r = repmat(q_k8, [1,1,p.nX]);

qx_k1(q_k1_r>0)=qx_k1_p(q_k1_r>0);
qx_k2(q_k2_r>0)=qx_k2_p(q_k2_r>0);
qx_k3(q_k3_r>0)=qx_k3_p(q_k3_r>0);
qx_k4(q_k4_r>0)=qx_k4_p(q_k4_r>0);
qx_k5(q_k5_r>0)=qx_k5_p(q_k5_r>0);
qx_k6(q_k6_r>0)=qx_k6_p(q_k6_r>0);
qx_k7(q_k7_r>0)=qx_k7_p(q_k7_r>0);
qx_k8(q_k8_r>0)=qx_k8_p(q_k8_r>0);


% Total mass flux in minerals
qx_tot = (qx_k1 + qx_k2 + qx_k3 + qx_k4 + qx_k5 + qx_k6 + qx_k7 + qx_k8);

% Insert safegard conservation of X here (for very, very low soil thicknesses qx_tot*dt can be larger than mass of X, leading to negative concentrations (bad)) 

qx_tot_dt = qx_tot.*p.dt;
qx_filter = ones(3);

for i=1:p.nX
    qx_tot_dt_1 = qx_tot_dt(:,:,i);
    qx_mass = -1*g.H*p.dx*p.dy*p.rhos.*X(:,:,i);
    qx_tot_dt_1(qx_tot_dt_1 <= qx_mass) = qx_mass(qx_tot_dt_1 <= qx_mass);
    qx_ratio = qx_tot_dt(:,:,i)./qx_tot_dt_1;
    qx_neighbors = conv2(-qx_mass, qx_filter, 'same');
	qx_tot_dt_1(qx_tot_dt_1 >= qx_neighbors) = qx_neighbors(qx_tot_dt_1 >= qx_neighbors);
    qx_tot_dt(:,:,i) = qx_tot_dt_1;
end

qx_tot = qx_tot_dt./p.dt;

% Get fluxes to lake sink before zeroing out 
g.qx_tot_lake = qx_tot;

% No flux zones
qx_tot(C==1)=0;
qx_tot(C2==1)=0;
qx_tot(M==0)=0;
qx_tot(F==1)=0;


g.qx_tot = qx_tot;

% Flux of 10Be atoms in all 8 direction; Quartz is fluxed using the same
% volumetric flux, so these can be modeled separately and still be
% conserved

qns_k1 = q_k1.*dx.*rhos.*X(:,:,host_min).*Ns;
qns_k2 = q_k2.*diag_k.*X(:,:,host_min).*rhos.*Ns;
qns_k3 = q_k3.*dy.*X(:,:,host_min).*rhos.*Ns;
qns_k4 = q_k4.*diag_k.*rhos.*X(:,:,host_min).*Ns;
qns_k5 = q_k5.*dx.*rhos.*X(:,:,host_min).*Ns;
qns_k6 = q_k6.*diag_k.*rhos.*X(:,:,host_min).*Ns;
qns_k7 = q_k7.*dy.*rhos.*X(:,:,host_min).*Ns;
qns_k8 = q_k8.*diag_k.*X(:,:,host_min).*rhos.*Ns;


% Correct for positive fluxes into grid cell

X_host = X(:,:,p.host_min);

qns_k1_p = q_k1.*dx.*rhos.*circshift(X_host, [1 0]).*circshift(Ns, [1 0]);
qns_k2_p = q_k2.*diag_k.*rhos.*circshift(X_host, [1 1]).*circshift(Ns, [1 1]);
qns_k3_p = q_k3.*dy.*rhos.*circshift(X_host, [0 1]).*circshift(Ns, [0 1]);
qns_k4_p = q_k4.*diag_k.*rhos.*circshift(X_host, [-1 1]).*circshift(Ns, [-1 1]);
qns_k5_p = q_k5.*dx.*rhos.*circshift(X_host, [-1 0]).*circshift(Ns, [-1 0]);
qns_k6_p = q_k6.*diag_k.*rhos.*circshift(X_host, [-1 -1]).*circshift(Ns, [-1 -1]);
qns_k7_p = q_k7.*dy.*rhos.*circshift(X_host, [0 -1]).*circshift(Ns, [0 -1]);
qns_k8_p = q_k8.*diag_k.*rhos.*circshift(X_host, [1 -1]).*circshift(Ns, [1 -1]);



qns_k1(q_k1>0)=qns_k1_p(q_k1>0);
qns_k2(q_k2>0)=qns_k2_p(q_k2>0);
qns_k3(q_k3>0)=qns_k3_p(q_k3>0);
qns_k4(q_k4>0)=qns_k4_p(q_k4>0);
qns_k5(q_k5>0)=qns_k5_p(q_k5>0);
qns_k6(q_k6>0)=qns_k6_p(q_k6>0);
qns_k7(q_k7>0)=qns_k7_p(q_k7>0);
qns_k8(q_k8>0)=qns_k8_p(q_k8>0);

qns_tot = (qns_k1 + qns_k2 + qns_k3 + qns_k4 + qns_k5 + qns_k6 + qns_k7 + qns_k8);

% threshold = 1e10;
% [maxNs, maxNs_ind] = max(qns_tot(:));
% if maxNs > threshold
%     disp('qns_tot has large value(s) (>1e10) after transport but before flux correction')
%     disp(maxNs)
%     [a,b] = ind2sub(size(qns_tot), maxNs_ind);
%     disp([a,b])
%     disp(g.H(a,b))
% end


% safegard conservation of Ns (eliminates unreasonably high
% fluxes) which are not eliminated by correcting the mass flux based on
% soil thickness availability (similar to above.. This doesn't happen but just to be safe)

qns_tot_dt = qns_tot.*p.dt;

qns_atoms = -Ns.*H*p.dx*p.dy*p.rhos.*X(:,:,p.host_min);

qns_tot_dt(qns_tot_dt < qns_atoms) = qns_atoms(qns_tot_dt < qns_atoms);

qns_filter = ones(3);
qns_atoms_neighbors = conv2(-qns_atoms, qns_filter, 'same');

qns_tot_dt(qns_tot_dt > qns_atoms_neighbors) = qns_atoms_neighbors(qns_tot_dt > qns_atoms_neighbors);

qns_tot = qns_tot_dt./p.dt;

g.qns_tot_lake = qns_tot;

% No flux zones
qns_tot(C==1)=0;
qns_tot(C2==1)=0;
qns_tot(M==0)=0;
qns_tot(F==1)=0;

threshold = 1e20;
[maxNs,maxNs_ind] = max(abs(qns_tot(:)));
if maxNs > threshold
    disp('qns_tot has large value (>1e20) after transport but before Erode.m')
    disp(qns_tot(maxNs_ind))
    [a,b] = ind2sub(size(qns_tot), maxNs_ind);
    disp([a,b])
    disp(g.H(a,b))
    disp(g.Q_k_tot(a,b))
end


g.qns_tot = qns_tot;



end