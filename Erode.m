function g = Erode(p,g)
%
% Update soil thickness H, bedrock elevation B, and mineral abundances X
% according to the perturbations to these quantities by soil erosion.
% This is part of a splitting method, whereby this function computes a
% portion of the temporal changes in H, B, and X, and other functions
% compute the other portions of those changes.


% NONLINEAR DIFFUSION  - calculate surface elevation
g = NL8(g, p);

g.S(g.U>g.S)=g.U(g.U>g.S);

% g.S(p.catch==1)=p.catch_elev(p.catch==1);

Hnplus1 = g.S - g.U; % calculate new soil thickness
Hnplus1(Hnplus1==0)=1e-20;

Hnplus1(isinf(Hnplus1))=10;
Hnplus1(Hnplus1>10)=10;

% Boundary conditions and the condition of no soil or extremely low soil
% thickness

C = g.Channels;
C2 = g.C2;
M = g.C;
F = p.F;

Hnplus1(Hnplus1<1e-20)=1e-20;
Hnplus1(M==0)=1e-20;
Hnplus1(C==1)=1e-20;
Hnplus1(C2==1)=1e-20;
Hnplus1(p.F==1)=1e-20;

Hnplus1(p.catch==1)=0.5;

if sum(isnan(Hnplus1(:)))>0
    fprintf('Erode.m: NaN(s) detected in Hnplus1 after transport\n')
end

%% Update mineral abundances X

% Use vector flux mineral flux to calculate change in concentration

X_old=g.X;

% Account for fluxes that want to zero out concentration (only occurs @
% very low soil thickness (<1 mm). 

% g.Q_k_tot(C==1)=0;
g.Q_k_tot(M==0)=0;
% g.Q_k_tot(C2==1)=0;
g.Q_k_tot(F==1)=0;


Xnplus1 =  ((g.H*p.dx*p.dy*p.rhos.*X_old + g.qx_tot*p.dt)./(g.H*p.dx*p.dy*p.rhos + g.Q_k_tot*p.rhos*p.dt + 1e-20));



% Only needed during spinup where some divide by zero situations arise



% Apply boundary conditions.
for i = 1:p.nX
    Xtemp = Xnplus1(:,:,i);
    Xtemp2 = g.X(:,:,i);
    Xtemp3 = g.Xr(:,:,i);
    Xtemp(C==1) = Xtemp3(C==1);
    Xtemp(C2==1) = Xtemp2(C2==1);
    % enforce boundary condition: no change
    Xtemp(M==0) = Xtemp2(M==0);
%     Xtemp(g.fl==1) = Xtemp2(g.fl==1);
    Xtemp(p.F==1) = Xtemp2(p.F==1);
    % in X at boundaries.
    Xnplus1(:,:,i) = Xtemp;
end

% %% Finish

H_before=g.H;
H = Hnplus1;
g.H = H;
g.H(p.F==1)=H_before(p.F==1);
g.H(M==0)=H_before(M==0);

% Xnplus1(Xnplus1<=0)=0; % no negative concentrations; not sure how it would arise, though
g.X = Xnplus1;


% % Make the concentrations sum to one as they don't 'know' how much of the
% % other mineral(s) where transported -- 3rd dim is number of mins; This


% X_ratio = sum(g.X,3);
% g.X = g.X./X_ratio;

for i = 1:p.nX
    Xtemp = g.X(:,:,i);
    Xtemp(Xtemp<=0) = p.xr(i);
    g.X(:,:,i) = Xtemp;
end

for i = 1:p.nX
    Xtemp = g.X(:,:,i);
    Xtemp(isnan(Xtemp(:))) = p.xr(i);
    g.X(:,:,i) = Xtemp;
end

X_ratio = sum(g.X,3);
g.X = g.X./X_ratio;

for i = 1:p.nX
    Xtemp = g.X(:,:,i);
    Xtemp(g.H<=0.002) = p.xr(i);
    g.X(:,:,i) = Xtemp;
end




% Ratios must be correct before doing Ns transport

% if anything became a 'channel' stripped to bedrock by weathering and soil
% transport

if sum(isnan(g.X(:)))>0
    fprintf('Erode.m: NaN(s) detected in X before Ns transport\n')
end

Ns_old = g.Ns;

% threshold = 1e10;
% [maxNs, maxNs_ind] = max(g.Ns(:));
% if maxNs > threshold
%     disp('g.Ns has large value (>1e10) in Erode.m before Ns calc')
%     disp(maxNs)
%     [a,b] = ind2sub(size(g.Ns), maxNs_ind);
%     disp([a,b])
%     disp(g.H(a,b))
% end


%g.Ns = ((g.Ns.*H_before*p.dx*p.dy*p.rhos.*X_old(:,:,p.host_min) + g.qns_tot*p.dt)./(H_before*p.dx*p.dy*p.rhos.*X_old(:,:,p.host_min) + g.Q_k_tot*p.rhos.*g.X(:,:,p.host_min)*p.dt));

g.Ns = ((g.Ns.*H_before*p.dx*p.dy*p.rhos.*X_old(:,:,p.host_min) + g.qns_tot*p.dt)./(H_before*p.dx*p.dy*p.rhos.*X_old(:,:,p.host_min) + g.qx_tot(:,:,p.host_min)*p.dt));

g.Ns(isnan(g.Ns(:))) = g.Nzb(isnan(g.Ns(:)));
g.Ns(isinf(g.Ns(:))) = g.Nzb(isinf(g.Ns(:)));

g.Ns(C==1)=g.Nzb(C==1);

g.Ns(M==0)=0; 
g.Ns(p.F==1)=0;
g.Ns(g.Ns<0)=g.Nzb(g.Ns<0);
g.Ns(g.H<=0.002)=g.Nzb(g.H<=0.002);
g.Ns(g.C2==1)=0;





if sum(isnan(g.Ns(:)))>0
    fprintf('Erode.m: NaN(s) detected in Ns after hillslope transport\n')
end



% % Set Ns to bedrock if channel forms so it can recover if it changes back
% to a soil-mantled zone; note that FT.m changed channelized areas back to
% hillslopes if the soil thickness is above 0.005 m, so this doesn't affect
% ongoing incision into soil

% 

g.Ns(C==1)=g.Nzb(C==1);

g.Ns(M==0)=0; 
g.Ns(p.F==1)=0;
g.Ns(g.Ns<0)=g.Nzb(g.Ns<0);
g.Ns(g.H<=0.002)=g.Nzb(g.H<=0.002);
g.Ns(g.C2==1)=0;





