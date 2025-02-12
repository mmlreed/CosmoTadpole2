function [p g] = Update(p,g)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE elevations using operator splitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sets up lake emplacement within model run for better output
if p.t == p.lake_time
    p.doLake = 1;
    p.lake_level = 260; %default emplacement level -- can be changed
    p.kx_orig = p.kx;
    g.S(g.S <= p.lake_level) = p.lake_level;
    
    g.C2(g.S == p.lake_level) = 1;
    g.Channels(g.S <= p.lake_level) = 0;
    g.U(g.C2 == 1) = p.lake_level;
    g.C2(g.C == 0) = 0;
    
end


g.U(g.U<0)=0;
g.S = g.U + g.H;

% Stream power incision
if p.doStreamPower
    
    g.Nsbeforeerosion_fluv = g.Ns;
    g.Nzbbeforeerosion_fluv = g.Nzb;
    g.Hbeforeerosion_fluv = g.H;
    g.Xbeforeerosion_fluv = g.X;
    g.Sbeforeerosion_fluv = g.S;
    
    [p g] = FT(p,g); % Forward-time explicit
    
    g.Haftererosion_fluv = g.H;
    g.Saftererosion_fluv = g.S;
    % Note that FT performs two additional steps:
    % 1. adjusts time step according to the Courant number
    % 2. records locations of channels
end

% Saving some fluxes for future calculations
g.Ract_fluv_s = (g.Hbeforeerosion_fluv - g.Haftererosion_fluv)/p.dt/(p.rhor/p.rhos);
g.Ract_fluv_b = (g.Sbeforeerosion_fluv - g.Saftererosion_fluv)/p.dt;

g.Ract_fluv(g.Channels==1) = g.Ract_fluv_b(g.Channels==1);
g.Ract_fluv(g.Channels==0) = g.Ract_fluv_s(g.Channels==0);

% Static eroder for unmoving stream network -- may be useful
if p.doStaticChannels
    g = Incise(p,g);
end


% SOURCE TERMS / PERTURBATIONS 
g = Source(p,g);

g.Hbeforesoilprod = g.H;
g.Xbeforesoilprod = g.X;

% Temperature time series

% Constant pre-lake temp
% p.T = 8.5*ones(size(g.H));

% Time-varying temperature
p.T = p.time_temp3(p.tf/1000-p.t/1000).*ones(size(g.H));


% Unused
p.eta_x = 1;
p.P0_x = 1;

% Soil production
% uses Rempel et al. (2016) frost-cracking model to determine
% frost-cracking intensity (FCI). This is then used to calculate a maximum
% % soil production rate

if p.SoilProd_FCI == 1 & p.t >= p.frost_time & p.t <= p.post_lgm_time & p.T < 8.5
    [g,p] = SoilProd_FCI(g,p);
else
%     p.P0=p.P0_orig*ones(size(g.H)).*p.P0_x;
    g = SoilProd(g,p); % normal Heimsath-style soil production 
end
    
g.Haftersoilprod = g.H;

g.S = g.U + g.H;

% some error checking (will halt model)
if sum(isnan(g.S(:)))>0
    fprintf('Update.m: NaN(s) detected in g.S before Cosmo.m\n')
    keyboard;
end

% Cosmogenic production 
g.cosmo_spal = ppval(p.P_spal, g.S)*1000;
g.cosmo_neg = ppval(p.P_neg, g.S)*1000;
g.cosmo_fast = ppval(p.P_fast, g.S)*1000;

% Cosmogenic function
g = Cosmo(p,g);

% Adjust kx (dissolution/alteration constants) based on the ref value at
% 284.15K and Ea for each rxn
% As written the model will adjust these unless a constant temperature is
% used

if p.SoilProd_FCI == 1
    T_kelvin = p.T(32,120) + 273.15; %random point (32,120)
    p.kx = p.kx_orig.*exp(p.Ea_x./p.R_gas.*(1/284.15 - 1/T_kelvin));
    p.kx1 = p.kx1_orig.*exp(p.Ea_x1./p.R_gas.*(1/284.15 - 1/T_kelvin)); % comment out if not using
else
    T_kelvin = p.T(32,120) + 273.15; %random point
    p.kx = p.kx_orig.*exp(p.Ea_x./p.R_gas.*(1/284.15 - 1/T_kelvin));
    p.kx1 = p.kx1_orig.*exp(p.Ea_x1./p.R_gas.*(1/284.15 - 1/T_kelvin)); % comment out if not using
%     p.kx = p.kx_orig;
end


%g = Weather(p,g); % original weathering scheme with no alteration

g = Weather_roering(p,g); % cascading weathering for plag->allophane->kaol

% Calculated chemicla erosion
temp = zeros(size(g.X));
for i = 1:length(p.kx)
    temp(:,:,i) = p.kx(i)*p.Ax(i)*g.X(:,:,i); % - p.sx(i)*p.wx(i)/p.rhos % this is needed if the old Weather() function is used
end
g.Wact = g.H .* sum(temp,3);  % Sum over all mineral phases (m/yr)
clear temp

g.S = g.U + g.H;
g.S(g.U>g.S)=g.U(g.U>g.S);

g.Hbeforeerosion_hill = g.H;
g.Xbeforeerosion = g.X;
g.Nsbeforeerosion = g.Ns;

% general non-linear, thickness-dependent transport
p.D = p.eta.*(1 - exp(-p.beta1.*g.H));

% Fully temperature based implementation with p.Tr obtained from the
% Andersen et al. (2015) code

if p.DasFrostHeave == 1 & p.t >= p.frost_time & p.t <= p.post_lgm_time & p.T < 8.5
    p.D = interp2(p.T0, p.Hs, p.Tr, p.T, g.H); % eta1
    p.D1 = p.eta.*(1 - exp(-p.beta1.*g.H));
    p.D(p.D<p.D1) = p.D1(p.D<p.D1);
    % very small diffusion coefficient so NL8() works
    p.D(g.C2==1)=0.0001;
    p.D(g.Channels==1)=0.0001;
    % if H gets over 6 m -- interp2 outputs NaNs as p.Hs is limited to 6 m
    DH6m = interp2(p.T0, p.Hs, p.Tr, p.T(32,120), 6); % random point to get single number
    p.D(g.H>=6) = DH6m;
else
    p.D = p.eta.*(1 - exp(-p.beta1.*g.H)).*p.eta_x;
end


Sbeforeerosion_hill = g.S;

g = Erode(p,g); % hillslope transport and all that comes with it

Saftererosion_hill = g.S;


g.Haftererosion_hill = g.H;

g.U(g.U<0)=0;
g.S = g.U + g.H;
disp(mean(g.S(p.F==0))) % turn this off it annoys you -- I like constant text output



% Physical erosion and denudation from fluvial incision & hillslope sediment transport (instantaneous)
if p.doStreamPower

    g.Ract_fluv_s = (g.Hbeforeerosion_fluv - g.Haftererosion_fluv)/p.dt/(p.rhor/p.rhos);
    g.Ract_fluv_b = (g.Sbeforeerosion_fluv - g.Saftererosion_fluv)/p.dt;
    
    g.Ract_fluv(g.Channels==1) = g.Ract_fluv_b(g.Channels==1);
    g.Ract_fluv(g.Channels==0) = g.Ract_fluv_s(g.Channels==0);
    
    g.Ract_hill = (g.Hbeforeerosion_hill - g.Haftererosion_hill)/p.dt;  % m/yr
    
    g.qs_tot_fluv_soil = p.dx.*p.dy.*g.Ract_fluv_s;
    g.qs_tot_fluv_soil(g.Channels==1 & g.H == 1e-20) = 0;
    g.qs_tot_fluv_bedrock = p.dx.*p.dy.*g.Ract_fluv_b;
    g.qs_tot_fluv_bedrock(g.Channels==0) = 0;
    
    g.qx_tot_fluv_soil = p.dx.*p.dy.*p.rhos.*g.Ract_fluv_s.*g.Xbeforeerosion_fluv;
    g.qx_tot_fluv_soil(g.Channels==1 & g.H == 1e-20) = 0;
    g.qx_tot_fluv_bedrock = p.dx.*p.dy.*p.rhos.*g.Ract_fluv_b.*g.Xbeforeerosion_fluv;
    g.qx_tot_fluv_bedrock(g.Channels==0) = 0;
    
    g.qns_tot_fluv_soil = p.dx.*p.dy.*p.rhos.*g.Ract_fluv_s.*g.Nsbeforeerosion_fluv;
    g.qns_tot_fluv_soil(g.Channels==1 & g.H == 1e-20) = 0;
    g.qns_tot_fluv_bedrock = p.dx.*p.dy.*p.rhor.*g.Ract_fluv_b.*g.Nzbbeforeerosion_fluv;
    g.qns_tot_fluv_bedrock(g.Channels==0) = 0;
    

    
    g.Dact = (g.Wact./(p.rhor/p.rhos) + g.Ract_hill./(p.rhor/p.rhos) + g.Ract_fluv); % bedrock lowering denudation rates [m yr-1]
    % g.Dact(g.Channels==1) = 0;
    g.Dact(g.C==0) = 0;
    
    g.Ract = g.Ract_hill./(p.rhor/p.rhos) + g.Ract_fluv;
    g.Ract(g.C==0)=0;
end

if p.doStaticChannels
    g.Ract = (g.Hbeforeerosion_hill - g.Haftererosion_hill)/p.dt;  % m/yr
    g.Ract(g.Channels == 1) = p.E;
    g.Dact = g.Ract + g.Wact;
    g.Dact(g.Channels == 1) = 0;
end

g.Wact(g.C==0)=0;




g.CDF = 1 - p.xr(p.zr_min)./g.Xbeforeerosion(:,:,p.zr_min);
g.Winf = g.Dinf.*g.CDF;
g.Einf = g.Dinf - g.Winf;


zr_rat = g.X(:,:,p.zr_min)./p.xr(p.zr_min);
Dinf_spal = (g.cosmo_spal.*p.L1./g.Ns).*(zr_rat.*(1-exp(-p.rhos.*g.H./p.L1)) + exp(-p.rhos.*g.H./p.L1)); 
Dinf_neg = (g.cosmo_neg.*p.L2./g.Ns).*(zr_rat.*(1-exp(-p.rhos.*g.H./p.L2)) + exp(-p.rhos.*g.H./p.L2)); 
Dinf_fast = (g.cosmo_fast.*p.L3./g.Ns).*(zr_rat.*(1-exp(-p.rhos.*g.H./p.L3)) + exp(-p.rhos.*g.H./p.L3));
Dinf_m = Dinf_spal + Dinf_neg + Dinf_fast;
Dinf = Dinf_m./p.rhor;

g.Dinf_alt = Dinf;
g.Dinf_alt(g.C==0)=0;
g.Dinf_alt(g.Channels==1)=0;

end