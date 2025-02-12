%  Miles Reed (2023)


function g = Cosmo(p,g)

% cosmogenic nuclide production
Ns_spal = ((g.cosmo_spal.*p.L1.*(1-exp(-p.rhos.*g.H./p.L1))));
Ns_neg = ((g.cosmo_neg.*p.L2.*(1-exp(-p.rhos.*g.H./p.L2)))); 
Ns_fast = ((g.cosmo_fast.*p.L3.*(1-exp(-p.rhos.*g.H./p.L3))));

Ns_prod_tot = (Ns_spal + Ns_neg + Ns_fast)./(p.rhos.*g.H); % If using two expoentials for neg muons add here

Ns = g.Ns + Ns_prod_tot*p.dt - ((p.xr(p.host_min).*g.Ns.*p.P0.*exp(-p.alpha.*g.H))./(p.rhos.*g.X(:,:,p.host_min).*g.H))*p.dt + ((p.xr(p.host_min).*g.Nzb.*p.P0.*exp(-p.alpha.*g.H))./(p.rhos.*g.X(:,:,p.host_min).*g.H))*p.dt - p.lambda_10Be*g.Ns*p.dt;

Ns(g.C==0) = g.Nzb(g.C==0);
Ns(g.Channels==1) = g.Nzb(g.Channels==1);
Ns(Ns<0)=0;
Ns(g.H<=0.002)=g.Nzb(g.H<=0.002);
Ns(g.C2==1)=0;

if sum(isnan(Ns(:)))>0
    fprintf('Cosmo.m: NaN(s) detected in Ns after Ns calculation\n')
end

g.Ns = Ns;

% error catching -- can be turned off
threshold = 1e12;
[maxNs,maxNs_ind] = max(g.Ns(:));
if maxNs > threshold
    disp('g.Ns has large value (>1e12) after production')
    [a,b] = ind2sub(size(g.Ns), maxNs_ind);
    disp([a,b])
    disp(g.H(a,b))
    disp(g.X(a,b,p.host_min))
    disp(sum(g.X(a,b,1) + g.X(a,b,2)))
    disp(Ns_prod_tot(a,b))
end

Nzb_prod_spal = g.cosmo_spal.*exp(-p.rhos.*g.H./p.L1)*p.dt;
Nzb_prod_neg =  g.cosmo_neg.*exp(-p.rhos.*g.H./p.L2)*p.dt;
Nzb_prod_fast = g.cosmo_fast.*exp(-p.rhos.*g.H./p.L3)*p.dt;


Nzb_spal_hill = g.Nzb_spal + Nzb_prod_spal - p.P0.*exp(-p.alpha.*g.Hbeforesoilprod).*(g.Nzb_spal./p.L1)*p.dt - g.Nzb_spal*p.lambda_10Be*p.dt;
Nzb_neg_hill = g.Nzb_neg + Nzb_prod_neg - p.P0.*exp(-p.alpha.*g.Hbeforesoilprod).*(g.Nzb_neg./p.L2)*p.dt - g.Nzb_neg*p.lambda_10Be*p.dt;
Nzb_fast_hill = g.Nzb_fast + Nzb_prod_fast - p.P0.*exp(-p.alpha.*g.Hbeforesoilprod).*(g.Nzb_fast./p.L3)*p.dt - g.Nzb_fast*p.lambda_10Be*p.dt;


Nzb_spal_chan = g.Nzb_spal + Nzb_prod_spal - g.Ract_fluv.*p.rhor.*(g.Nzb_spal./p.L1)*p.dt - g.Nzb_spal*p.lambda_10Be*p.dt;
Nzb_neg_chan = g.Nzb_neg + Nzb_prod_neg - g.Ract_fluv.*p.rhor.*(g.Nzb_neg./p.L2)*p.dt - g.Nzb_neg*p.lambda_10Be*p.dt;
Nzb_fast_chan = g.Nzb_fast + Nzb_prod_fast - g.Ract_fluv.*p.rhor.*(g.Nzb_fast./p.L3)*p.dt - g.Nzb_fast*p.lambda_10Be*p.dt;

g.Nzb_spal(g.Channels==0) = Nzb_spal_hill(g.Channels==0);
g.Nzb_neg(g.Channels==0) = Nzb_neg_hill(g.Channels==0);
g.Nzb_fast(g.Channels==0) = Nzb_fast_hill(g.Channels==0);

g.Nzb_spal(g.Channels==1) = Nzb_spal_chan(g.Channels==1);
g.Nzb_neg(g.Channels==1) = Nzb_neg_chan(g.Channels==1);
g.Nzb_fast(g.Channels==1) = Nzb_fast_chan(g.Channels==1);

Nzb = g.Nzb_spal + g.Nzb_neg + g.Nzb_fast;
Nzb(g.C==0) = 0;
Nzb(g.C2==1) = 0;

g.Nzb = Nzb;

zr_rat = g.X(:,:,p.zr_min)./p.xr(p.zr_min);

Dinf_spal = (g.cosmo_spal.*p.L1./g.Ns).*(zr_rat.*(1-exp(-p.rhos.*g.H./p.L1)) + exp(-p.rhos.*g.H./p.L1)); 
Dinf_neg = (g.cosmo_neg.*p.L2./g.Ns).*(zr_rat.*(1-exp(-p.rhos.*g.H./p.L2)) + exp(-p.rhos.*g.H./p.L2)); 
Dinf_fast = (g.cosmo_fast.*p.L3./g.Ns).*(zr_rat.*(1-exp(-p.rhos.*g.H./p.L3)) + exp(-p.rhos.*g.H./p.L3));
Dinf_m = Dinf_spal + Dinf_neg + Dinf_fast;
Dinf = Dinf_m./p.rhor;

Dinf(g.C==0) = 0;
Dinf(g.Channels==1) = 0;
Dinf(g.C2==1) = 0;
g.Dinf = Dinf;

Dinf_spal_lal = (g.cosmo_spal.*p.L1./g.Ns);
Dinf_neg_lal = (g.cosmo_neg.*p.L2./g.Ns);
Dinf_fast_lal = (g.cosmo_fast.*p.L3./g.Ns);
Dinf_m_lal = Dinf_spal_lal + Dinf_neg_lal + Dinf_fast_lal;
Dinf_lal = Dinf_m_lal./p.rhor;

Dinf_lal(g.C==0) = 0;
Dinf_lal(g.Channels==1) = 0;
Dinf_lal(g.C2==1) = 0;
g.Dinf_lal = Dinf_lal;

Dinf_spal_zb = (g.cosmo_spal.*p.L1./g.Nzb).*exp(-p.rhos.*g.H./p.L1);
Dinf_neg_zb = (g.cosmo_neg.*p.L2./g.Nzb).*exp(-p.rhos.*g.H./p.L2);
Dinf_fast_zb = (g.cosmo_fast.*p.L3./g.Nzb).*exp(-p.rhos.*g.H./p.L3);
Dinf_m_zb = Dinf_spal_zb + Dinf_neg_zb + Dinf_fast_zb;
Dinf_zb = Dinf_m_zb/p.rhor;

Dinf_zb(g.C==0) = 0;
g.Dinf_zb = Dinf_zb;

end