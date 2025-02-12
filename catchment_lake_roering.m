output = g_roering_just_lake;
lake_mask = output.output_C2(:,:,52);

clear r

for i = 1:(p.tf/p.saveint/p.dt+1)
    % quartz/host mineral flux to sink(s)
    if i < (p.lake_time + p.dt)/(p.saveint*p.dt) + 1
        qx_lake = output.output_qx_tot_lake(:,:,p.host_min,i)*p.dt;
        
        qx_fluv = output.output_qx_tot_fluv_soil(:,:,p.host_min,i)*p.dt; % + output.output_qx_tot_fluv_bedrock(:,:,p.host_min,i)*p.dt;
        
        
        qns_lake = output.output_qns_tot_lake(:,:,i)*p.dt;
        qns_fluv = output.output_qns_tot_fluv_soil(:,:,i)*p.dt; % + outp
        
        C = output.output_C(:,:,i);
        C2 = output.output_C2(:,:,i);
        
        mass = output.output_qs_tot_lake(:,:,i)*p.dt*p.rhos;
        mass_fluv = output.output_qs_tot_fluv_soil(:,:,i)*p.dt*p.rhos; % + output.output_qs_tot_fluv_bedrock(:,:,i)*p.dt*p.rhor;
        
        total_mass_stream = sum(mass(C==1));
        total_mass_lake = sum(mass(C2==1));
        
        total_mass_fluv = sum(mass_fluv(:));
        
        total_mass = total_mass_stream + total_mass_lake;
        
        r.mass_frac_stream(i) = sum(mass(C==1))/total_mass;
        r.mass_frac_lake(i) = sum(mass(C2==1))/total_mass;
        
        r.mass_fluv_frac_total(i) = total_mass_fluv/total_mass;
        
        
        
        % this is for when zr and qtz are not the 'same' (i.e., two-mineral
        % lithology)
        zr_mass = output.output_qx_tot_lake(:,:,p.zr_min,i)*p.dt;
        zr_mass_fluv = output.output_qx_tot_fluv_soil(:,:,p.zr_min,i)*p.dt;% + output.output_qx_tot_fluv_bedrock(:,:,p.zr_min,i)*p.dt;
        
        plag_mass = output.output_qx_tot_lake(:,:,2,i)*p.dt;
        plag_mass_fluv = output.output_qx_tot_fluv_soil(:,:,2,i)*p.dt;
        
        kspar_mass =  output.output_qx_tot_lake(:,:,3,i)*p.dt;
        kspar_mass_fluv = output.output_qx_tot_fluv_soil(:,:,3,i)*p.dt;
        
        bio_mass =  output.output_qx_tot_lake(:,:,4,i)*p.dt;
        bio_mass_fluv = output.output_qx_tot_fluv_soil(:,:,4,i)*p.dt;
        
        horn_mass =  output.output_qx_tot_lake(:,:,5,i)*p.dt;
        horn_mass_fluv = output.output_qx_tot_fluv_soil(:,:,5,i)*p.dt;
        
        allo_mass =  output.output_qx_tot_lake(:,:,6,i)*p.dt;
        allo_mass_fluv = output.output_qx_tot_fluv_soil(:,:,6,i)*p.dt;
        
        kaol_mass =  output.output_qx_tot_lake(:,:,7,i)*p.dt;
        kaol_mass_fluv = output.output_qx_tot_fluv_soil(:,:,7,i)*p.dt;
        
        % instantaneous denudation rate
        Dact = output.output_Dact(:,:,i);
        Dact(C == 1) = 0;
        
        r.t(i) = output.t(i);
        
        r.lake_be(i) = sum(qns_lake(C==1)) + sum(qns_fluv(:));
        
        % tracks lake area -- may be useful in m^2
        r.lake_area(i) = sum(C(C==1))*p.dx*p.dy;
        
        % total sediment mass in kg
        r.lake_mass(i) = sum(mass(C==1)) + sum(mass_fluv(:));
        
        % total lake quartz/host mineral in kg
        r.lake_qx(i) = sum(qx_lake(C==1)) + sum(qx_fluv(:));
        
        % total Zr flux in kg
        r.lake_zr(i) = sum(zr_mass(C==1)) + sum(zr_mass_fluv(:));
        
        r.lake_qtz(i) = r.lake_qx(i);
        r.lake_plag(i) = sum(plag_mass(C==1)) + sum(plag_mass_fluv(:));
        r.lake_kspar(i) = sum(kspar_mass(C==1)) + sum(kspar_mass_fluv(:));
        r.lake_bio(i) = sum(bio_mass(C==1)) + sum(bio_mass_fluv(:));
        r.lake_horn(i) = sum(horn_mass(C==1)) + sum(horn_mass_fluv(:));
        r.lake_allo(i) = sum(allo_mass(C==1)) + sum(allo_mass_fluv(:));
        r.lake_kaol(i) = sum(kaol_mass(C==1)) + sum(kaol_mass_fluv(:));
        
        r.lake_conc(i) = r.lake_be(i)./r.lake_qx(i);
        
        % Lake Zr conc in kg/kg
        r.lake_zr_conc(i) = (r.lake_zr(i)/r.lake_mass(i));
        
        r.lake_qtz_conc(i) = (r.lake_qtz(i)/r.lake_mass(i));
        r.lake_plag_conc(i) = (r.lake_plag(i)/r.lake_mass(i));
        r.lake_kspar_conc(i) = (r.lake_kspar(i)/r.lake_mass(i));
        r.lake_bio_conc(i) = (r.lake_bio(i)/r.lake_mass(i));
        r.lake_horn_conc(i) = (r.lake_horn(i)/r.lake_mass(i));
        r.lake_allo_conc(i) = (r.lake_allo(i)/r.lake_mass(i));
        r.lake_kaol_conc(i) = (r.lake_kaol(i)/r.lake_mass(i));
        
        % Lake CDF - chemical depletion factor aka W/D
        r.lake_cdf(i) = 1 - p.xr(p.zr_min)/r.lake_zr_conc(i);
        
        % Get production rates from elevations -- ppval outputs to atoms g
        % hence the *1000 for atoms kg
        cosmo_prod_spal = ppval(p.P_spal, (output.output_U(:,:,i) + output.output_H(:,:,i)))*1000;
        cosmo_prod_neg = ppval(p.P_neg, (output.output_U(:,:,i) + output.output_H(:,:,i)))*1000;
        cosmo_prod_fast = ppval(p.P_fast, (output.output_U(:,:,i) + output.output_H(:,:,i)))*1000;
        
        % Get hyposometrically-averaged rates for each production pathway --
        % this could be median just as well (may be better) -- get rid of
        % Channels == 0 to include rates from bedrock channels
        lake_spal = median(cosmo_prod_spal(C==0 & p.F==0));
        lake_neg = median(cosmo_prod_neg(C==0 & p.F==0));
        lake_fast = median(cosmo_prod_fast(C==0 & p.F==0));
        
        lake_zr_hill1 = output.output_Xqtz(:,:,i);
        lake_zr_hill = mean(lake_zr_hill1(C==0 & p.F==0))./p.xr(p.host_min);
        lake_H1 = output.output_H(:,:,i);
        lake_H = mean(lake_H1(C==0 & p.F==0));
       
        lake_zr_rat = r.lake_zr_conc(i)./p.xr(p.zr_min);   

        
        % calculate Dinf w/ CEF for each pathway
        lake_dinf_spal = (lake_spal.*p.L1./r.lake_conc(i)).*(lake_zr_hill.*(1-exp(-p.rhos.*lake_H./p.L1)) + exp(-p.rhos.*lake_H./p.L1));
        lake_dinf_neg = (lake_neg.*p.L2./r.lake_conc(i)).*(lake_zr_hill.*(1-exp(-p.rhos.*lake_H./p.L2)) + exp(-p.rhos.*lake_H./p.L2));
        lake_dinf_fast = (lake_fast.*p.L3./r.lake_conc(i)).*(lake_zr_hill.*(1-exp(-p.rhos.*lake_H./p.L3)) + exp(-p.rhos.*lake_H./p.L3));
        r.lake_dinf(i) = (lake_dinf_spal + lake_dinf_neg + lake_dinf_fast)/p.rhor;
        
%         lake_dinf_spal = (lake_spal.*p.L1./r.lake_conc(i)).*(lake_zr_rat);
%         lake_dinf_neg = (lake_neg.*p.L2./r.lake_conc(i)).*(lake_zr_rat);
%         lake_dinf_fast = (lake_fast.*p.L3./r.lake_conc(i)).*(lake_zr_rat);
%         r.lake_dinf(i) = (lake_dinf_spal + lake_dinf_neg + lake_dinf_fast)/p.rhor;
        
        syms Dinf_lake
        Ps_vpa = lake_spal;
        Pmu_vpa =  lake_neg;
        Pfast_vpa = lake_fast;
        lambda_10Be_vpa = p.lambda_10Be;
        L1_vpa = p.L1;
        L2_vpa = p.L2;
        L4_vpa = p.L3;
        Ns_vpa = r.lake_conc(i);
        % can use sediment Zr_s/Zr_r (lake_zr_conc) or soil Zr_s/Zr_r
        % (lake_zr_hill) -- this is futuristic shit
        Zr_s_vpa = lake_zr_hill*p.xr(p.zr_min);
        Zr_r_vpa = p.xr(p.zr_min);
        H_vpa = lake_H;
        eqn = (Ps_vpa/(lambda_10Be_vpa + Dinf_lake/L1_vpa))*((Zr_s_vpa/Zr_r_vpa)*(1-exp(-p.rhos*H_vpa/L1_vpa)) + exp(-p.rhos*H_vpa/L1_vpa)) + (Pmu_vpa/(lambda_10Be_vpa + Dinf_lake/L2_vpa))*((Zr_s_vpa/Zr_r_vpa)*(1-exp(-p.rhos*H_vpa/L2_vpa)) + exp(-p.rhos*H_vpa/L2_vpa)) + (Pfast_vpa/(lambda_10Be_vpa + Dinf_lake/L4_vpa))*((Zr_s_vpa/Zr_r_vpa)*(1-exp(-p.rhos*H_vpa/L4_vpa)) + exp(-p.rhos*H_vpa/L4_vpa));
        Y_vpa = vpasolve(eqn == Ns_vpa, Dinf_lake);
        r.lake_dinf_vpa(i) = double(Y_vpa(3))/p.rhor;
        
        lake_dinf_spal_marsh = lake_spal*p.L1./r.lake_conc(i);
        lake_dinf_neg_marsh = lake_neg*p.L2./r.lake_conc(i);
        lake_dinf_fast_marsh = lake_fast*p.L3./r.lake_conc(i);
        r.lake_dinf_marsh(i) = (lake_dinf_spal_marsh + lake_dinf_neg_marsh + lake_dinf_fast_marsh)/p.rhor;
        
        % spallation only Dinf
        r.lake_dinf_marsh_spal_only(i) = lake_dinf_spal_marsh/p.rhor;
        
        % Winf from Schatchman et al. (2019)-esque CDF calc
        r.lake_w_marsh(i) = r.lake_dinf_marsh(i)*r.lake_cdf(i);
        
        % Mean actual denudation rate
        r.lake_dact(i) = mean(Dact(C == 0 & p.F == 0));
        H = output.output_H(:,:,i);
        r.mean_H(i) = mean(H(lake_mask==0 & p.F==0 & C==0));
        r.median_H(i) = median(H(lake_mask==0 & p.F==0 & C==0));
        r.mean_H_erosive(i) = r.mean_H(i);
        r.mean_H_old(i) = mean(H(C2==0 & p.F==0 & C==0));
        r.mean_H_num(i) = numel(H(lake_mask==0 & p.F==0 & C==0));
        
        P0 = output.output_P0(:,:,i);
        
        r.mean_P0(i) = mean(P0(lake_mask==0 & p.F==0 & C == 0)); 
        r.min_P0(i) = min(P0(lake_mask==0 & p.F==0 & C == 0));
        r.max_P0(i) = max(P0(lake_mask==0 & p.F==0 & C == 0));
        
    end
    
    if i >= (p.lake_time + p.dt)/(p.saveint*p.dt) + 1
    
        r.t(i) = output.t(i);
        
        qx_lake = output.output_qx_tot_lake(:,:,p.host_min,i)*p.dt;
        
        % can combined this or just do from soil -- most of the time bedrock
        % flux isn't differentiable
        qx_fluv = output.output_qx_tot_fluv_soil(:,:,p.host_min,i)*p.dt; % + output.output_qx_tot_fluv_bedrock(:,:,p.host_min,i)*p.dt;
        
        % total cosmogenic atoms to sink(s)
        qns_lake = output.output_qns_tot_lake(:,:,i)*p.dt;
        qns_fluv = output.output_qns_tot_fluv_soil(:,:,i)*p.dt; % + output.output_qns_tot_fluv_bedrock(:,:,i)*p.dt;
        
        % lake evolves over time so it needs updated each save point
        C = output.output_C(:,:,i);
        C2 = output.output_C2(:,:,i);
        
        % mass fluxes
        mass = output.output_qs_tot_lake(:,:,i)*p.dt*p.rhos;
        mass_fluv = output.output_qs_tot_fluv_soil(:,:,i)*p.dt*p.rhos; % + output.output_qs_tot_fluv_bedrock(:,:,i)*p.dt*p.rhor;
        
        total_mass_stream = sum(mass(C==1));
        total_mass_lake = sum(mass(C2==1));
        
        total_mass_fluv = sum(mass_fluv(:));
        
        total_mass = total_mass_stream + total_mass_lake;
        
        r.mass_frac_stream(i) = sum(mass(C==1))/total_mass;
        r.mass_frac_lake(i) = sum(mass(C2==1))/total_mass;
        
        H = output.output_H(:,:,i);
        
        r.qns_fluv(i) = sum(qns_fluv(H>0.002));
        
        % this is for when zr and qtz are not the 'same' (i.e., two-mineral
        % lithology)
        zr_mass = output.output_qx_tot_lake(:,:,p.zr_min,i)*p.dt;
        zr_mass_fluv = output.output_qx_tot_fluv_soil(:,:,p.zr_min,i)*p.dt; % + output.output_qx_tot_fluv_bedrock(:,:,p.zr_min,i)*p.dt;
        
        plag_mass = output.output_qx_tot_lake(:,:,2,i)*p.dt;
        plag_mass_fluv = output.output_qx_tot_fluv_soil(:,:,2,i)*p.dt;
        
        kspar_mass =  output.output_qx_tot_lake(:,:,3,i)*p.dt;
        kspar_mass_fluv = output.output_qx_tot_fluv_soil(:,:,3,i)*p.dt;
        
        bio_mass =  output.output_qx_tot_lake(:,:,4,i)*p.dt;
        bio_mass_fluv = output.output_qx_tot_fluv_soil(:,:,4,i)*p.dt;
        
        horn_mass =  output.output_qx_tot_lake(:,:,5,i)*p.dt;
        horn_mass_fluv = output.output_qx_tot_fluv_soil(:,:,5,i)*p.dt;
        
        allo_mass =  output.output_qx_tot_lake(:,:,6,i)*p.dt;
        allo_mass_fluv = output.output_qx_tot_fluv_soil(:,:,6,i)*p.dt;
        
        kaol_mass =  output.output_qx_tot_lake(:,:,7,i)*p.dt;
        kaol_mass_fluv = output.output_qx_tot_fluv_soil(:,:,7,i)*p.dt;
        
        
        % instantaneous denudation rate
        Dact = output.output_Dact(:,:,i);
        Dact(C2 == 1) = 0;
        
        % Assumes lake g.C2==1 is ultimate sink for everything -- atoms
        r.lake_be(i) = sum(qns_lake(C2==1)) + sum(qns_lake(C==1)) + sum(qns_fluv(:));
        
        % tracks lake area -- may be useful in m^2
        r.lake_area(i) = sum(C2(C2==1))*p.dx*p.dy;
        
        % total sediment mass in kg
        r.lake_mass(i) = sum(mass(C2==1)) + sum(mass(C==1)) + sum(mass_fluv(:));
        
        % total lake quartz/host mineral in kg
        r.lake_qx(i) = sum(qx_lake(C2==1)) + sum(qx_lake(C==1)) + sum(qx_fluv(:));
        
        % total Zr flux in kg
        r.lake_zr(i) = sum(zr_mass(C2==1)) + sum(zr_mass(C==1)) + sum(zr_mass_fluv(:));
        
        r.lake_qtz(i) = r.lake_qx(i);
        r.lake_plag(i) = sum(plag_mass(C2==1)) + sum(plag_mass(C==1)) + sum(plag_mass_fluv(:));
        r.lake_kspar(i) = sum(kspar_mass(C2==1)) + sum(kspar_mass(C==1)) + sum(kspar_mass_fluv(:));
        r.lake_bio(i) = sum(bio_mass(C2==1)) + sum(bio_mass(C==1)) + sum(bio_mass_fluv(:));
        r.lake_horn(i) = sum(horn_mass(C2==1)) + sum(horn_mass(C==1)) + sum(horn_mass_fluv(:));
        r.lake_allo(i) = sum(allo_mass(C2==1)) + sum(allo_mass(C==1)) + sum(allo_mass_fluv(:));
        r.lake_kaol(i) = sum(kaol_mass(C2==1)) + sum(kaol_mass(C==1)) + sum(kaol_mass_fluv(:));
        
        % Lake cosmo conc in atoms/kg -- usually in atoms/g but kg is easier to
        % work with here
        %r.lake_conc(i) = sum(qns_lake(C2==1))./sum(qx_lake(C2==1));
        r.lake_conc(i) = r.lake_be(i)./r.lake_qx(i);
        
        
        
        % Lake Zr conc in kg/kg
        r.lake_zr_conc(i) = (r.lake_zr(i)/r.lake_mass(i));
        r.lake_qtz_conc(i) = (r.lake_qtz(i)/r.lake_mass(i));
        r.lake_plag_conc(i) = (r.lake_plag(i)/r.lake_mass(i));
        r.lake_kspar_conc(i) = (r.lake_kspar(i)/r.lake_mass(i));
        r.lake_bio_conc(i) = (r.lake_bio(i)/r.lake_mass(i));
        r.lake_horn_conc(i) = (r.lake_horn(i)/r.lake_mass(i));
        r.lake_allo_conc(i) = (r.lake_allo(i)/r.lake_mass(i));
        r.lake_kaol_conc(i) = (r.lake_kaol(i)/r.lake_mass(i));
        
        lake_min_ratio(i) = r.lake_qtz_conc(i) + r.lake_plag_conc(i) + r.lake_kspar_conc(i) + r.lake_bio_conc(i) ...
        + r.lake_horn_conc(i) + r.lake_allo_conc(i) + r.lake_kaol_conc(i);
        
        
    
        % Lake CDF - chemical depletion factor aka W/D
        r.lake_cdf(i) = 1 - p.xr(p.zr_min)/r.lake_zr_conc(i);
        
        r.mean_H(i) = mean(H(C2==0 & p.F==0 & C == 0));
        r.median_H(i) = median(H(C2==0 & p.F==0 & C == 0));
        r.mean_H_old(i) = mean(H(C2==0 & p.F==0 & C == 0));
        r.mean_H_num(i) = numel(H(C2==0 & p.F==0 & C == 0));
        
        P0 = output.output_P0(:,:,i);
        
        r.mean_P0(i) = mean(P0(C2==0 & p.F==0 & C == 0)); 
        r.min_P0(i) = min(P0(C2==0 & p.F==0 & C == 0));
        r.max_P0(i) = max(P0(C2==0 & p.F==0 & C == 0));
        
        % Get production rates from elevations -- ppval outputs to atoms g
        % hence the *1000 for atoms kg
        cosmo_prod_spal = ppval(p.P_spal, (output.output_U(:,:,i) + output.output_H(:,:,i)))*1000;
        cosmo_prod_neg = ppval(p.P_neg, (output.output_U(:,:,i) + output.output_H(:,:,i)))*1000;
        cosmo_prod_fast = ppval(p.P_fast, (output.output_U(:,:,i) + output.output_H(:,:,i)))*1000;
        
        % Get hyposometrically-averaged rates for each production pathway --
        % this could be median just as well (may be better) -- get rid of
        % Channels == 0 to include rates from bedrock channels
        lake_spal = median(cosmo_prod_spal(C2==0 & p.F==0 & C == 0));
        lake_neg = median(cosmo_prod_neg(C2==0 & p.F==0 & C == 0));
        lake_fast = median(cosmo_prod_fast(C2==0 & p.F==0 & C==0));
        
        % all local values for calculation of CEF in inferred denudation rate -- think of it as a worker sampling zr enrichment
        % in the soil and the thickness of soils throughout the basin
        lake_zr_hill1 = output.output_Xqtz(:,:,i);
        lake_zr_hill = mean(lake_zr_hill1(C2==0 & p.F==0 & C==0))./p.xr(p.host_min);
        lake_H1 = output.output_H(:,:,i);
        lake_H = mean(lake_H1(C2==0 & p.F==0 & C==0));
        
        % calculate Dinf w/ CEF for each pathway
        lake_dinf_spal = (lake_spal.*p.L1./r.lake_conc(i)).*(lake_zr_hill.*(1-exp(-p.rhos.*lake_H./p.L1)) + exp(-p.rhos.*lake_H./p.L1));
        lake_dinf_neg = (lake_neg.*p.L2./r.lake_conc(i)).*(lake_zr_hill.*(1-exp(-p.rhos.*lake_H./p.L2)) + exp(-p.rhos.*lake_H./p.L2));
        lake_dinf_fast = (lake_fast.*p.L3./r.lake_conc(i)).*(lake_zr_hill.*(1-exp(-p.rhos.*lake_H./p.L3)) + exp(-p.rhos.*lake_H./p.L3));
        r.lake_dinf(i) = (lake_dinf_spal + lake_dinf_neg + lake_dinf_fast)/p.rhor;
        
        % calculate vpasolve'd Dinf with CEF -- must have symbolic toolkit for
        % MATLAB -- probably could be similarly done w/ OCTAVE's symbolic
        % package -- comment out for speed
        
        syms Dinf_lake
        Ps_vpa = lake_spal;
        Pmu_vpa =  lake_neg;
        Pfast_vpa = lake_fast;
        lambda_10Be_vpa = p.lambda_10Be;
        L1_vpa = p.L1;
        L2_vpa = p.L2;
        L4_vpa = p.L3;
        Ns_vpa = r.lake_conc(i);
        % can use sediment Zr_s/Zr_r (lake_zr_conc) or soil Zr_s/Zr_r
        % (lake_zr_hill) -- this is futuristic shit
        Zr_s_vpa = lake_zr_hill*p.xr(p.host_min);
        Zr_r_vpa = p.xr(p.host_min);
        H_vpa = lake_H;
        eqn = (Ps_vpa/(lambda_10Be_vpa + Dinf_lake/L1_vpa))*((Zr_s_vpa/Zr_r_vpa)*(1-exp(-p.rhos*H_vpa/L1_vpa)) + exp(-p.rhos*H_vpa/L1_vpa)) + (Pmu_vpa/(lambda_10Be_vpa + Dinf_lake/L2_vpa))*((Zr_s_vpa/Zr_r_vpa)*(1-exp(-p.rhos*H_vpa/L2_vpa)) + exp(-p.rhos*H_vpa/L2_vpa)) + (Pfast_vpa/(lambda_10Be_vpa + Dinf_lake/L4_vpa))*((Zr_s_vpa/Zr_r_vpa)*(1-exp(-p.rhos*H_vpa/L4_vpa)) + exp(-p.rhos*H_vpa/L4_vpa));
        Y_vpa = vpasolve(eqn == Ns_vpa, Dinf_lake);
        r.lake_dinf_vpa(i) = double(Y_vpa(3))/p.rhor;
        
        % calculate Dinf w/o CEF
        lake_dinf_spal_marsh = lake_spal*p.L1./r.lake_conc(i);
        lake_dinf_neg_marsh = lake_neg*p.L2./r.lake_conc(i);
        lake_dinf_fast_marsh = lake_fast*p.L3./r.lake_conc(i);
        r.lake_dinf_marsh(i) = (lake_dinf_spal_marsh + lake_dinf_neg_marsh + lake_dinf_fast_marsh)/p.rhor;
        
        % spallation only Dinf
        r.lake_dinf_marsh_spal_only(i) = lake_dinf_spal_marsh/p.rhor;
        
        % Winf from Schatchman et al. (2019)-esque CDF calc
        r.lake_w_marsh(i) = r.lake_dinf_marsh(i)*r.lake_cdf(i);
        
        % Mean actual denudation rate
        r.lake_dact(i) = mean(Dact(C2 == 0 & p.F == 0 & C==0));
        
        % Mean actual denudation rates in net erosive zones (some go
        % depositional along lake edge)
        lake_dact_e = Dact(C2==0 & p.F==0 & C==0);
        r.lake_dact_erosive(i) = mean(lake_dact_e(lake_dact_e>0));
        r.mean_H_erosive(i) = mean(H(lake_dact_e>0));
    
    end
    
end

cmap = mako(7);
figure;
plot(r.t,r.lake_qtz_conc,'color',cmap(1,:));
hold on
plot(r.t,r.lake_plag_conc,'color',cmap(2,:));
plot(r.t,r.lake_kspar_conc,'color',cmap(3,:));
plot(r.t,r.lake_bio_conc,'color',cmap(4,:));
plot(r.t,r.lake_horn_conc,'color',cmap(5,:));
plot(r.t,r.lake_allo_conc,'color',cmap(6,:));
plot(r.t,r.lake_kaol_conc,'color',cmap(7,:));


for i=1:length(r.lake_qtz_conc)
    lake_min_sum(i) = r.lake_qtz_conc(i) + r.lake_plag_conc(i) + r.lake_kspar_conc(i) + r.lake_bio_conc(i) ...
        + r.lake_horn_conc(i) + r.lake_allo_conc(i) + r.lake_kaol_conc(i);
end



