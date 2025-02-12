function [p,g] = TadpoleRun(p,g)

% TadpoleRun.m
%
% Performs main iteration loop of Tadpole

n=0;

while p.t < p.tf
    
    n = n + 1;  
    
    %%%%%%%%%%%%%%%%%%%%%% UPDATE ELEVATIONS %%%%%%%%%%%%%%%%%%%%%%%

    [p g] = Update(p,g);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%% INCREMENT TIME %%%%%%%%%%%%%%%%%%%%%%%%%
    
    p.t = p.t + p.dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if it is time to redraw the plot, do so.
    if p.doDrawPlot
        if ~rem(n,p.plotint)
            DrawPlot(n,p,g)
        end        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % All these are needed to compare basin-averaged Dinf to Dact
    
    if p.doSaveOutput
        if ~rem(n,p.saveint)
            p.lastsave = n/p.saveint + 1;
            g.output_U(:,:,p.lastsave) = g.U;
            g.output_H(:,:,p.lastsave) = g.H;
%             g.output_H_beforefluv(:,:,p.lastsave) = g.Hbeforeerosion_fluv;
%             g.output_H_beforehill(:,:,p.lastsave) = g.Hbeforeerosion_hill;
            g.output_Xqtz(:,:,p.lastsave) = g.X(:,:,p.host_min);
%             g.output_Xqtz_beforefluv(:,:,p.lastsave) = g.Xbeforeerosion_fluv(:,:,p.host_min);
%             g.output_Xqtz_beforehill(:,:,p.lastsave) = g.Xbeforeerosion(:,:,p.host_min);
            g.output_Ns(:,:,p.lastsave) = g.Ns;
%             g.output_Ns_beforefluv(:,:,p.lastsave) = g.Nsbeforeerosion_fluv;
%             g.output_Ns_beforehill(:,:,p.lastsave) = g.Nsbeforeerosion;
            g.output_Dact(:,:,p.lastsave) = g.Dact;
            g.output_Dinf(:,:,p.lastsave) = g.Dinf;
            g.output_Ract(:,:,p.lastsave) = g.Ract;
%             g.output_Ract_fluv(:,:,p.lastsave) = g.Ract_fluv;
%             g.output_Ract_hill(:,:,p.lastsave) = g.Ract_hill;
            g.output_C2(:,:,p.lastsave) = g.C2;
            g.output_C(:,:,p.lastsave) = g.Channels;
            g.output_qs_tot_lake(:,:,p.lastsave) = g.Q_k_tot;
            g.output_qx_tot_lake(:,:,:,p.lastsave) = g.qx_tot_lake;
            g.output_qns_tot_lake(:,:,p.lastsave) = g.qns_tot_lake;
            g.output_qns_tot_fluv_soil(:,:,p.lastsave) = g.qns_tot_fluv_soil;
            g.output_qs_tot_fluv_soil(:,:,p.lastsave) = g.qs_tot_fluv_soil;
            g.output_qx_tot_fluv_soil(:,:,:,p.lastsave) = g.qx_tot_fluv_soil;
            g.output_P0(:,:,p.lastsave) = g.P0;
            g.output_CDF(:,:,p.lastsave) = g.CDF;
            
            g.t(p.lastsave) = p.t;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

p.iterations = n;