function [p g] = TadpoleInitialize(initial,p,g)

% TadpoleInitialize.m
%
% Performs initialization steps for Tadpole


%%%%%%%%%%%%%%%%%%
% TIME VARIABLES % 
%%%%%%%%%%%%%%%%%%
                           
p.t = 0; % time in yr
p.dt = p.dtmax;
if p.doAdaptiveTimeStep && ~isfield(p,'Courant')
    p.Courant = 0.5; % default maximum Courant number
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL AND BOUNDARY CONDITIONS % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[p.K p.J] = size(initial);

g.p_save = p;

% g.H = ones(size(g.U))*0.5;

g.S = g.U + g.H;
% g.Uprev = g.U; % elevations at previous time step

[p g] = BoundaryMat(p,g);

g.H(g.C==0)=1e-20;
% g.S(g.C==0)=g.U(g.C==0);

% if p.doStreamPower %&& ~p.doChannelDiffusion
%     g.Channels = ones(p.K,p.J);
%     [p g] = MarkChannels(p,g);
%     g.H(g.Channels==1)=1e-20;
% else
%     g.Channels = zeros(p.K,p.J);
% end

% g.S(g.Channels==1) = g.U(g.Channels==1); % only w/ threshold


%%%%%%%%%%%%%%%%%
% SET UP OUTPUT % 
%%%%%%%%%%%%%%%%%

if p.doSaveOutput
    
    % assign a default name for the output if it was not specified
    if ~isfield(p,'runname')
        p.runname = 'default_run_name';
    end
   
    p.saveint = ceil(p.saveint);
    N = round(p.tf/p.dt);
    
    g.output_U = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    g.output_H = zeros(p.K,p.J,ceil(N/p.saveint)+1);
%     g.output_H_beforefluv = zeros(p.K,p.J,ceil(N/p.saveint)+1);
%     g.output_H_beforehill = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    g.output_Xqtz = zeros(p.K,p.J,ceil(N/p.saveint)+1);
%     g.output_Xqtz_beforefluv = zeros(p.K,p.J,ceil(N/p.saveint)+1);
%     g.output_Xqtz_beforehill = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    g.output_Ns = zeros(p.K,p.J,ceil(N/p.saveint)+1);
%     g.output_Ns_beforefluv = zeros(p.K,p.J,ceil(N/p.saveint)+1);
%     g.output_Ns_beforehill = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    g.output_Dact = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    g.output_Dinf = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    g.output_Ract = zeros(p.K,p.J,ceil(N/p.saveint)+1);
%     g.output_Ract_fluv = zeros(p.K,p.J,ceil(N/p.saveint)+1);
%     g.output_Ract_hill = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    g.output_C2 = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    g.output_C = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    g.output_qs_tot_lake = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    g.output_qx_tot_lake = zeros(p.K,p.J,p.nX,ceil(N/p.saveint)+1);
    g.output_qns_tot_lake = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    g.output_qs_tot_fluv_soil = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    g.output_qs_tot_fluv_bedrock = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    g.output_qx_tot_fluv_soil = zeros(p.K,p.J,p.nX,ceil(N/p.saveint)+1);
    g.output_qx_tot_fluv_bedrock = zeros(p.K,p.J,p.nX,ceil(N/p.saveint)+1);
    g.output_qns_tot_fluv_soil = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    g.output_qns_tot_fluv_bedrock = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    g.output_P0 = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    g.output_CDF = zeros(p.K,p.J,ceil(N/p.saveint)+1);
    
    g.output_U(:,:,1) = g.U; % the initial surface
    g.output_H(:,:,1) = g.H;
%     g.output_H_beforefluv(:,:,1) = g.Hbeforeerosion_fluv;
%     g.output_H_beforehill(:,:,1) = g.Hbeforeerosion_hill;
    g.output_Xqtz(:,:,1) = g.X(:,:,p.host_min);
%     g.output_Xqtz_beforefluv(:,:,1) = g.Xbeforeerosion_fluv(:,:,p.host_min);
%     g.output_Xqtz_beforefluv(:,:,1) = g.Xbeforeerosion(:,:,p.host_min);
    g.output_Ns(:,:,1) = g.Ns;
%     g.output_Ns_beforefluv(:,:,1) = g.Nsbeforeerosion_fluv;
%     g.output_Ns_beforehill(:,:,1) = g.Nsbeforeerosion;
    g.output_Dact(:,:,1) = g.Dact;
    g.output_Dinf(:,:,1) = g.Dinf;
    g.output_Ract(:,:,1) = g.Ract;
%     g.output_Ract_fluv(:,:,1) = g.Ract_fluv;
%     g.output_Ract_hill(:,:,1) = g.Ract_hill;
    g.output_C2(:,:,1) = g.C2;
    g.output_C(:,:,1) = g.Channels;
    g.output_qs_tot_lake(:,:,1) = g.Q_k_tot;
    g.output_qx_tot_lake(:,:,:,1) = g.qx_tot_lake;
    g.output_qns_tot_lake(:,:,1) = g.qns_tot_lake;
    g.output_qs_tot_fluv_soil(:,:,1) = g.qs_tot_fluv_soil;
    g.output_qs_tot_fluv_bedrock(:,:,1) = g.qs_tot_fluv_bedrock;
    g.output_qx_tot_fluv_soil(:,:,:,1) = g.qx_tot_fluv_soil;
    g.output_qx_tot_fluv_bedrock(:,:,:,1) = g.qx_tot_fluv_bedrock;
    g.output_qns_tot_fluv_soil(:,:,1) = g.qns_tot_fluv_soil;
    g.output_qns_tot_fluv_bedrock = g.qs_tot_fluv_bedrock;
    g.output_P0(:,:,1) = p.P0*ones(size(g.H));
    g.output_CDF(:,:,1) = g.CDF;
    g.t = 0; % vector that will hold the times corresponding to the saved grids
end


%%%%%%%%%%%%%%%%%%%%%
% MEMORY ALLOCATION % 
%%%%%%%%%%%%%%%%%%%%%
                           
g.Ract_fluv = zeros(size(g.U));


%%%%%%%%%%%%%%%%%%%%%
%   SET UP PLOT     % 
%%%%%%%%%%%%%%%%%%%%%

if p.doDrawPlot
    p.plotint = round(p.plotint);
    [p g] = SetUpPlot(p,g);
end
