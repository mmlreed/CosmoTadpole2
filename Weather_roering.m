function g = Weather(p,g)
%
% Update soil thickness H and mineral abundances X according to the 
% perturbations to these quantities by mineral weathering. This is part of
% a splitting method, whereby this function computes a portion of the 
% temporal changes in H and X, and other functions compute the other 
% portions of those changes.


%% Update mineral abundances X

Xnplus1 = zeros(size(g.X));

% p.kx - dissolution ; p.kx1 - transformation loss (all but 2 is zero: plag->allo and allow->kaol); p.kx2 - transformation gain (all but 2 is zero: allo<-plag and kaol<-allo) -- these should match
% calculate the summation term
BigSigma = 0; 
for i=1:p.nX
    BigSigma = BigSigma + (p.kx(i)*p.Ax(i)*g.X(:,:,i)); % - p.sx(i)*p.wx(i)/p.rhos)  || - gains + losses
end

% 2 - plag; 6 - allophane; 7 - kaolinite -- hardcoded for now -- could revive p.xname

% incongruent weathering loss
Loss = zeros(p.K,p.J,p.nX);
Loss(:,:,2) = (p.kx1(2)*p.Ax(2)*g.X(:,:,2)); 
Loss(:,:,6) = (p.kx1(6)*p.Ax(6)*g.X(:,:,6)); 


% incongruent weathering gain 
Gain = zeros(p.K,p.J,p.nX);
Gain(:,:,6) = Loss(:,:,2);
Gain(:,:,7) = Loss(:,:,6);

% calculate Xn+1 for each mineral species
    
for i=1:p.nX
    Xtemp2 = g.X(:,:,i); 
    Xtemp = g.X(:,:,i);
    DeltaX = p.dt * (-p.kx(i)*p.Ax(i)*g.X(:,:,i) + ...
        g.X(:,:,i).*BigSigma) + Gain(:,:,i) - g.X(:,:,i).*Loss(:,:,i); 
    % DeltaX = p.dt * (-p.kx(i)*p.Ax(i)*g.X(:,:,i) + p.sx(i)*p.wx(i)/p.rhos + ...
    % g.X(:,:,i).*BigSigma);
    Xtemp = Xtemp + DeltaX;
    % Apply boundary conditions.
    Xtemp(g.C==0) = Xtemp2(g.C==0);% enforce boundary condition: no change
    Xtemp(g.Channels==1) = Xtemp2(g.Channels==1);
    Xtemp(p.F==1) = Xtemp2(p.F==1);
    Xtemp(p.catch==1)=p.xr(i);
    Xnplus1(:,:,i) = Xtemp;
end

X_ratio = sum(Xnplus1,3);
Xnplus1 = Xnplus1./X_ratio;



%% Update soil thickness H

% calculate the summation term at time n+1/2
BigSigma = 0; 
for i=1:p.nX
    BigSigma = BigSigma + (p.kx(i)*p.Ax(i)*0.5*(g.X(:,:,i)+Xnplus1(:,:,i)));
    % BigSigma = BigSigma + (p.kx(i)*p.Ax(i)*0.5*(g.X(:,:,i)+Xnplus1(:,:,i)) - ...
       % p.sx(i)*p.wx(i)/p.rhos);
end

Hnplus1 = g.H.*exp(-BigSigma*p.dt);

Hnplus1(Hnplus1<1e-20) = 1e-20;
% Boundary condition: no change in H due to weathering at boundaries.
Hnplus1(g.C==0) = g.H(g.C==0);
Hnplus1(g.Channels==1) = g.H(g.Channels==1); % ALL GOOD

Hnplus1(p.F==1) = g.H(p.F==1);




%% Finish
% assign soil thickness and mineral concentrations at time t + dt to output
% arguments
g.H = Hnplus1;

% g.H(isnan(g.H))=0.002;

if sum(isnan(g.H(:)))>0
    fprintf('Weather.m: NaN(s) detected in H\n')
end


g.X = Xnplus1;

end

