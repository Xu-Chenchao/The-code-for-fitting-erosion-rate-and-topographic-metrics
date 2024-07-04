function [S_mod,ksn_mod,ksnq_mod] = river_nonlinear_hillslope_forward_model(E,K,Kp,m,n,Sc,D,B)

% calculate some parameters
C1 = -Sc^2./(2.*B.*E);

% run the forward model to predict the data
S_mod = nan(length(E),1);
ksn_mod = nan(length(E),1);
ksnq_mod = nan(length(E),1);
for i = 1:length(E)
    Lh = findLh_fzero(E(i),D,Sc,K,n,m,B,1000);
%     if isnan(Lh)
%         Lh = findLh_fzero(E(i),D,Sc,K,n,m,B,1000);
%     end
%     if isnan(Lh)
%         Lh = findLh_fzero(E(i),D,Sc,K,n,m,B,5000);
%     end
    zo = C1(i).*(sqrt(D^2 + ((2.*B.*E(i).*0)./Sc).^2)...
        - D.*log((sqrt(D^2 + ((2.*B.*E(i).*0)./Sc).^2)+D)./((2.*B.*E(i))./Sc)));

    zlh = C1(i).*(sqrt(D^2 + ((2.*B.*E(i).*Lh)./Sc).^2)...
        - D.*log((sqrt(D^2 + ((2.*B.*E(i).*Lh)./Sc).^2)+D)./((2.*B.*E(i))./Sc)));
    S_mod(i) = (zo-zlh)/Lh;
    ksn_mod(i) = (E(i)./K).^(1/n);
    ksnq_mod(i) = (E(i)./(Kp)).^(1/n);
end
end

function [Lh] = findLh_fzero(E,D,Sc,K,n,m,B,Lhi)
% this function will numerically solve for the charateristic hillslope
% length, Lh, by combining the steady-state model for non-linear hillslope
% erosion of Roering et al., (2007) with the stream power incision model
% via Forte et al., (2016) (see supplemental section), which follows Ouimet
% et al., (2009), Howard et al. (1997), and Perron et al. (2008)
% the inputs are erosion rate (E), diffusivity (D), erodability constant 
% (K), slope exponent (n), the area exponent, (m), critical slope angle as 
% rise over run (Sc), and the ratio of soil density to bedrock density (B).

%%
% lines below are based on the assumption that the erosion rate in a
% landscape at all points is a function of both advective processes and
% diffustive process such that E = E(s) + E(h), where E is the total
% erosion, E(s) is advective of stream channel erosion, and E(h) is
% diffusive or hillslope erosion. In this case the channel head is defined
% as the location where the channel gradient and hillslope gradient are
% equal AND E(s) = E(h) such that both are equal to E/2. This methodology
% follows the rational of Howard et al. (1997)

f_all = @(x) ((E./(K.*x.^(2.*m))).^(1/n) - ...
    ((D.*Sc.^2)./(2.*B.*E.*x)).*(sqrt(1+(((2.*B.*E.*x)./(D.*Sc)).^2))-1));

%opts = optimset('display', 'none');
%%
% find minimum misfit between Sr and Sh
%Lh = fzero(@(x) f_all(x),Lhi);
opts = optimset('display', 'none');
Lh = fzero(@(x) f_all(x),Lhi,opts);
end