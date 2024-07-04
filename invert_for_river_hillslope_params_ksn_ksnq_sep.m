clear; clc

%% load the data
load data_Adams2020.mat

%% pull out data
E = data.E_rate./1000;
Eer = data.E_error./1000;

% precip
P = data.mean_trmm;
Per = data.std_trmm;

% ksn
ksn = data.mean_ksn;
ksner = data.std_ksn;

% ksnq
ksnq = data.mean_ksn_q;
ksnqer = data.std_ksn_q;

% slope
slo = data.mean_gradient;
sloer = data.std_gradient;

%% bin data by precip and pack up
p_bins = [1.2,2.0];

% calculate mean precip for each bin
Pm1 = mean(P(P<p_bins(1)));
Pm2 = mean(P(P>=p_bins(1) & P<=p_bins(2)));
Pm3 = mean(P>p_bins(2));

% low precip
E1 = E(P<p_bins(1));
slo1 = slo(P<p_bins(1));
slo1er = sloer(P<p_bins(1));
ksn1 = ksn(P<p_bins(1));
ksn1er = ksner(P<p_bins(1));
ksnq1 = ksnq(P<p_bins(1));
ksnq1er = ksnqer(P<p_bins(1));
dat1 = [E1,slo1,slo1er,ksn1,ksn1er,ksnq1,ksnq1er];
col1 = [0.8 0.1 0.1];

% mid precip
E2 = E(P>=p_bins(1) & P<=p_bins(2));
slo2 = slo(P>=p_bins(1) & P<=p_bins(2));
slo2er = sloer(P>=p_bins(1) & P<=p_bins(2));
ksn2 = ksn(P>=p_bins(1) & P<=p_bins(2));
ksn2er = ksner(P>=p_bins(1) & P<=p_bins(2));
ksnq2 = ksnq(P>=p_bins(1) & P<=p_bins(2));
ksnq2er = ksnqer(P>=p_bins(1) & P<=p_bins(2));
dat2 = [E2,slo2,slo2er,ksn2,ksn2er,ksnq2,ksnq2er];
col2 = [0.95 0.75 0.0];

% high precip
E3 = E(P>p_bins(2));
slo3 = slo(P>p_bins(2));
slo3er = sloer(P>p_bins(2));
ksn3 = ksn(P>p_bins(2));
ksn3er = ksner(P>p_bins(2));
ksnq3 = ksnq(P>p_bins(2));
ksnq3er = ksnqer(P>p_bins(2));
dat3 = [E3,slo3,slo3er,ksn3,ksn3er,ksnq3,ksnq3er];
col3 = [0.1 0.1 0.8];

%% define priors fr parameters interms of mean and std and pack up
K_p = [-25,-5];  %log10 K
%Kp_p = log10(10.^K_p./(mean(P)^(3*0.45))); % log10 Kp
n_p = [1,15];
Sc_p = tand([20,50]); % Sc in degrees
D_p = [-4,-0.1]; % log10 D
priors = [K_p;n_p;Sc_p;D_p];

%% define constants in the model
mn = 0.45; % concavity index used in adams et al. (2020)
pr = 2700; % rock density
ps = 2000; % soil density;
B = pr/ps;



%% run low precip inversion
% pick starting place for the parameters
%param_start = [mean(priors,2);mn;B];

% here near MAP based on ealier model runs
param_start = [-11.2207    3.2775    0.6279   -1.6872 mn B]';

burn_in = 5e3;
n_iter = 5e4;
[params,like_data,mMAP,accrate] = HM_MCMC_river_hillslope_ksn(dat1,param_start,priors,burn_in,n_iter);
save('lowP_params_ksn.mat','params');
save('lowP_like_ksn.mat','like_data');

param_start = [-10.2355    2.8278    0.6086   -1.7405 mn B]';
[params,like_data,mMAPq,accrate] = HM_MCMC_river_hillslope_ksnq(dat1,param_start,priors,burn_in,n_iter);
save('lowP_params_ksnq.mat','params');
save('lowP_like_ksnq.mat','like_data');

%% plot the data
param_names = {'K','n','Sc','D'};
% % plot the chains
% [r,c] = size(params);
% for i = 1:c
%     subplot(c,1,i)
%     if i == 1 || i ==4
%         plot(10.^params(:,i));
%     elseif i == 3
%         plot(atand(params(:,i)));
%     else
%         plot(params(:,i));
%     end
%     xlabel('iterations');
%     ylabel(param_names{i})
% end


f_num = 99;
plot_data(dat1,col1,f_num)
E = logspace(log10(0.001),log10(100),500)./1e3;
paramss = [mMAP,mn,B];
paramssq = [mMAPq,mn,B];
plot_model_ksn_ksnq(paramss,paramssq,E,col1,f_num)
disp('low P results:')
disp(param_names)
disp(mMAP)
disp(mMAPq)
disp(' ')

%% mid precip inversion

param_start = [-14.6166    4.9519    0.5859   -2.1456 mn  B]';

[params,like_data,mMAP,accrate] = HM_MCMC_river_hillslope_ksn(dat2,param_start,priors,burn_in,n_iter);
save('midP_params_ksn.mat','params');
save('midP_like_ksn.mat','like_data');

param_start = [-16.0257    5.3935    0.5894   -2.0059 mn  B]';
[params,like_data,mMAPq,accrate] = HM_MCMC_river_hillslope_ksnq(dat2,param_start,priors,burn_in,n_iter);
save('midP_params_ksnq.mat','params');
save('midP_like_ksnq.mat','like_data');

% plot results
f_num = 99;
plot_data(dat2,col2,f_num)
E = logspace(log10(0.001),log10(100),500)./1e3;
paramss = [mMAP,mn,B];
paramssq = [mMAPq,mn,B];
plot_model_ksn_ksnq(paramss,paramssq,E,col2,f_num)
disp('mid P results:')
disp(param_names)
disp(mMAP)
disp(mMAPq)
disp(' ')


%% high precip inversion
param_start = [-15.5182    5.7144    0.5817   -1.8293 mn  B];
[params,like_data,mMAP,accrate] = HM_MCMC_river_hillslope_ksn(dat3,param_start,priors,burn_in,n_iter);
save('highP_params_ksn.mat','params');
save('highP_like_ksn.mat','like_data');

param_start = [-12.2276    3.7809    0.5776   -1.6665 mn  B];
[params,like_data,mMAPq,accrate] = HM_MCMC_river_hillslope_ksnq(dat3,param_start,priors,burn_in,n_iter);
save('highP_params_ksnq.mat','params');
save('highP_like_ksnq.mat','like_data');

% plot results
f_num = 99;
plot_data(dat3,col3,f_num)
E = logspace(log10(0.001),log10(100),500)./1e3;
paramss = [mMAP,mn,B];
paramssq = [mMAPq,mn,B];
plot_model_ksn_ksnq(paramss,paramssq,E,col3,f_num)
disp('high P results:')
disp(param_names)
disp(mMAP)
disp(mMAPq)
disp(' ')


% set axes
subplot(3,2,1)
xlim([0 5]); ylim([15 40])

subplot(3,2,3)
xlim([0 5]); ylim([0 350])

subplot(3,2,4)
xlim([0 5]); ylim([0 600])

subplot(3,2,5)
xlim([0 500]); ylim([15 45])

subplot(3,2,6)
xlim([0 600]); ylim([15 45])


%% take a look at the acceptance rates and chains
par_name = '_params_ksn';
like_name = '_like_ksn';

pre_names = {'lowP','midP','highP'};
cols =  [0.8 0.1 0.1; 0.95 0.75 0.0; 0.1 0.1 0.8];

for i = 1:length(pre_names)
    
    figure(101)
    subplot(1,2,1)
    load([pre_names{i},like_name,'.mat'])
    plot(like_data(:,3),'-','color',cols(i,:)); hold on

    subplot(1,2,2)
    load([pre_names{i},like_name,'q.mat'])
    plot(like_data(:,3),'-','color',cols(i,:)); hold on

    figure(102)
    load([pre_names{i},par_name,'.mat']);
    n=0;
    for j = 1:4
        subplot(4,2,j+n)
        plot(params(:,j),'-','color', cols(i,:)); hold on
        n = n +1;
    end
    
    load([pre_names{i},par_name,'q.mat']);
    n=1;
    for j = 1:4
        subplot(4,2,j+n)
        plot(params(:,j),'-','color', cols(i,:)); hold on
        n = n+1;
    end
end

figure(101)
subplot(1,2,1)
xlabel('iteration')
ylabel('acceptance rate (%)')
subplot(1,2,2)
xlabel('iteration')

figure(102)
ylabs = {'log10(K)','n','Sc','log10(D)'};
    n=0;
    for j = 1:4
        subplot(4,2,j+n)
        ylabel(ylabs{j})
        n = n +1;
    end
    xlabel('iteration')
subplot(4,2,8)
xlabel('iteration')

%% matrix plots and determine parameters and uncertainties

par_name = '_params_ksn';
like_name = '_like_ksn';

pre_names = {'lowP','midP','highP'};
cols =  [0.8 0.1 0.1; 0.95 0.75 0.0; 0.1 0.1 0.8];

ksn_tab = nan(3,8);
ksnq_tab = nan(3,8);

fign = 200;
for i = 1:length(pre_names)
    
    figure(fign)
    load([pre_names{i},par_name,'.mat']);
    [S,ax]=plotmatrix(params(burn_in+1:end,:));
    ax(1,1).YLabel.String='log10(K)';
    ax(2,1).YLabel.String='n';
    ax(3,1).YLabel.String='Sc';
    ax(4,1).YLabel.String='log10(D)';
    ax(4,1).XLabel.String='log10(K)';
    ax(4,2).XLabel.String='n';
    ax(4,3).XLabel.String='Sc';
    ax(4,4).XLabel.String='log10(D)';
    title([pre_names{i},' k_{sn}'])

    ksn_tab(i,1) = mean(params(burn_in+1:end,1));
    ksn_tab(i,2) = std(params(burn_in+1:end,1));
    ksn_tab(i,3) = mean(params(burn_in+1:end,2));
    ksn_tab(i,4) = std(params(burn_in+1:end,2));
    ksn_tab(i,5) = mean(params(burn_in+1:end,3));
    ksn_tab(i,6) = std(params(burn_in+1:end,3));
    ksn_tab(i,7) = mean(params(burn_in+1:end,4));
    ksn_tab(i,8) = std(params(burn_in+1:end,4));

    fign = fign+1;
    figure(fign)
    load([pre_names{i},par_name,'q.mat']);
    [S,ax]=plotmatrix(params(burn_in+1:end,:));
    ax(1,1).YLabel.String='log10(K)';
    ax(2,1).YLabel.String='n';
    ax(3,1).YLabel.String='Sc';
    ax(4,1).YLabel.String='log10(D)';
    ax(4,1).XLabel.String='log10(K)';
    ax(4,2).XLabel.String='n';
    ax(4,3).XLabel.String='Sc';
    ax(4,4).XLabel.String='log10(D)';
    title([pre_names{i},' k_{snQ}'])

    ksnq_tab(i,1) = mean(params(burn_in+1:end,1));
    ksnq_tab(i,2) = std(params(burn_in+1:end,1));
    ksnq_tab(i,3) = mean(params(burn_in+1:end,2));
    ksnq_tab(i,4) = std(params(burn_in+1:end,2));
    ksnq_tab(i,5) = mean(params(burn_in+1:end,3));
    ksnq_tab(i,6) = std(params(burn_in+1:end,3));
    ksnq_tab(i,7) = mean(params(burn_in+1:end,4));
    ksnq_tab(i,8) = std(params(burn_in+1:end,4));

    fign = fign+1;

end

% save tables with the estimated parameters

att = {'mean_log10(K)','std_log10(K)','mean_n','std_n',...
    'mean_Sc','std_Sc','mean_log10(D)','std_log10(D)'};
T = array2table(ksn_tab,'VariableNames',att);
writetable(T,'ksn_results.xlsx');

att = {'mean_log10(Kp)','std_log10(Kp)','mean_n','std_n',...
    'mean_Sc','std_Sc','mean_log10(D)','std_log10(D)'};
T = array2table(ksnq_tab,'VariableNames',att);
writetable(T,'ksnq_results.xlsx');