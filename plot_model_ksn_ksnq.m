function plot_model_ksn_ksnq(params,paramsq,E,col,f_num)

% update free parameters
K = 10^params(1);
%Kp = 10^params(2);
n = params(2);
Sc = params(3);
D = 10^params(4);
m = n*params(5);
B = params(6);

% run model with candidate parameters      %plot E-Slope   ---by xcc
[S_mod,ksn_mod,ksnq_mod] = river_nonlinear_hillslope_forward_model(E,K,K,m,n,Sc,D,B);
S_mod = atand(S_mod);
E = E*1000;
figure(f_num)
subplot(3,2,1);
plot(E,S_mod,'k-','color',col); hold on
xlabel('Erosion rate (mm/yr)'); ylabel(['Slope (',char(176),')'])

subplot(3,2,3);                            %plot E-ksn   ---by xcc
plot(E,ksn_mod,'k-','color',col); hold on
xlabel('Erosion rate (mm/yr)'); ylabel('k_{sn}');



subplot(3,2,5);                            %plot slope-ksn   ---by xcc
plot(ksn_mod,S_mod,'k-','color',col); hold on
xlabel('k_{sn}'); ylabel(['Slope (',char(176),')'])



% update free parameters
K = 10^paramsq(1);     %Kp!!!!!!!!!!!!!!!
%Kp = 10^params(2);
n = paramsq(2);
Sc = paramsq(3);
D = 10^paramsq(4);
m = n*paramsq(5);
B = paramsq(6);
% 
E = E./1000;
[S_mod,ksn_mod,ksnq_mod] = river_nonlinear_hillslope_forward_model(E,K,K,m,n,Sc,D,B);
S_mod = atand(S_mod);
E = E*1000;
subplot(3,2,1);                                           %plot E-slope
plot(E,S_mod,'k--','color',col); hold on
%%虚线是用ksnq的数据拟合出来的 E-Slope 的关系 
xlabel('Erosion rate (mm/yr)');ylabel(['Slope (',char(176),')'])

subplot(3,2,4);
plot(E,ksnq_mod,'k-','color',col); hold on                %plot E-ksnQ
xlabel('Erosion rate (mm/yr)'); ylabel('k_{snQ}');

subplot(3,2,6);
plot(ksnq_mod,S_mod,'k-','color',col); hold on            %plot slope-ksnQ
xlabel('k_{snQ}'); ylabel(['Slope (',char(176),')'])
end
