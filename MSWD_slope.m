%% calculate MSWD 用散点约束的函数曲线和数据点本身做计算
function MSWD = MSWD(E_pred,topo_metrics_pred,E,E_error,topo_metrics,topo_metrics_error,data5)

% E_pred and topo_metrics_pred stand for the variables x and y
% E and topo_metrics stand for the observed point. E_error and topo_metrics_error
% are the error of the E and topo_metrics
% data5 is all data (including all bins)

n = length(E); % 残差个数
d = 2; % 自由度设置为2
r2w = []; %用于存储带权重的残差
% emax = mean(data5(:,1)) + 10*std(data5(:,1));%用观测侵蚀速率的平均值加上10倍的标准差作为Sc以上的erosion rate的预测值
% emax = (mean(data5(:,1)) + 4*std(data5(:,1)))/1000;%用观测侵蚀速率的平均值加上10倍的标准差作为Sc以上的erosion rate的预测值 m/yr
emax = (mean(E) + 4*std(E))/1000;%用观测侵蚀速率的平均值加上10倍的标准差作为Sc以上的erosion rate的预测值 m/yr 用分组bin中的侵蚀速率平均值
for k = 1:length(E)
    % 计算拟合曲线与实际数据之间的残差
    if topo_metrics(k) < topo_metrics_pred(end-40)
        residuals = interp1(topo_metrics_pred,E_pred,topo_metrics(k)) - E(k);  %计算MSWD的残差  erosion rate 的MSWD
%         interp_func = scatteredInterpolant(topo_metrics_pred,E_pred');
%         E_value = interp_func(topo_metrics(k));
%         residuals = E_value - E(k);
    elseif topo_metrics(k) > topo_metrics_pred(end-40)
%     TF = isnan(residuals);
%     if TF == 1
%         residuals = emax - E(k);
          residuals = 0;
%     end
    end
   % residuals = interp1(E_pred,topo_metrics_pred,E(k)) - topo_metrics(k);  %计算MSWD的残差  erosion rate 的MSWD
%     slo = topo_metrics_error(k)/E_error(k);
%     x = E_pred;
%     y = E - slo*(x-topo_metrics);  %过观测点的直线的表达式，x范围给的和E_pred一样
%     a = [];
%     for i = 1:length(E)
%         res = (topo_metrics_pred(i) - y(i))^2;
%         a(i) = res;
%     end
%     idx = find(min(a));
%     % 交点坐标
%     x_intersect = E_pred(idx); % 获取交点的 x 坐标
%     y_intersect = topo_metrics_pred(idx); % 获取交点的 y 坐标
%     % 计算曲线在交点处的切线斜率（即交点处的斜率）
%     if idx == 1
%        slope_at_intersection = (topo_metrics_pred(idx+1) - topo_metrics_pred(idx)) / (E_pred(idx+1) - E_pred(idx));
%     elseif idx == length(x_fit)
%        slope_at_intersection = (topo_metrics_pred(idx) - topo_metrics_pred(idx-1)) / (E_pred(idx) - E_pred(idx-1));
%     else
%        slope_at_intersection = (topo_metrics_pred(idx+1) - topo_metrics_pred(idx-1)) / (E_pred(idx+1) - E_pred(idx-1));
%     end
% 
%     % weights = ones(size(residuals)); % 设置权重（这里简单地设置为1）
%     % weights = 1/((slope_at_intersection^2)*(E_error(k)^2)+(1^2)*(topo_metrics(k)^2));
%     weights = 1/((E_error(k)^2)+(slope_at_intersection^2)*(topo_metrics_error(k)^2)); %参考adams的代码写的
%     % 带权重的残差
%     if k == 1
%         weights = 1/(E_error(k)^2 + (E_pred(k+1)-E_pred(k))/(topo_metrics_pred(k+1)-topo_metrics_pred(k))*topo_metrics_error(k)^2);
%     elseif k == length(E)
%         weights = 1/(E_error(k)^2 + (E_pred(k)-E_pred(k-1))/(topo_metrics_pred(k)-topo_metrics_pred(k-1))*topo_metrics_error(k)^2);
%     else
%         weights = 1/(E_error(k)^2 + (E_pred(k+1)-E_pred(k-1))/(topo_metrics_pred(k+1)-topo_metrics_pred(k-1))*topo_metrics_error(k)^2);
%     end
% E unit m/Myr
    if k == 1
        weights = 1/((E_error(k)*1000000)^2 + (E_pred(k+1)-E_pred(k))*1000000/(topo_metrics_pred(k+1)-topo_metrics_pred(k))*topo_metrics_error(k)^2);
    elseif k == length(E)
        weights = 1/((E_error(k)*1000000)^2 + (E_pred(k)-E_pred(k-1))*1000000/(topo_metrics_pred(k)-topo_metrics_pred(k-1))*topo_metrics_error(k)^2);
    else
        weights = 1/((E_error(k)*1000000)^2 + (E_pred(k+1)-E_pred(k-1))*1000000/(topo_metrics_pred(k+1)-topo_metrics_pred(k-1))*topo_metrics_error(k)^2);
    end
    r2w(k) = ((residuals*10^6)^2)*weights;
end
% mean_squared_residuals = sum(r2w) / (n-d); % 计算平均平方残差 %计算除NAN以外的值的和
mean_squared_residuals = sum(r2w) / (nnz(r2w)-d); % 计算平均平方残差 %计算除NAN以外的值的和
% mean_squared_residuals = nansum(r2w) / (nnz(~isnan(r2w))-d); % 计算平均平方残差 %计算除NAN以外的值的和
MSWD = mean_squared_residuals; % MSWD即为平均平方残差
end