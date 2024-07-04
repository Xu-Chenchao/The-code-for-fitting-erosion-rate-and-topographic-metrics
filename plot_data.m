%%%%%把数据散点投出来
function plot_data(data,col,f_num)
data(:,1) = data(:,1).*1000;
data(:,2) = atand(data(:,2));
figure(f_num)
subplot(3,2,1);
plot(data(:,1),data(:,2),'ko','MarkerSize',7,'markerfacecolor',col); hold on
xlabel('Erosion rate (mm/yr)'); ylabel(['Slope (',char(176),')'])   %char字符代表的ASCII码，char（176）代表°

subplot(3,2,3);
plot(data(:,1),data(:,4),'ko','MarkerSize',7,'markerfacecolor',col); hold on
xlabel('Erosion rate (mm/yr)'); ylabel('k_{sn}'); 

subplot(3,2,4);
plot(data(:,1),data(:,6),'ko','MarkerSize',7,'markerfacecolor',col); hold on
xlabel('Erosion rate (mm/yr)'); ylabel('k_{snQ}'); 

subplot(3,2,5);
plot(data(:,4),data(:,2),'ko','MarkerSize',7,'markerfacecolor',col); hold on
xlabel('k_{sn}'); ylabel(['Slope (',char(176),')'])

subplot(3,2,6);
plot(data(:,6),data(:,2),'ko','MarkerSize',7,'markerfacecolor',col); hold on
xlabel('k_{snQ}'); ylabel(['Slope (',char(176),')'])  
end