gn=importdata('data_to_anomalies (1).txt');
fi=gn.data(:,1);
la=gn.data(:,2);
H=gn.data(:,3);
g=gn.data(:,4);

sigma=2.25

for i=1:length(fi)
    gamma(i,1)=978032.66*(1+0.0053024*(sin(fi(i)*pi/180)^2-0.00000585*(sin(2*fi(i)*pi/180)^2)));
end

for i=1:length(H)
    Rfa(i,1)=0.30855*H(i)
    RBg(i,1)=-0.04192*sigma*H(i)
    RPP(i,1)=(0.3086-0.08384*sigma)*H(i)
    
end

for i=1:length(H)
    Afa(i,1)=g(i)+Rfa(i)-gamma(i) 
    ABg(i,1)=g(i)+RBg(i)+Rfa(i)-gamma(i) 
    APP(i,1)=g(i)+RPP(i)-gamma(i) 

end

figure('Name','gravity anomalies vs height')
subplot(2,2,1)
plot(fi,gamma,'x')
title('gamma = f(fi)')
xlabel(['\phi','[deg'])
ylabel('gravity [mGals]')


subplot(2,2,2)
p1=polyfit(H,Afa,1);
y1=polyval(p1,H);
plot(H,Afa,'o')
hold on
plot(H,y1)
hold off
title([' Afa = f(H) ',num2str(p1(1,1),'%.4f'),' mGal/m'])
xlabel('H')
ylabel('Afa')

subplot(2,2,3)
p2=polyfit(H,ABg,1);
y2=polyval(p2,H);
plot(H,ABg,'o')
hold on
plot(H,y2)
hold off
title([' ABg = f(H) ',num2str(p2(1,1),'%.4f'),' mGal/m'])
xlabel('H')
ylabel('ABg')

subplot(2,2,4)
p3=polyfit(H,APP,1);
y3=polyval(p3,H);
plot(H,APP,'o')
hold on
plot(H,y2)
hold off
title([' APP = f(H) ',num2str(p3(1,1),'%.4f'),' mGal/m'])
xlabel('H')
ylabel('APP')

figure('Name','gravity anomalies vs lambda')
subplot(1,1,1)
plot(la,gamma,'x')
title('gamma = f(fi)')
xlabel(['\lambda','[deg]'])
ylabel('gravity [mGals]')


%maps - countur plot and 3D plot
x = min(la(:,1)):0.1:max(la(:,1)); 
y = min(fi(:,1)):0.1:max(fi(:,1)); 

[XI,YI] = meshgrid(x,y);

Zfa =  griddata(la(:,1),fi(:,1),Afa(:,1),XI,YI,'cubic');
Zb =  griddata(la(:,1),fi(:,1),ABg(:,1),XI,YI,'cubic');
ZPP =  griddata(la(:,1),fi(:,1),APP(:,1),XI,YI,'cubic');

figure('Name','Free-air anomalies')
subplot(1,2,1)
contourf(XI,YI,Zfa,20), colorbar, hold on
xlabel(['\lambda', '[deg]'])
ylabel(['\phi', '[deg]'])

subplot(1,2,2)
surfc(XI,YI,Zfa)
xlabel(['\lambda', '[deg]'])
ylabel(['\phi', '[deg]'])

figure('Name','Bouguer anomalies')
subplot(1,2,1)
contourf(XI,YI,Zb,20), colorbar, hold on
xlabel(['\lambda', '[deg]'])
ylabel(['\phi', '[deg]'])

subplot(1,2,2)
surfc(XI,YI,Zb)
xlabel(['\lambda', '[deg]'])
ylabel(['\phi', '[deg]'])

figure('Name','Poincare-prey anomalies')
subplot(1,2,1)
contourf(XI,YI,ZPP,20), colorbar, hold on
xlabel(['\lambda', '[deg]'])
ylabel(['\phi', '[deg]'])

subplot(1,2,2)
surfc(XI,YI,ZPP)
xlabel(['\lambda', '[deg]'])
ylabel(['\phi', '[deg]'])


fii=52.33720
lai=19.20506
Hi=122.9
gi=981233.733

Afa_i1=p1(1,1)*Hi+p1(1,2)
Ab_i1=p2(1,1)*Hi+p2(1,2)
App_i1=p3(1,1)*Hi+p3(1,2)

Afa_i2=griddata(la(:,1),fi(:,1),Afa(:,1),lai,fii,'linear')
Ab_i2=griddata(la(:,1),fi(:,1),ABg(:,1),lai,fii,'linear')
App_i2=griddata(la(:,1),fi(:,1),APP(:,1),lai,fii,'linear')


% Calculate Test point data

gn_test=importdata('test points (1).txt');
fi_test=gn_test.data([9,59,109,189,239],2);
la_test=gn_test.data([9,59,109,189,239],3);
H_test=gn_test.data([9,59,109,189,239],4);
g_test=gn_test.data([9,59,109,189,239],5);

sigma=2.25

for i=1:length(fi_test)
    gamma_test(i,1)=978032.66*(1+0.0053024*sin(fi_test(i)*pi/180)^2-0.00000585*(sin(2*fi_test(i)*pi/180)^2));
end

for i=1:length(H_test)
    Rfa_test(i,1)=0.30855*H_test(i)
    RBg_test(i,1)=-0.04192*sigma*H_test(i)
    RPP_test(i,1)=(0.3086-0.08384*sigma)*H_test(i)
    
end

for i=1:length(H_test)
    Afa_test(i,1)=g_test(i)+Rfa_test(i)-gamma_test(i) 
    ABg_test(i,1)=g_test(i)+Rfa_test(i)+RBg_test(i)-gamma_test(i) 
    APP_test(i,1)=g_test(i)+RPP_test(i)-gamma_test(i) 

end

p1_test=polyfit(H_test,Afa_test,1);
p2_test=polyfit(H_test,ABg_test,1);
p3_test=polyfit(H_test,APP_test,1);

fii=52.33720
lai=19.20506
Hi=122.9
gi=981233.733

Afa_i1_RAW=p1_test(1,1)*Hi+p1_test(1,2)
Ab_i1_RAW=p2_test(1,1)*Hi+p2_test(1,2)
App_i1_RAW=p3_test(1,1)*Hi+p3_test(1,2)

Afa_i2_RAW=griddata(la_test(:,1),fi_test(:,1),Afa_test(:,1),lai,fii,'linear')
Ab_i2_RAW=griddata(la_test(:,1),fi_test(:,1),ABg_test(:,1),lai,fii,'linear')
App_i2_RAW=griddata(la_test(:,1),fi_test(:,1),APP_test(:,1),lai,fii,'linear')


% The interpolation based on H

res_fa_1=Afa_i1_RAW-Afa_i1
res_b_1=Ab_i1_RAW-Ab_i1
res_pp_1=App_i1_RAW-App_i1

%The interpolation based on fi, lambda
res_fa_2=Afa_i2_RAW-Afa_i2
res_b_2=Ab_i2_RAW-Ab_i2
res_pp_2=App_i2_RAW-App_i2;

Afai=gi+0.30855*Hi-(978032.66*(1+0.0053024*((sin(fii*pi/180))^2-0.00000585*((sin(2*fii*pi/180)))^2)));

figure('Name','error values (RMSE)')
subplot(2,2,1)
bar([res_fa_1,res_fa_2])
xlabel(['residual free anomalies based on H, (fi, lambda)'])
ylabel(['RMSE'])

subplot(2,2,2)
bar([res_b_1,res_b_2])
xlabel(['residual bouguer anomalies based on H, (fi, lambda)'])
ylabel(['RMSE'])

subplot(2,2,3)
bar([res_pp_1,res_pp_2])
xlabel(['residual Poincare-prey anomalies based on H, (fi, lambda)'])
ylabel(['RMSE'])






