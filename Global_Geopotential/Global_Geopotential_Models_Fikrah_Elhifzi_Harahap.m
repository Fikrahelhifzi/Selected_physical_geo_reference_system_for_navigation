% Computation gravity anomaly using data anomalies
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
    
end

for i=1:length(H)
    Afa(i,1)=g(i)+Rfa(i)-gamma(i) 
    ABg(i,1)=g(i)+RBg(i)+Rfa(i)-gamma(i) 

end



% S (Grace) for free air gravity anomalies. 'ITSG-Grace2018 S in 2019 with the degree 200'
only_s_grace_free_air=fopen('Only_S_TSG-Grace2018s_Free Air Gravity Anomaly.txt');
for i=1:32
fgetl(only_s_grace_free_air);
end

a_1=textscan(only_s_grace_free_air, '%s %f %f %f %f');
long_s_grace_free_air=a_1{2};
lat_s_grace_free_air=a_1{3};
h_s_grace_free_air=a_1{4};
afa_s_grace_free_air=a_1{5};
fclose(only_s_grace_free_air);


% S (Grace) for Bouguer gravity anomalies. 'ITSG-Grace2018 S in 2019 with the degree 200'
only_s_grace_Bouguer=fopen('Only_S_ITSG-Grace2018s_Bougar.txt');
for i=1:36
fgetl(only_s_grace_Bouguer);
end

a_2=textscan(only_s_grace_Bouguer, '%s %f %f %f %f');
long_s_grace_Bouguer=a_2{2};
lat_s_grace_Bouguer=a_2{3};
h_s_grace_Bouguer=a_2{4};
afa_s_grace_Bouguer=a_2{5};
fclose(only_s_grace_Bouguer);


% AGS Up to 720 for free air gravity anomalies. 'EIGEN-6C in 2011 with the degree 1420'
AGS_up_to_720_free_air=fopen('AGS_Up_To_720_EIGEN-6C_Free Air.txt');
for i=1:32
fgetl(AGS_up_to_720_free_air);
end

d_1=textscan(AGS_up_to_720_free_air, '%s %f %f %f %f');
long_AGS_up_to_720_free_air=d_1{2};
lat_AGS_up_to_720_free_air=d_1{3};
h_AGS_up_to_720_free_air=d_1{4};
afa_AGS_up_to_720_free_air=d_1{5};
fclose(AGS_up_to_720_free_air);


% AGS Up to 720 for Bouguer gravity anomalies. 'EIGEN-6C in 2011 with the degree 1420'
AGS_up_to_720_Bouguer=fopen('AGS_Up_To_720_EIGEN-6C_Bougar.txt');
for i=1:36
fgetl(AGS_up_to_720_Bouguer);
end

d_2=textscan(AGS_up_to_720_Bouguer, '%s %f %f %f %f');
long_AGS_up_to_720_Bouguer=d_2{2};
lat_AGS_up_to_720_Bouguer=d_2{3};
h_AGS_up_to_720_Bouguer=d_2{4};
afa_AGS_up_to_720_Bouguer=d_2{5};
fclose(AGS_up_to_720_Bouguer);


% AGS Up to 2100 for free air gravity anomalies. 'EIGEN-6C4 in 2014 with the degree 2190'
AGS_up_to_2100_free_air=fopen('AGS_Up_To_2100_EIGEN-6C4_Free Air.txt');
for i=1:32
fgetl(AGS_up_to_2100_free_air);
end

e_1=textscan(AGS_up_to_2100_free_air, '%s %f %f %f %f');
long_AGS_up_to_2100_free_air=e_1{2};
lat_AGS_up_to_2100_free_air=e_1{3};
h_AGS_up_to_2100_free_air=e_1{4};
afa_AGS_up_to_2100_free_air=e_1{5};
fclose(AGS_up_to_2100_free_air);



% AGS Up to 2100 for Bouguer gravity anomalies. 'EIGEN-6C4 in 2014 with the degree 2190'
AGS_up_to_2100_Bouguer=fopen('AGS_Up_To_2100_EIGEN-6C4_Bougar.txt');
for i=1:36
fgetl(AGS_up_to_2100_Bouguer);
end

e_2=textscan(AGS_up_to_2100_Bouguer, '%s %f %f %f %f');
long_AGS_up_to_2100_Bouguer=e_2{2};
lat_AGS_up_to_2100_Bouguer=e_2{3};
h_AGS_up_to_2100_Bouguer=e_2{4};
afa_AGS_up_to_2100_Bouguer=e_2{5};
fclose(AGS_up_to_2100_Bouguer);



% Residuum FA1 'only s satellite' for Free air gravity anomalies.'ITSG-Grace2018 S in 2019 with the degree 200'
Residuum_FA1=Afa-afa_s_grace_free_air

% Residuum FA2 AGS Up to 720 for free air gravity anomalies.'EIGEN-6C in 2011 with the degree 1420'
Residuum_FA2=Afa-afa_AGS_up_to_720_free_air

% Residuum FA2 AGS Up to 2100 for free air gravity anomalies. 'EIGEN-6C4 in 2014 with the degree 2190'
Residuum_FA3=Afa-afa_AGS_up_to_2100_free_air


% Residuum BgA1 'only s satellite' for Bouguer gravity anomalies.'ITSG-Grace2018 S in 2019 with the degree 200'
Residuum_BgA1=Afa-afa_s_grace_Bouguer

% Residuum FA2 AGS Up to 720 for Bouguer gravity anomalies.'EIGEN-6C in 2011 with the degree 1420'
Residuum_BgA2=Afa-afa_AGS_up_to_720_Bouguer

% Residuum FA2 AGS Up to 2100 for Bouguer gravity anomalies. 'EIGEN-6C4 in 2014 with the degree 2190'
Residuum_BgA3=Afa-afa_AGS_up_to_2100_Bouguer

% Statistical
RMSE_Residuum_FA1=sqrt((Residuum_FA1.^2)/length(Afa))
RMSE_Residuum_FA2=sqrt((Residuum_FA2.^2)/length(Afa))
RMSE_Residuum_FA3=sqrt((Residuum_FA3.^2)/length(Afa))

RMSE_Residuum_BgA1=sqrt((Residuum_BgA1.^2)/length(Afa))
RMSE_Residuum_BgA2=sqrt((Residuum_BgA2.^2)/length(Afa))
RMSE_Residuum_BgA3=sqrt((Residuum_BgA3.^2)/length(Afa))

Residuum_FA1
Residuum_FA2
Residuum_FA3

Residuum_BgA1
Residuum_BgA2
Residuum_BgA3


mean_Residuum_FA1=mean(Residuum_FA1) 
mean_Residuum_FA2=mean(Residuum_FA2) 
mean_Residuum_FA3=mean(Residuum_FA3) 

mean_Residuum_BgA1=mean(Residuum_BgA1) 
mean_Residuum_BgA2=mean(Residuum_BgA2) 
mean_Residuum_BgA3=mean(Residuum_BgA3) 

max_Residuum_FA1=max(Residuum_FA1) 
max_Residuum_FA2=max(Residuum_FA2) 
max_Residuum_FA3=max(Residuum_FA3) 

max_Residuum_BgA1=max(Residuum_BgA1) 
max_Residuum_BgA2=max(Residuum_BgA2) 
max_Residuum_BgA3=max(Residuum_BgA3) 


min_Residuum_FA1=min(Residuum_FA1) 
min_Residuum_FA2=min(Residuum_FA2) 
min_Residuum_FA3=min(Residuum_FA3) 

min_Residuum_BgA1=min(Residuum_BgA1) 
min_Residuum_BgA2=min(Residuum_BgA2) 
min_Residuum_BgA3=min(Residuum_BgA3) 


mean_RMSE_Residuum_FA1=mean(RMSE_Residuum_FA1)
mean_RMSE_Residuum_FA2=mean(RMSE_Residuum_FA2)
mean_RMSE_Residuum_FA3=mean(RMSE_Residuum_FA2)

mean_RMSE_Residuum_BgA1=mean(RMSE_Residuum_BgA1)
mean_RMSE_Residuum_BgA2=mean(RMSE_Residuum_BgA2)
mean_RMSE_Residuum_BgA3=mean(RMSE_Residuum_BgA3)

fii=52.33720
lai=19.20506
Hi=122.9
gi=981233.733
Afai=gi+0.30855*Hi-(978032.66*(1+0.0053024*((sin(fii*pi/180))^2-0.00000585*((sin(2*fii*pi/180)))^2)));

% Visualization statistical free air gravity anomalies
figure('Name','Statistical_Residuum free air gravity anomalies')
subplot(2,2,1)
bar([mean_Residuum_FA1,mean_Residuum_FA2,mean_Residuum_FA3])
xlabel(['Residuum free air gravity anomalies in every model'])
ylabel(['Mean'])

subplot(2,2,2)
bar([max_Residuum_FA1,max_Residuum_FA2,max_Residuum_FA3])
xlabel(['Residuum free air gravity anomalies in every model'])
ylabel(['Max'])


subplot(2,2,3)
bar([min_Residuum_FA1,min_Residuum_FA2,min_Residuum_FA3])
xlabel(['Residuum free air gravity anomalies in every model'])
ylabel(['Min'])

subplot(2,2,4)
bar([mean_RMSE_Residuum_FA1,mean_RMSE_Residuum_FA2,mean_RMSE_Residuum_FA3])
xlabel(['Residuum free air gravity anomalies in every model'])
ylabel(['Mean RMSE'])


% Visualization statistical Bouguer gravity anomalies
figure('Name','Statistical_Residuum Bouguer gravity anomalies')
subplot(2,2,1)
bar([mean_Residuum_BgA1,mean_Residuum_BgA2,mean_Residuum_BgA3])
xlabel(['Residuum  Bouguer gravity anomalies in every model'])
ylabel(['Mean'])

subplot(2,2,2)
bar([max_Residuum_BgA1,max_Residuum_BgA2,max_Residuum_BgA3])
xlabel(['Residuum  Bouguer gravity anomalies in every model'])
ylabel(['Max'])

subplot(2,2,3)
bar([min_Residuum_BgA1,min_Residuum_BgA2,min_Residuum_BgA3])
xlabel(['Residuum  Bouguer gravity anomalies in every model'])
ylabel(['Min'])

subplot(2,2,4)
bar([mean_RMSE_Residuum_BgA1,mean_RMSE_Residuum_BgA2,mean_RMSE_Residuum_BgA3])
xlabel(['Residuum free air gravity anomalies in every model'])
ylabel(['Mean RMSE'])















