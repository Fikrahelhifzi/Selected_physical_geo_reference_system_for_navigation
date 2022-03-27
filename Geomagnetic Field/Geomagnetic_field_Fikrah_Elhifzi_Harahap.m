%%
clear

%loading data from igrf13coeffs_data.txt 
igrf=importdata('first_dipole_igrf13coeffs.txt');
year=igrf(1,:);
g1_0=igrf(2,:);
g1_1=igrf(3,:);
h1_1=igrf(4,:);

% I have taken the coordiants of Jakarta Special Capital Region, Indonesia
lat=-6.1753942
co_lat=90-lat
lon=106.827183

%compute the magnetic intensity for the point latitude and longitude 

for i=1:length(g1_0)
X(i,1)=-g1_0(i)*sin(co_lat*pi/180)+(((g1_1(i)*cos(lon*pi/180))+(h1_1(i)*sin(lon*pi/180)))*cos(co_lat*pi/180))
Y(i,1)=(g1_1(i)*sin(lon*pi/180))-(h1_1(i)*cos(lon*pi/180))
Z(i,1)=-2*((g1_0(i)*cos(co_lat*pi/180))+(((g1_1(i)*cos(lon*pi/180))+(h1_1(i)*sin(lon*pi/180)))*sin(co_lat*pi/180)))
end

%declination
for i=1:length(year)
D(i,1)=atand(Y(i)/X(i))
end

%inclination
for i=1:length(year)
I(i,1)=atand(Z(i)/sqrt((X(i)^2)+(Y(i)^2)))
end


%visualisation of magnetic intensity and declination and the inclination

figure('Name','magnetic intensity plots')
 subplot(1,3,1)
 plot(year,X);
 title(' X = f(Year)')
 xlabel(' year ')
 ylabel(' X [nanoTesla] ') 
 grid on
  subplot(1,3,2)
 plot(year,Y);
 title(' Y = f(Year) ')
 xlabel(' year ')
 ylabel(' Y [nanoTesla] ') 
 grid on
  subplot(1,3,3)
 plot(year,Z); 
 title(' Z = f(Year)')
 xlabel(' year ')
 ylabel('Z [nanoTesla]') 
 grid on

figure('Name','declination and inclination plots')
subplot(1,2,1)
plot(year,D)
title(' Declination = f(Year) ')
xlabel(['year'])
ylabel('D [°]')
grid on

subplot(1,2,2)
plot(year,I)
title(' Inclination = f(Year) ')
xlabel('year')
ylabel('I [°]')
grid on
  