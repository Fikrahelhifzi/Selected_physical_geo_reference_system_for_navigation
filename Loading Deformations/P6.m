%%  
clear

%load the data from mydata.txt file containing 9th section of data.txt, we need ϕ, λ and M
data=importdata('data.txt');
phi=data(:,2);
la=data(:,3);
M=data(:,4);

%lets suppose ϕs=52 and λs=21
phi_s=52
la_s=21

%compute ∆ϕ and ∆λ for each value of ϕ and λ
for i=1:length(phi)
delta_phi(i,1)=abs(phi_s-phi(i))
delta_la(i,1)=abs(la_s-la(i))
end

%compute ψ "by both methods of calculation" and compute A
for i=1:length(phi)
result_PHI_1(i,1)= acosd(((sin(phi_s*pi/180))*sin(phi(i)*pi/180))+(cos(phi_s*pi/180)*cos(phi(i)*pi/180)*cos(delta_la(i)*pi/180)))
result_PHI_2(i,1)=2*asind(sqrt(((sin(delta_phi(i)*pi/(180*2)))^2)+(cos(phi_s*pi/180)*cos(phi(i)*pi/180)*((sin(delta_la(i)*pi/(180*2))^2)))))
A(i,1)=asind(cos(phi(i)*pi/180)*sin(delta_la(i)*pi/180)/sin(result_PHI_1(i)*pi/180))
end

%Sort the results of result_PHI_1 by ascending
result_PHI_1_sorted = sort(result_PHI_1)

%load the data from grn.txt file containing ψ, Gr and Gh values for the gridded data set
grh=importdata('grn.txt');
val_phi=grh(:,1);
Val_Gr=grh(:,2);
Val_Gh=grh(:,3)

%convert degree to rad and divided by 10^12
for i=1:length(Val_Gr)
Val_Gr_new(i,1)= deg2rad(Val_Gr(i)/10^12);
Val_Gh_new(i,1)=deg2rad(Val_Gh(i)/10^12);
end



%Create a grid to perform interpolation
[~, index] = sort(val_phi);
F_Gr = griddedInterpolant(val_phi(index), Val_Gr_new(index));
F_Gh = griddedInterpolant(val_phi(index), Val_Gh_new(index));

%compute Gr and Gh 
for i=1:length(phi)
Gr(i,1)=F_Gr(result_PHI_1_sorted(i))
Gh(i,1)=F_Gh(result_PHI_1_sorted(i))
end


%compute ∆r
delta_r = 0.0;
for i=1:length(M)
  delta_r = delta_r + Gr(i)*M(i);
end

%compute ∆Hn and  ∆He
delta_Hn = 0.0;
delta_He = 0.0;
for i=1:length(M)
  delta_Hn = delta_Hn + Gh(i)*M(i)*(-cos(A(i)*pi/180));
  delta_He = delta_He + Gh(i)*M(i)*(-sin(A(i)*pi/180));
end


%visualisations of the Green functions 

figure('Name','Plots of Greens functions with respect to ψ')
 subplot(1,2,1)
 plot(val_phi,Val_Gr_new,val_phi,Val_Gh_new);
 title(' G = f(ψ) with linear scale of X-axis')
 xlabel(' ψ[°] ')
 ylabel(' G[·R [m]·ψ[rad]/kg·10^12]') 
 grid on
 legend('Gr','Gh')
 subplot(1,2,2)
 semilogx(val_phi,Val_Gr_new,val_phi,Val_Gh_new); 
 title(' G = f(ψ) with logarithmic scale of X-axis')
 xlabel(' ψ[°] ')
 ylabel(' G[·R [m]·ψ[rad]/kg·10^12]') 
 grid on
 legend('Gr','Gh')


%visualisations of the Computed Green functions 

 figure('Name','Plots of the computed Greens functions with respect to ψ')
 subplot(1,2,1)
 plot(result_PHI_1_sorted,Gr,result_PHI_1_sorted,Gh);
  title(' G = f(ψ) with linear scale of X-axis')
 xlabel(' ψ[°] ')
 ylabel(' G[·R [m]·ψ[rad]/kg·10^12]') 
 grid on
 legend('Gr','Gh')
 subplot(1,2,2)
 semilogx(result_PHI_1_sorted,Gr,result_PHI_1_sorted,Gh); 
 title(' G = f(ψ) with logarithmic scale of X-axis')
 xlabel(' ψ[°] ')
 ylabel(' G[·R [m]·ψ[rad]/kg·10^12]') 
 grid on
 legend('Gr','Gh')