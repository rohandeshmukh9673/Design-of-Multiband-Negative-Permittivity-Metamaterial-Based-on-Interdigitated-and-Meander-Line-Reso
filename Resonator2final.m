%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retriveal of effective metamaterial parameters from transmission reflection data.
%
% Author: Zsolt Szabó
%
% The details of the algorithm are described in the paper:
%
% Zsolt Szabó, Park Gi-Ho, Ravi Hedge and Er-Ping Li, A unique extraction
% of metamaterial parameters based on Kramers-Kronig relationship, 
% IEEE Transactions on Microwave Theory and Techniques, 
% Vol. 58, Nr. 10, pp. 2646-2653, October 2010.
%
% The input data to the code are the the magnitude and the phase of the S parameters
% and the effective thickness of the metamaterial.
%
% Some modifications are carried out in order to achieve desired result for proposed multiband negative metamaterial resonator 
% such as manual thresholding (th)
% Deshmukh R, Marathe D, Kulat KD. Design of Multiband Negative Permittivity Metamaterial Based on Interdigitated and Meander Line Resonator. 
%In 2019 National Conference on Communications (NCC) 2019 Feb 20 (pp. 1-6). IEEE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
th=0.4e-1;
%load the S parameters
load_choose = 0;
if load_choose == 0 %the S parameters of the homogeneous slab with 40 nm thickness
    path = 'C:\Users\Administrator\Desktop\Rohan\kk data\8mm\';
    %the effective thickness of the metamaterial
    d_eff = 8e-3;  %[m]
else                %the S parameters of the homogeneous slab with 200 nm thickness
    path = 'd:\Home\Cikk\EffParam_KK\Code\200nm\';
    %the effective thickness of the metamaterial
    d_eff = 200e-9; %[m]
end;
    
S11_abs   = load([path 'S11_abs.csv']);
S11_phase = load([path 'S11_phase.csv']);
S21_abs   = load([path 'S21_abs.csv']);
S21_phase = load([path 'S21_phase.csv']);

UnitFrec = 1e9;

%for plotting
f_min   = 2;     
f_max   = 12;
delta_f = 1;

TextUnitFrec = 'GHz';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%permittivity of free space
eps_0 = 8.85418e-12; % [F/m]
%permeability of free space 
mu_0  = 1.2566e-6;   % [H/m]
%the speed of light
speed_c = 1.0/sqrt(eps_0*mu_0); 

n_S11_abs = length(S11_abs);

%search for frequencies larger than zero
ind_frec = 1;
for i = 1:n_S11_abs
    if(S11_abs(i,1) > 1.0e-16)
        ind_frec = i;
        break;
    end;
end;

frec   = S11_abs(ind_frec:n_S11_abs,1);
n_frec = length(frec);

%frequency in rad/s
omega  = 2.0*pi*frec*UnitFrec;
%wavenumber in 1/m
k0     = omega/speed_c;

%calculate the complex S parameters
S11 = S11_abs(ind_frec:n_S11_abs,2).*(cos(S11_phase(ind_frec:n_S11_abs,2)*pi/180) - sqrt(-1)*sin(S11_phase(ind_frec:n_S11_abs,2)*pi/180));
S21 = S21_abs(ind_frec:n_S11_abs,2).*(cos(S21_phase(ind_frec:n_S11_abs,2)*pi/180) - sqrt(-1)*sin(S21_phase(ind_frec:n_S11_abs,2)*pi/180));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Step - 4 calculate the wave impedance
Z_eff = sqrt( ((1.0 + S11).^2 - S21.^2)./((1.0 - S11).^2 - S21.^2) );
exp_ink0d = S21./( 1.0 - S11.*((Z_eff - 1.0)./(Z_eff + 1.0)));

%choose the proper sign for re(Z_eff) and im(n)
for i = 1:n_frec
    re_Z_eff = real(Z_eff(i));
    if (abs(re_Z_eff) > th) %if the magnitude of Z_eff is large enough
        if ( re_Z_eff < 0 )
            Z_eff(i) = -Z_eff(i);
            exp_ink0d(i) = S21(i)./( 1.0 - S11(i).*((Z_eff(i) - 1.0)./(Z_eff(i) + 1.0)));
        end;
    else
        if (abs(exp_ink0d(i)) > 1.0)
            Z_eff(i) = -Z_eff(i);
            exp_ink0d(i) = S21(i)./( 1.0 - S11(i).*((Z_eff(i) - 1.0)./(Z_eff(i) + 1.0)));
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step - 5 calculate the complex refractive index 
logExp_ink0d = log(exp_ink0d);
n_eff = (imag(logExp_ink0d) - sqrt(-1)*real(logExp_ink0d))./(k0*d_eff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step - 6 Apply Kramers-Kronig relations to approximate the real part of
%the refractive index from the imaginary part
%auxiliary variable
imag_n = imag(n_eff);

%stores the real part of N
n_re_KK = zeros(n_frec,1);

%it is supposed that the frec is equispaced
delta_omega = omega(2)-omega(1);

%Evaluate the Kramers Kronig integral with trapeze rule
%Calculate the first element
term_a = imag_n(2)*omega(2)/(omega(2)^2 - omega(1)^2);
for i = 2:n_frec-1
    %trapeze rule for integration
    term_b        = imag_n(i+1)*omega(i+1)/(omega(i+1)^2 - omega(1)^2);
    n_re_KK(1)  = n_re_KK(1) + term_a + term_b;
    term_a        = term_b;
end;
n_re_KK(1) = 1.0 + delta_omega/pi*n_re_KK(1);

%Calculate the last element
term_a = imag_n(1)*omega(1)/(omega(1)^2-omega(n_frec)^2);
for i = 1:n_frec-2
    %trapeze rule for integration
    term_b = imag_n(i+1)*omega(i+1)/(omega(i+1)^2-omega(n_frec)^2);
    n_re_KK(n_frec)  = n_re_KK(n_frec) + term_a + term_b;
    term_a = term_b;
end;
n_re_KK(n_frec) = 1.0 + delta_omega/pi*n_re_KK(n_frec);

%Calculate the middle elements
for i = 2:n_frec-1; 
    n_re_KK(i) = 0.0;
    term_a = imag_n(1)*omega(1)/(omega(1)^2 - omega(i)^2);
    for j = 1:i-2
        %trapeze rule for integration
        term_b     = imag_n(j+1)*omega(j+1)/(omega(j+1)^2 - omega(i)^2);
        n_re_KK(i) = n_re_KK(i) + term_a + term_b;
        term_a     = term_b;
    end;
    term_a = imag_n(i+1)*omega(i+1)/(omega(i+1)^2-omega(i)^2);
    for j = i+1:n_frec-1
         %trapeze rule for integration
         term_b     = imag_n(j+1)*omega(j+1)/(omega(j+1)^2-omega(i)^2);
         n_re_KK(i) = n_re_KK(i) +  term_a + term_b;
         term_a     = term_b;
     end;
    n_re_KK(i) = 1.0 + delta_omega/pi*n_re_KK(i);
end;
%End apply Kramers-Kronig relations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step - 7 calculate the branch number m
m_branch = round( (n_re_KK - real(n_eff)).*k0*d_eff/(2.0*pi) );

%store the -2, -1, 0, 1, 2 branches for plotting
n_eff_0    = real(n_eff);
n_eff_1    = n_eff_0 + 2.0*pi*1./(k0*d_eff);
n_eff_min1 = n_eff_0 + 2.0*pi*(-1)./(k0*d_eff);
n_eff_2    = n_eff_0 + 2.0*pi*2./(k0*d_eff);
n_eff_min2 = n_eff_0 + 2.0*pi*(-2)./(k0*d_eff);

%Step - 8 calculate the real part of n with the branch
n_eff   = n_eff + 2.0*pi*m_branch./(k0*d_eff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step - 9 calculate eps and mu
eps_eff = conj(n_eff./Z_eff);
mu_eff  = conj(n_eff.*Z_eff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step - 10 plot the results
%plot the S parameters amplitude
figure
     set(gcf,'Color','w');
     plot(frec,S11_abs(ind_frec:n_S11_abs,2),'Color','r','Linewidth',2);
     hold on;
     plot(frec,S21_abs(ind_frec:n_S11_abs,2),'Color','b','Linewidth',2);
     axis([f_min f_max 0 1])
     set(gca,'XTick',f_min:delta_f:f_max);
     set(gca,'YTick',0:0.2:1.0);
     grid on;
     set(gca,'FontSize',14);
     xlabel(['Frequency (',TextUnitFrec, ')'],'FontSize',16);
 	 ylabel('Magnitude of S ','FontSize',16);
     legend('\itS\rm_1_1','\itS\rm_2_1');
     
 %plot the S parameters phase
 figure
     set(gcf,'Color','w');
     plot(frec,S11_phase(ind_frec:n_S11_abs,2),'r','Linewidth',2);
     hold on;
     plot(frec,S21_phase(ind_frec:n_S11_abs,2),'b','Linewidth',2);
     axis([f_min f_max -180 180])
     set(gca,'XTick',f_min:delta_f:f_max);
     set(gca,'YTick',-180:60:180);
     grid on;
     set(gca,'FontSize',14);
     xlabel(['Frequency (',TextUnitFrec, ')'],'FontSize',16);
 	 ylabel('Phase of S ','FontSize',16);
     legend('\itS\rm_1_1','\itS\rm_2_1');

 %plot the wave impedance Z
 figure
     set(gcf,'Color','w');
     plot(frec,real(Z_eff),'r','Linewidth',2);
     hold on;
     plot(frec,imag(Z_eff),'b','Linewidth',2);
     grid on;
     axis([f_min f_max -7 7])
     set(gca,'XTick',f_min:delta_f:f_max);
     set(gca,'YTick',-7:3.5:7);
     set(gca,'FontSize',14);
     xlabel(['Frequency (',TextUnitFrec, ')'],'FontSize',16);
     ylabel('\itZ\rm_{eff}','FontSize',16);
     legend('Re(\itZ\rm_{eff})','Im(\itZ\rm_{eff})');
 
%plot the refractive index N
figure
     set(gcf,'Color','w');
     plot(frec,imag(n_eff),'b','Linewidth',2); %imaginary part
     hold on;
     plot(frec,n_re_KK,'k','Linewidth',2);     %the Kramers Kronig approximation 
     plot(frec,real(n_eff),'r','Linewidth',2); %the real part
     %several branches to see what is happening
     plot(frec,n_eff_min2,':','Color',[0 0.7 0],'Linewidth',2);
     plot(frec,n_eff_min1,':','Color',[0.5 0.5 0.5],'Linewidth',2);
     plot(frec,n_eff_0,':','Color',[0.078431372549020   0.168627450980392   0.549019607843137],'Linewidth',2);
     plot(frec,n_eff_1,':','Color',[0.7 0 0],'Linewidth',2);
     plot(frec,n_eff_2,':','Color',[1 0.6 0.2],'Linewidth',2);
     plot(frec,real(n_eff),'r','Linewidth',2); %the real part again
     axis([f_min f_max -10 10])
     set(gca,'XTick',f_min:delta_f:f_max);
     set(gca,'YTick',-10:5:10);
     grid on;
     set(gca,'FontSize',14);
     xlabel(['Frequency (',TextUnitFrec, ')'],'FontSize',16);
     ylabel('\itN\rm_{eff} ','FontSize',16);
     legend('\kappa_{eff}','\itn\rm^K^K','\itn_{eff}');
  
%plot the branch number m
figure
     set(gcf,'Color','w');
     plot(frec,m_branch,'r','Linewidth',2);
     axis([f_min f_max -2 2])
     grid on;
     set(gca,'XTick',f_min:delta_f:f_max);
     set(gca,'YTick',-2:2);
     set(gca,'FontSize',14);
     xlabel(['Frequency (',TextUnitFrec, ')'],'FontSize',16);
     ylabel('\itm\rm_{branch}','FontSize',16);
 
%plot the effective electric permittivity and magnetic permeability
figure;
    set(gcf,'Color','w');
    %eps
    plot(frec,real(eps_eff),'r','Linewidth',2);
    hold on;
    plot(frec,-imag(eps_eff),'b','Linewidth',2);
    %mu
%     plot(frec,real(mu_eff),'k','Linewidth',2);
%     plot(frec,-imag(mu_eff),'Color',[0 0.7 0],'Linewidth',2);
    axis([f_min f_max -45 40]);
    set(gca,'XTick',f_min:delta_f:f_max);
    set(gca,'YTick',-45:10:40);
    grid on;
    set(gca,'FontSize',14);
    xlabel(['Frequency (',TextUnitFrec, ')'],'FontSize',16);
	ylabel('\epsilon_{eff} , \mu_{eff} ','FontSize',18);
	legend('Re(\epsilon_{eff})','Im(\epsilon_{eff})','Re(\mu_{eff})','Im(\mu_{eff})');
figure;
set(gcf,'Color','w');
 plot(frec,real(mu_eff),'k','Linewidth',2);
 hold on;
    plot(frec,-imag(mu_eff),'Color',[0 0.7 0],'Linewidth',2);
    axis([f_min f_max -5 5]);
    set(gca,'XTick',f_min:delta_f:f_max);
    set(gca,'YTick',-5:1:5);
    grid on;
    set(gca,'FontSize',14);
    xlabel(['Frequency (',TextUnitFrec, ')'],'FontSize',16);
	ylabel('\mu_{eff} ','FontSize',18);
	legend('Re(\mu_{eff})','Im(\mu_{eff})');
