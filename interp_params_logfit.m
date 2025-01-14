%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use linear interpolated parameters to predict treatment response at   %
% dose of 0.1 and 1 muM.                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load in the time and data
clear all; close all; clc;
data=readtable('pcbi.1007688_reformat_ALL.xlsx');
data=table2array(data);
time=data(1:750,1); %time 0 to 100
fulltime=data(:,1);
A_1=data(1:750,18:21);
A_01=data(1:750,10:13);

for i=1:length(time)
    Amean_1(i) = mean(A_1(i,:)/A_1(1,:));
    Astd_1(i) = std(A_1(i,:)./A_1(1,:))/2;

    Amean_01(i)=mean(A_01(i,:)/A_01(1,:));
    Astd_01(i) = std(A_01(i,:)./A_01(1,:))/2;
end


%ORDER OF PARAMS: r_s,d_s, a, r_r, e
%PIECEWISE LINEAR
interp_01 = [0.0626	0.2922	0.2244	0.04741025475	0.2481177276];
interp_1  = [0.08	0.5453	0.151	0.02442767425	0.07720520814];

[t_01,y_01]=ode23s(@(t,x) two_pop_model(t,x,interp_01),time,[1,0]);
[t_1,y_1]=ode23s(@(t,x) two_pop_model(t,x,interp_1),time,[1,0]);

z=1:length(time);
allpops_01(z)=y_01(z,1)+y_01(z,2);
allpops_1(z)=y_1(z,1)+y_1(z,2);

figure
set(groot,'defaultAxesFontSize',14) % axes font size
set(groot,'defaultAxesLabelFontSize',14) % axes label font size
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.65, 0.55]);

subplot(1,2,1)
errorbar(time,Amean_01,Astd_01,Astd_01); hold on;
plot(time,allpops_01,'LineWidth',3, 'Color', 'r');
plot(time, y_01(:,1),'--g','LineWidth',2)
plot(time, y_01(:,2),':k','LineWidth',2)
xlabel('Time (hr)','FontSize',16)
ylabel('Mean of Normalized Cell Count','FontSize',16)
title('0.1 \muM','FontSize',16)
legend('Data: Mean \pm Std Deviation','Model: S+R','Model: S','Model: R',...
    'Location','NorthWest','FontSize',16)
hold off

subplot(1,2,2)
hold on;
errorbar(time,Amean_1,Astd_1,Astd_1)
plot(time,allpops_1,'LineWidth',3, 'Color', 'r');
plot(time, y_1(:,1),'--g','LineWidth',2)
plot(time, y_1(:,2),':k','LineWidth',2)
xlabel('Time (hr)','FontSize',16)
ylabel('Mean of Normalized Cell Count','FontSize',16)
title('1 \muM','FontSize',16)
hold off
box on;  % Ensure a black box surrounds the plot


%% FUNCTIONS
function xp=two_pop_model(t,x,p)
    %ORDER OF PARAMS: r_s,d_1, rev, r_r, d_2 (e)
    r_s=p(1);
    d=p(2);
    a=p(3);
    r_r=p(4);
    e=p(5);
    gamma_1=0.01;
    gamma_2=0.01;
    
    S=x(1);
    R=x(2);
    
    xp(1)=r_s*S-d*S*(1-exp(-1*gamma_1*t))-a*S*(1-exp(-1*gamma_2*t));
    xp(2)=a*S*(1-exp(-1*gamma_2*t))+r_r*R-e*d*R*(1-exp(-1*gamma_1*t));
    
    xp=xp';
end