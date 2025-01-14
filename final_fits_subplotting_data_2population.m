%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Plots experimental data at dose of 0.0032, 0.32, 3.2 in top left.     %
% Then plots best-fit for the 2-population model at these 3 doses.      %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all; clc;

%% Read in and store experimental data
data=readtable('pcbi.1007688_reformat_ALL.xlsx');
data=table2array(data);
time=data(1:750,1);
A_32=data(1:750,22:25);
A_032=data(1:750,14:17);
A_0032=data(1:750,6:9);
for i=1:750
    Amean_32(i) = mean(A_32(i,:)/A_32(1,:));
    Amean_032(i) = mean(A_032(i,:)/A_032(1,:));
    Amean_0032(i) = mean(A_0032(i,:)/A_0032(1,:));

    Astd_32(i) = std(A_32(i,:)./A_32(1,:));
    Astd_032(i) = std(A_032(i,:)./A_032(1,:));
    Astd_0032(i) = std(A_0032(i,:)./A_0032(1,:));
end
Amean_32=Amean_32';
Amean_032=Amean_032';
Amean_0032=Amean_0032';
time=time';


%% best parameters
%ORDER OF PARAMS: r_s,d_s, a, r_r, e
best_0032=[0.0532	0.1855	0.2520	0.0532	0.4052];
best_032=[0.0722	0.4012	0.1963	0.0415	0.1742];
best_32=[0.0880	0.6924	0.1048	0.0070	0.0198];

%% running all 4 best fits
[t,y_0032]=ode23s(@(t,x) two_pop_model(t,x,best_0032),time,[1,0]);
[t,y_032]=ode23s(@(t,x) two_pop_model(t,x,best_032),time,[1,0]);
[t,y_32]=ode23s(@(t,x) two_pop_model(t,x,best_32),time,[1,0]);

z=1:length(time);
allpops_0032(z)=y_0032(z,1)+y_0032(z,2);
allpops_032(z)=y_032(z,1)+y_032(z,2);
allpops_32(z)=y_32(z,1)+y_32(z,2);

figure;
set(groot,'defaultAxesLabelFontSize',14) % axes label font size
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.65, 0.85]);
subplot(2,2,1)
plot(time,Amean_0032,'LineWidth',3); hold on; 
plot(time,Amean_032,'--','LineWidth',3);
plot(time,Amean_32,':','LineWidth',3); hold off;
xlabel('Time (hr)','FontSize',16)
ylabel('Mean of Normalized Cell Count','FontSize',16)
legend('0.032 \muM','0.32 \muM','3.2 \muM','FontSize',16,'Location','Northwest')

subplot(2,2,2)
errorbar(time,Amean_0032, Astd_0032,Astd_0032); hold on;
plot(time,allpops_0032,'LineWidth',3, 'Color', 'r');
plot(time, y_0032(:,1),'--g','LineWidth',2)
plot(time, y_0032(:,2),':k','LineWidth',2)
xlabel('Time (hr)','FontSize',16)
ylabel('Normalized Cell Count','FontSize',16)
title('0.032 \muM','FontSize',16)
hold off

subplot(2,2,3)
errorbar(time,Amean_032, Astd_032,Astd_032); hold on;
plot(time,allpops_032,'LineWidth',3, 'Color', 'r');
plot(time, y_032(:,1),'--g','LineWidth',2)
plot(time, y_032(:,2),':k','LineWidth',2)
xlabel('Time (hr)','FontSize',16)
ylabel('Normalized Cell Count','FontSize',16)
title('0.32 \muM','FontSize',16)
hold off

subplot(2,2,4)
errorbar(time,Amean_32, Astd_32,Astd_32); hold on;
plot(time,allpops_32,'LineWidth',3, 'Color', 'r');
plot(time, y_32(:,1),'--g','LineWidth',2)
plot(time, y_32(:,2),':k','LineWidth',2)
xlabel('Time (hr)','FontSize',16)
ylabel('Normalized Cell Count','FontSize',16)
ylim([0,inf])
title('3.2 \muM','FontSize',16)
hold off
legend('Data: Mean \pm Std Deviation','Model: S+R','Model: S','Model: R','FontSize',16)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS                                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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