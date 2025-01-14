%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Display best-fits and suboptimal parameter sets from 3-population     %
% induced resistance model.                                             %
% Authors: Jana Gevertz and Samantha Prosperi                           %
% Updated: 11/19/2024                                                   %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% Read in and store experimental data
data   = readtable('pcbi.1007688_reformat_ALL.xlsx');
data   = table2array(data);
time   = data(1:750,1);
A_32   = data(1:750,22:25);
A_032  = data(1:750,14:17);
A_0032 = data(1:750,6:9);
num_timePts = size(A_32,1);
Amean_32    = zeros(1,num_timePts);
Amean_032   = zeros(1,num_timePts);
Amean_0032  = zeros(1,num_timePts);
for i=1:750
    Amean_32(i)   = mean(A_32(i,:)/A_32(1,:));
    Amean_032(i)  = mean(A_032(i,:)/A_032(1,:));
    Amean_0032(i) = mean(A_0032(i,:)/A_0032(1,:));

    Astd_32(i)   = std(A_32(i,:)./A_32(1,:));
    Astd_032(i)  = std(A_032(i,:)./A_032(1,:));
    Astd_0032(i) = std(A_0032(i,:)./A_0032(1,:));
end
Amean_32   = Amean_32';
Amean_032  = Amean_032';
Amean_0032 = Amean_0032';
time=time';


%% Load best fit parameter values
load("best_params_3pop.mat");
best_0032 = best_pars_0032(1,:);
best_032  = best_pars_032(1,:);
best_32   = best_pars_32(1,:);

%% Run and then display all best fits on one plot
[t,y_0032] = ode23s(@(t,x) three_pop_model(t,x,best_0032),time,[1,0,0]);
[t, y_032] = ode23s(@(t,x) three_pop_model(t,x,best_032),time,[1,0,0]);
[t,  y_32] = ode23s(@(t,x) three_pop_model(t,x,best_32),time,[1,0,0]);
z=1:length(time);
allpops_0032(z)=y_0032(z,1)+y_0032(z,2)+y_0032(z,3);
allpops_032(z)=y_032(z,1)+y_032(z,2)+y_032(z,3);
allpops_32(z)=y_32(z,1)+y_32(z,2)+y_32(z,3);

figure;
subplot(3,1,1)
errorbar(time,Amean_0032, Astd_0032,Astd_0032); hold on;
plot(time,allpops_0032,'LineWidth',2, 'Color', 'r');
plot(time, y_0032(:,1), 'Color', 'g')
plot(time, y_0032(:,2), 'Color','k')
plot(time, y_0032(:,3), 'Color','magenta')
xlabel('time (hr)')
ylabel('Norm Cell Count')
title('0.032uM')
hold off
subplot(3,1,2)
errorbar(time,Amean_032, Astd_032,Astd_032); hold on;
plot(time,allpops_032,'LineWidth',2, 'Color', 'r');
plot(time, y_032(:,1), 'Color', 'g')
plot(time, y_032(:,2), 'Color','k')
plot(time, y_032(:,3), 'Color','magenta')
xlabel('time (hr)')
ylabel('Norm Cell Count')
title('0.32uM')
hold off
subplot(3,1,3)
errorbar(time,Amean_32, Astd_32,Astd_32); hold on;
plot(time,allpops_32,'LineWidth',2, 'Color', 'r');
plot(time, y_32(:,1), 'Color', 'g')
plot(time, y_32(:,2), 'Color','k')
plot(time, y_32(:,3), 'Color','magenta')
xlabel('time (hr)')
ylabel('Norm Cell Count')
title('3.2uM')
hold off

%% Histograms of best parameters
% Convert epsilon parameter to d_R: d_r = d_S*epsilon
best_pars_0032(:,6) = best_pars_0032(:,2).*best_pars_0032(:,6);
best_pars_032(:,6)  = best_pars_032(:,2).*best_pars_032(:,6);
best_pars_32(:,6)   = best_pars_32(:,2).*best_pars_32(:,6);

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.65, 0.95]);
for i = 1:6
    binEdges0032_rS = linspace(min(best_pars_0032(:,i)), max(best_pars_0032(:,i)), 30);
    binEdges032_rS = linspace(min(best_pars_032(:,i)), max(best_pars_032(:,i)), 30);
    binEdges32_rS = linspace(min(best_pars_32(:,i)), max(best_pars_32(:,i)), 30);

    subplot(3,2,i)
    histogram(best_pars_0032(:,i),binEdges0032_rS,'Normalization', 'probability'); hold on;
    histogram(best_pars_032(:,i),binEdges032_rS,'Normalization', 'probability')
    histogram(best_pars_32(:,i),binEdges32_rS,'Normalization', 'probability')
    hold off;

    if i == 1
        ylabel('r_S','fontsize', 16) 
        legend('0.032\muM','0.32\muM','3.2\muM')
    elseif i == 2
        ylabel("d_S",'fontsize', 16) 
    elseif i == 3
        ylabel("q",'fontsize', 16) 
    elseif i == 4
        ylabel("\beta",'fontsize', 16) 
    elseif i == 5
        ylabel("r_R",'fontsize', 16) 
    else
        ylabel("d_R",'fontsize', 16) 
    end
end

%% 3 POP MODEL
function xp = three_pop_model(t,x,p)
    %ORDER OF PARAMS: r_s, d_S, q, beta, r_r, d_R
    r_s=p(1);
    d=p(2);
    q=p(3);
    beta=p(4);
    r_r=p(5);
    e=p(6);
    gamma_1=0.01;
    gamma_2=0.01;
    
    S=x(1);
    Q=x(2);
    R=x(3);    
    
    xp(1) = r_s*S-d*S*(1-exp(-1*gamma_1*t))-q*S*(1-exp(-1*gamma_2*t));
    xp(2) = q*S*(1-exp(-1*gamma_2*t))-beta*Q;
    xp(3) = beta*Q+r_r*R-e*d*R*(1-exp(-1*gamma_1*t));
    
    xp=xp';
end
