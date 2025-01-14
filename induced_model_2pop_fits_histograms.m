%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Display best-fits and suboptimal parameter sets from 2-population     %
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
load("best_params.mat");
best_0032 = best_pars_0032(1,:); 
best_032  = best_pars_032(1,:);
best_32   = best_pars_32(1,:);

%% Run and then display all best fits on one plot
[t,y_0032] = ode23s(@(t,x) two_pop_model(t,x,best_0032),time,[1,0]);
[t, y_032] = ode23s(@(t,x) two_pop_model(t,x,best_032),time,[1,0]);
[t,  y_32]  = ode23s(@(t,x) two_pop_model(t,x,best_32),time,[1,0]);
z = 1:length(time);
allpops_0032(z) = y_0032(z,1)+y_0032(z,2);
allpops_032(z)  = y_032(z,1)+y_032(z,2);
allpops_32(z)   = y_32(z,1)+y_32(z,2);

figure;
subplot(3,1,1)
errorbar(time,Amean_0032, Astd_0032,Astd_0032); hold on;
plot(time,allpops_0032,'LineWidth',2, 'Color', 'r');
plot(time, y_0032(:,1), 'Color', 'g')
plot(time, y_0032(:,2), 'Color','magenta')
xlabel('time (hr)')
ylabel('Norm Cell Count')
title('0.032uM')
hold off
subplot(3,1,2)
errorbar(time,Amean_032, Astd_032,Astd_032); hold on;
plot(time,allpops_032,'LineWidth',2, 'Color', 'r');
plot(time, y_032(:,1), 'Color', 'g')
plot(time, y_032(:,2), 'Color','magenta')
xlabel('time (hr)')
ylabel('Norm Cell Count')
title('0.32uM')
hold off
subplot(3,1,3)
errorbar(time,Amean_32, Astd_32,Astd_32); hold on;
plot(time,allpops_32,'LineWidth',2, 'Color', 'r');
plot(time, y_32(:,1), 'Color', 'g')
plot(time, y_32(:,2), 'Color','magenta')
xlabel('time (hr)')
ylabel('Norm Cell Count')
title('3.2uM')
hold off

%% Histograms of best parameters
% Convert epsilon parameter to d_R: d_r = d_S*epsilon
best_pars_0032(:,5) = best_pars_0032(:,2).*best_pars_0032(:,5);
best_pars_032(:,5)  = best_pars_032(:,2).*best_pars_032(:,5);
best_pars_32(:,5)   = best_pars_32(:,2).*best_pars_32(:,5);
for param = 1:2 % These parameters are in increasing order
    figure;
    set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.45, 0.35]);
    t = tiledlayout(1,3,'TileSpacing','compact');
    bgAx = axes(t,'XTick',[],'YTick',[],'Box','off');
    bgAx.Layout.TileSpan = [1 3];
    % First histogram
    bin_min = 0.99*min(best_pars_0032(:,param));
    bin_max = 1.01*max(best_pars_0032(:,param));
    bin_mean = mean(best_pars_0032(:,param));
    coeff_var = std(best_pars_0032(:,param))./mean(best_pars_0032(:,param));
    bin_edges = linspace(bin_min, bin_max, 30);
    ax1 = axes(t);
    histogram(best_pars_0032(:,param),bin_edges,'Normalization','probability');     
    ax1.Box = 'off';
    xlim(ax1,[bin_min bin_max])
    %bin_middle = (bin_min+bin_max)/2;
    xticks([bin_min bin_mean]);
    my_labels = {};
    my_labels{2} = sprintf('%.3f\\newline(CV=%.3f)', bin_mean,coeff_var);
    my_labels{1} = sprintf('%.3f', bin_min);
    xticklabels(my_labels)

    ax1 = gca();
    if param == 1
        %set(get(ax1,'YLabel'),'String','r_S','FontSize',18);
        ylabel('r_S','FontSize',18)
    elseif param == 2
        ylabel('d_S','FontSize',18)
    elseif param == 3
        ylabel('\alpha','FontSize',18)
    elseif param == 4
        ylabel('r_R','FontSize',18)
    else
        ylabel('d_R','FontSize',18)
    end
    if param == 1
        leg = legend('0.032\muM');
    end
    
    % Second histogram
    bin_min = 0.99*min(best_pars_032(:,param));
    bin_max = 1.01*max(best_pars_032(:,param));
    bin_mean = mean(best_pars_032(:,param));
    coeff_var = std(best_pars_032(:,param))./mean(best_pars_032(:,param));
    bin_edges = linspace(bin_min, bin_max, 30);
    ax2 = axes(t);
    ax2.Layout.Tile = 2;
    histogram(best_pars_032(:,param),bin_edges,'Normalization','probability',...
        'FaceColor',[0.85 0.33 0.10]);
    ax2.YAxis.Visible = 'off';
    ax2.Box = 'off';
    xlim(ax2,[bin_min bin_max])
    %bin_middle = (bin_min+bin_max)/2;
    xticks(bin_mean);
    my_labels = {};
    my_labels{1} = sprintf('%.3f\\newline(CV=%.3f)', bin_mean,coeff_var);
    xticklabels(my_labels)
    if param == 1
        leg = legend('0.32\muM');
    end
    
    % Third histogram
    bin_min = 0.99*min(best_pars_32(:,param));
    bin_max = 1.01*max(best_pars_32(:,param));
    bin_mean = mean(best_pars_32(:,param));
    coeff_var = std(best_pars_32(:,param))./mean(best_pars_32(:,param));
    bin_edges = linspace(bin_min, bin_max, 30);
    ax3 = axes(t);
    ax3.Layout.Tile = 3;
    histogram(best_pars_32(:,param),bin_edges,'Normalization','probability',...
        'FaceColor',[0.93 0.69 0.13]);
    ax3.YAxis.Visible = 'off';
    ax3.Box = 'off';
    xlim(ax3,[bin_min bin_max])
    %bin_middle = (bin_min+bin_max)/2;
    xticks([bin_mean bin_max]);
    my_labels = {};
    my_labels{1} = sprintf('%.3f\\newline(CV=%.3f)', bin_mean,coeff_var);
    my_labels{2} = sprintf('%.3f', bin_max);
    xticklabels(my_labels)
    if param == 1
        leg = legend('3.2\muM');
    end
    
    % Link the axes
    linkaxes([ax1 ax2 ax3], 'y')

    % Set ticks and limits
    ylim([0,0.85]);
    yticks([0.25 0.5 0.75]);
    ax = findall(gcf, 'type', 'axes');
    set(ax,'fontsize', 16) 

    %fname_fig = ['histogram' num2str(param)];
    %saveas(gcf,[fname_fig,'.fig'])
    %saveas(gcf,[fname_fig,'.png']);
end

for param = 3:5 % These parameters are in decreasing order
    figure;
    set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.45, 0.35]);
    t = tiledlayout(1,3,'TileSpacing','compact');
    bgAx = axes(t,'XTick',[],'YTick',[],'Box','off');
    bgAx.Layout.TileSpan = [1 3];
    % First histogram
    bin_min = 0.99*min(best_pars_32(:,param));
    bin_max = 1.01*max(best_pars_32(:,param));
    bin_mean = median(best_pars_32(:,param));
    coeff_var = std(best_pars_32(:,param))./mean(best_pars_32(:,param));

    bin_edges = linspace(bin_min, bin_max, 30);
    ax1 = axes(t);
    histogram(best_pars_32(:,param),bin_edges,'Normalization','probability',...
        'FaceColor',[0.93 0.69 0.13]);     
    ax1.Box = 'off';
    xlim(ax1,[bin_min bin_max])
    xticks([bin_min bin_mean]);
    my_labels = {};
    my_labels{2} = sprintf('%.3f\\newline(CV=%.3f)', bin_mean,coeff_var);
    my_labels{1} = sprintf('%.3f', bin_min);
    xticklabels(my_labels)
    ax1 = gca();
    if param == 1
        %set(get(ax1,'YLabel'),'String','r_S','FontSize',18);
        ylabel('r_S','FontSize',18)
    elseif param == 2
        ylabel('d_S','FontSize',18)
    elseif param == 3
        ylabel('\alpha','FontSize',18)
    elseif param == 4
        ylabel('r_R','FontSize',18)
    else
        ylabel('d_R','FontSize',18)
    end

    % Second histogram
    bin_min = 0.99*min(best_pars_032(:,param));
    bin_max = 1.01*max(best_pars_032(:,param));
    bin_mean = median(best_pars_032(:,param));
    coeff_var = std(best_pars_032(:,param))./mean(best_pars_032(:,param));

    bin_edges = linspace(bin_min, bin_max, 30);
    ax2 = axes(t);
    ax2.Layout.Tile = 2;
    histogram(best_pars_032(:,param),bin_edges,'Normalization','probability',...
        'FaceColor',[0.85 0.33 0.10]);
    ax2.YAxis.Visible = 'off';
    ax2.Box = 'off';
    xlim(ax2,[bin_min bin_max])
    % bin_middle = (bin_min+bin_max)/2;
    xticks(bin_mean);
    my_labels = {};
    my_labels{1} = sprintf('%.3f\\newline(CV=%.3f)', bin_mean,coeff_var);
    xticklabels(my_labels)

    % Third histogram
    bin_min = 0.99*min(best_pars_0032(:,param));
    bin_max = 1.01*max(best_pars_0032(:,param));
    bin_mean = median(best_pars_0032(:,param));
    coeff_var = std(best_pars_0032(:,param))./mean(best_pars_0032(:,param));
   
    bin_edges = linspace(bin_min, bin_max, 30);
    ax3 = axes(t);
    ax3.Layout.Tile = 3;
    histogram(best_pars_0032(:,param),bin_edges,'Normalization','probability');
    ax3.YAxis.Visible = 'off';
    ax3.Box = 'off';
    xlim(ax3,[bin_min bin_max])
    %bin_middle = (bin_min+bin_max)/2;
    xticks([bin_mean bin_max]);
    my_labels = {};
    my_labels{1} = sprintf('%.3f\\newline(CV=%.3f)', bin_mean,coeff_var);
    my_labels{2} = sprintf('%.3f', bin_max);
    xticklabels(my_labels)

    % Link the axes
    linkaxes([ax1 ax2 ax3], 'y')

    % Set ticks and limits
    ylim([0,0.85]);
    yticks([0.25 0.5 0.75]);
    ax = findall(gcf, 'type', 'axes');
    set(ax,'fontsize', 16) 

    %fname_fig = ['histogram' num2str(param)];
    %saveas(gcf,[fname_fig,'.fig'])
    %saveas(gcf,[fname_fig,'.png']);
end

%% 2 POP MODEL
function xp=two_pop_model(t,x,p)
    %ORDER OF PARAMS: r_s,d_1, rev, r_r, d_2
    r_s=p(1);
    d=p(2);
    res=p(3);
    r_r=p(4);
    e=p(5);
    gamma_1=0.01;
    gamma_2=0.01;
    
    S=x(1);
    R=x(2);
    xp(1)=r_s*S-d*S*(1-exp(-1*gamma_1*t))-res*S*(1-exp(-1*gamma_2*t));
    xp(2)=res*S*(1-exp(-1*gamma_2*t))+r_r*R-e*d*R*(1-exp(-1*gamma_1*t));
    
    xp=xp';
end
