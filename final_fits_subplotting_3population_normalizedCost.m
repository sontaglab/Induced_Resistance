%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Plots best-fit to dose of 0.0032, 0.32, 3.2 for 3-population model.   %
% Also computes normalized value of absolute error for both the 2-pop   %
% and 3-pop model. For the 2-pop model, this value is also computed at  %
% the interpolated doses of 0.1 and 1.                                  %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% Read in data
data=readtable('pcbi.1007688_reformat_ALL.xlsx');
data=table2array(data);
time=data(1:750,1);
A_32=data(1:750,22:25);
A_032=data(1:750,14:17);
A_0032=data(1:750,6:9);
A_1=data(1:750,18:21);
A_01=data(1:750,10:13);

for i=1:750
    Amean_32(i) = mean(A_32(i,:)/A_32(1,:));
    Amean_032(i) = mean(A_032(i,:)/A_032(1,:));
    Amean_0032(i) = mean(A_0032(i,:)/A_0032(1,:));
    Amean_1(i) = mean(A_1(i,:)/A_1(1,:));
    Amean_01(i)=mean(A_01(i,:)/A_01(1,:));

    Astd_32(i) = std(A_32(i,:)./A_32(1,:));
    Astd_032(i) = std(A_032(i,:)./A_032(1,:));
    Astd_0032(i) = std(A_0032(i,:)./A_0032(1,:));
    Astd_1(i) = std(A_1(i,:)./A_1(1,:))/2;
    Astd_01(i) = std(A_01(i,:)./A_01(1,:))/2;
end

Amean_32   = Amean_32';
Amean_032  = Amean_032';
Amean_0032 = Amean_0032';
Amean_1    = Amean_1';
Amean_01   = Amean_01';
time       = time';

%% 2-population model: best fit parameters
%ORDER OF PARAMS: r_s,d_s, a, r_r, e
best_0032 = [0.0532	0.1855	0.2520	0.0532	0.4052];
best_032  = [0.0722	0.4012	0.1963	0.0415	0.1742];
best_32   = [0.0880	0.6924	0.1048	0.0070	0.0198];
[t_0032, y_0032, cost_0032_2pop] = objective(best_0032,time,Amean_0032,2);
[t_032, y_032, cost_032_2pop]    = objective(best_032,time,Amean_032,2);
[t_32, y_32,  cost_32_2pop]      = objective(best_32,time,Amean_32,2);

%% 3-population model: best fit parameters
%ORDER OF PARAMS: r_s,d_s, a, r_r, e
best_0032_3pop = [0.0530	0.1650	0.2691	0.5310	0.0530	0.4535];
best_032_3pop  = [0.0716	0.3709	0.2330	0.1689	0.0414	0.1877];
best_32_3pop   = [0.08812307266	 0.692579592	0.106663255 	0.320499141	...
                  0.006861658324	0.01934840644];
[t_0032_3pop, y_0032_3pop, cost_0032_3pop] = objective(best_0032_3pop,time,Amean_0032,3);
[t_032_3pop, y_032_3pop, cost_032_3pop]    = objective(best_032_3pop,time,Amean_032,3);
[t_32_3pop, y_32_3pop,  cost_32_3pop]      = objective(best_32_3pop,time,Amean_32,3);

figure;
set(groot,'defaultAxesFontSize',14) % axes font size
set(groot,'defaultAxesLabelFontSize',14) % axes label font size
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.65, 0.85]);
subplot(2,2,1)
errorbar(time,Amean_0032, Astd_0032,Astd_0032); hold on;
plot(t_0032_3pop,y_0032_3pop(:,1)+y_0032_3pop(:,2)+y_0032_3pop(:,3),...
    'LineWidth',3, 'Color', 'r');
plot(t_0032_3pop, y_0032_3pop(:,1),'--g','LineWidth',2)
plot(t_0032_3pop, y_0032_3pop(:,2),':k','LineWidth',2)
plot(t_0032_3pop, y_0032_3pop(:,3),'-.','Color','magenta','LineWidth',2)
xlabel('Time (hr)','FontSize',16)
ylabel('Normalized Cell Count','FontSize',16)
title('0.032\muM','FontSize',16)
hold off

subplot(2,2,2)
errorbar(time,Amean_032, Astd_032,Astd_032); hold on;
plot(t_032_3pop,y_032_3pop(:,1)+y_032_3pop(:,2)+y_032_3pop(:,3),...
    'LineWidth',3, 'Color', 'r');
plot(t_032_3pop, y_032_3pop(:,1),'--g','LineWidth',2)
plot(t_032_3pop, y_032_3pop(:,2),':k','LineWidth',2)
plot(t_032_3pop, y_032_3pop(:,3),'-.','Color','magenta','LineWidth',2)
xlabel('Time (hr)','FontSize',16)
ylabel('Normalized Cell Count','FontSize',16)
title('0.32\muM','FontSize',16)
hold off

subplot(2,2,3)
errorbar(time,Amean_32, Astd_32,Astd_32); hold on;
plot(t_32_3pop,y_32_3pop(:,1)+y_32_3pop(:,2)+y_32_3pop(:,3),...
    'LineWidth',3, 'Color', 'r');
plot(t_32_3pop, y_32_3pop(:,1),'--g','LineWidth',2)
plot(t_32_3pop, y_32_3pop(:,2),':k','LineWidth',2)
plot(t_32_3pop, y_32_3pop(:,3),'-.','Color','magenta','LineWidth',2)
xlabel('Time (hr)','FontSize',16)
ylabel('Normalized Cell Count','FontSize',16)
title('3.2\muM','FontSize',16)
hold off
legend('Data: Mean \pm Std Deviation','Model: S+R','Model: S','Model: Q','Model: R','FontSize',16)

subplot(2,2,4)
cost_all = [cost_0032_2pop cost_0032_3pop; cost_032_2pop cost_032_3pop; ...
    cost_32_2pop cost_32_3pop];
cost_all = round(cost_all,2)
str = {'0.032 \muM','0.32 \muM','3.2 \muM'};
b = bar(cost_all);
set(gca,'XTickLabel',str,'XTick',1:numel(str)); %,'FontSize',14)
xlabel('Dose','FontSize',16)
ylabel('Sum of Relative Absolute Error','FontSize',16)
% Will include numerical height of graph for dataset 1
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',14) 
% Will include numerical height of graph for dataset 1
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',14) 
legend('2 pop model','3 pop model','FontSize',16,'Location','NorthWest')

%% 2-population model: extrapolated parameters
%ORDER OF PARAMS: r_s,d_s, a, r_r, e
interp_01 = [0.0626	0.2922	0.2244	0.04741025475	0.2481177276];
interp_1  = [0.08	0.5453	0.151	0.02442767425	0.07720520814];
cost_01_2pop = objective(interp_01,time,Amean_01,2)
cost_1_2pop  = objective(interp_1,time,Amean_1,2)


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

function xp=three_pop_model(t,x,p)
    %ORDER OF PARAMS: r_s,d_1, rev, r_r, d_2
    r_s=p(1);
    d=p(2);
    q=p(3);
    res=p(4);
    r_r=p(5);
    e=p(6);
    gamma_1=0.01;
    gamma_2=0.01;
    
    S=x(1);
    Q=x(2);
    R=x(3);
    
    xp(1) = r_s*S-d*S*(1-exp(-1*gamma_1*t))-q*S*(1-exp(-1*gamma_2*t));
    xp(2) = q*S*(1-exp(-1*gamma_2*t))-res*Q;
    xp(3) = res*Q+r_r*R-e*d*R*(1-exp(-1*gamma_1*t));
    
    xp=xp';
end

function [t, y, cost] = objective(p,time,ydata,flag) % Cost to minimize
    % Solves model assuming all cells initiall sensitive
    z=1:length(time);
    if flag == 2 % 2-pop
        [t, y] = ode23s(@(t,x) two_pop_model(t,x,p),time,[1,0]);
        modeldata(z)=y(z,1)+y(z,2); 
    elseif flag == 3 %3-pop
        [t, y] = ode23s(@(t,x) three_pop_model(t,x,p),time,[1,0,0]);
        modeldata(z)=y(z,1)+y(z,2)+y(z,3); 
    end
    %[t y(:,1) y(:,2)]
    % size(t)
    % stop

    cost =sum(abs(ydata-modeldata')./ydata);
end

