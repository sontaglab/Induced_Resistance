%This code will model fit a specfiic treatment of COLO858 cells treated
%with 1uM Venurafenib

%Evaluating the possibility of modeling inducible resistance in this cell
%line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description of the model:
%
% - All cells start as sensitive to treatment (S)
% -After exposure to treatment, 2 possible fates:
% 1. Cell death (sensitivity)-->goes first to a “dying” compartment at rate d_1, and goes to null at
%    some rate d_2
% 3. Inducibly resistance (R) (continues to divide, this fate is reversible)
%
% Fitting X= S+ Q + R
%
% S' = r_s*S-d*S*(1-exp(-1*gamma_1*t))-q*S*(1-exp(-1*gamma_2*t));
% Q' = q*S*(1-exp(-1*gamma_2*t))-res*Q;
% R' = res*Q+r_r*R-e*d*R*(1-exp(-1*gamma_1*t));
%
%
% with these variables to fit:
% r_s=growth rate of sensitive cells
% d=rate that cell death
% q= rate of S --> Q
% res=rate of Q-->R
% r_r=growth rate of sensitive cells (should be less than r_s)
% e=epsilon: from 0 to 1 of how fast R cells can die compared to S
% gamma_1: rate of delayed drug onset's impact on death (fixed at 0.01)
% gamma_2: rate of delayed drug onset's imapct on S-->R (fixed at 0.01)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% load in the time and data
clear all;
close all;
clc;

tic;

%% 
data=readtable('pcbi.1007688_reformat_ALL.xlsx');
data=table2array(data);

time=data(1:750,1); %time 0 to 100

fulltime=data(:,1);

prompt2= "What treatment do you want to fit?\n Enter enter 3.2, 1, 0.32, .1, .032, or 0 for uM treatment: ";
treat_flag = input(prompt2); % Stores value answered to prompt
switch treat_flag
    case 3.2
        A=data(1:750,22:25);
        A_all=data(:,22:25);
    case 1
        A=data(1:750,18:21);
        A_all=data(:,18:21);
    case 0.32
        A=data(1:750,14:17);
        A_all=data(:,14:17);
    case 0.1
        A=data(1:750,10:13);
        A_all=data(:,10:13);
    case 0.032
        A=data(1:750, 6:9);
        A_all=data(:,6:9);
    case 0
        A=data(1:750, 2:5);
        A_all=data(:,2:5);
    otherwise
        fprintf('Entered a treatment that is unavailable');
        stop
end

for i=1:length(A)
    Amean(i) = mean(A(i,:)/A(1,:));
    Astd(i) = std(A(i,:)/A(1,:));
end

for i=1:length(A_all)
    all_mean(i)=mean(A_all(i,:));
end


num_timePts=length(A);
time=transpose(time);


%% setup output path
cost_flag=1;
%set up exporting all of the figure files created
rng(1); % fixed seed for debugging
cl = clock; clN = 0;
for ii = 2:5
    clN = floor(100*clN + cl(ii));
end
if cost_flag == 1
    path = ['Output_fmincon_', num2str(treat_flag), 'uM_L1_AutoTracking_3pop_' , num2str(clN)];
elseif cost_flag == 2
    path = ['Output_fmincon_', num2str(treat_flag), 'uM_L2_AutoTracking_3pop' , num2str(clN)];
end
if exist(path, 'dir') ~= 7
    mkdir(path)
end
diary_file = [path '/diary_output'];
diary(diary_file)
%% Plot the raw data

figure(1)
%plot(time,Amean)%,Astd)
hold on
plot(time,Amean,'ro')
xlabel("time(hrs)")
ylabel("raw cell count")
title(treat_flag+ "uM treatment of COLO858 cells with Vem")
legend('errorbars','mean','location','best')
exportgraphics(gcf,[path '/', 'Visualization.png']);
saveas(gcf,[path '/','Visualization.fig']);

%% initialize fmincon

options = optimoptions(@fmincon,'MaxFunctionEvaluations', 3*10000);

%ALIGN THESE


%bounds for multi start
%    r_s   d  q   res   r_r    e 
lb=[0.001 1e-4 1e-4 1e-4 1e-4 1e-4];
ub=[0.15 0.7 0.6 0.6 0.6 1];


%multi start setup
numParams = length(lb);
numSobols = 1000;% how many random points to sample
n_skip = 1000; n_leap = 0; % Parameters needed by sobolset - not very important
uniform_sobol = sobolset(numParams,'Skip',n_skip,'Leap',n_leap);
uniform_sobol = net(uniform_sobol,numSobols);
% Rescale parameters to be in range [lb,ub]
for i = 1:numParams
    uniform_sobol_scaled(:,i) = (ub(i) - lb(i))*uniform_sobol(:,i)+lb(i);
end
fits= zeros(1,numSobols);

%% par for to run fmincon
%set inequality bounds on the parameters
A_ineq=[-1 0  0 0 1 0 ];
b=0;
%should evaluate so Ax<=b is a bound on solving the system
%simplifies to (r_r)-(r_s)<0, which is true when r_s>r_r

% bounds for fmincon

%    r_s   d     q   res   r_r    e 
lb2=[0.001 1e-6 1e-6 1e-6 1e-6 1e-6];
ub2=[0.15 0.7 0.6 0.6 0.6 1];

%call set up function before loop
std=1;
fun = @(z)objective(z,time,Amean,cost_flag,Astd);

M=10; %number of computers to send to--> 16 is the max for me locally

parfor (i = 1:numSobols, M)
    %call fmincon
    [param(i,:), fits(i),exitflag(i)] = fmincon(fun,uniform_sobol_scaled(i,:),A_ineq,b,[],[],lb2,ub2,[],options);
    exitflag(i)
end
%saving ALL param fits to path
save([path '/', 'param_fits_all.mat'],"param")
%%
[fit_sorted, index] = sort(fits); % Sort by best fit
%saving sorted objective functions and indexes to path
save([path '/', 'sorted_obj_fxn.mat'],"fit_sorted")
save([path '/', 'sorting_indexes.mat'],"index")
best_fit = zeros(1,numParams);
%% 
best_counter=2;
%setting top of best 20% as the best fit
best_fit_objective_list(1)=fit_sorted(1);
best_fit_params_list(1,:)=param(index(1),:);

perc_dif(1)=0;

for q=2:length(fit_sorted)
    %percent difference from best objective function
    perc_dif(q)=(fit_sorted(q)-fit_sorted(1))/fit_sorted(1);

    if perc_dif(q) < 0.05 %if less than 20%, save it separately
        best_fit_objective_list(best_counter)=fit_sorted(q);
        best_fit_params_list(best_counter,:)=param(index(q),:);
        save([path '/', 'top_0.05_param.mat'],"best_fit_params_list")
        save([path '/', 'top_0.05_obj_fxn.mat'],"best_fit_objective_list")
        best_counter=best_counter+1;
    end
end

%% 
best_fit_objective=fit_sorted(1);
best_fit_params=param(index(1),:);

save([path '/', 'top_param.mat'],"best_fit_params")
save([path '/', 'top_obj_fxn.mat'],"best_fit_objective")

%% Plotting the best fit model

figure()
[t,y] = ode23s(@(t,x) three_pop_model(t,x,best_fit_params),time,[1 0 0]);
z=1:length(time);
allpops(z)=y(z,1)+y(z,2)+y(z,3);
allpops=allpops;
plot(time,Amean, 'Color', 'b', 'LineWidth',1.5); hold on;
%errorbar(time,Amean,Astd/sqrt(4),"o"); hold on;
plot(time,allpops,'LineWidth',2);%, 'Color', 'r');
plot(time, y(:,1))
plot(time, y(:,2))
plot(time, y(:,3))
%plot(time, replicates, 'Color','k')
hold off;
xlabel('time (hr)')
ylabel('Normalized Cell Count')
title("3 Population model of resistance of " +treat_flag + "uM Vemurafenib in COLO858")
legend('Data Average','S+R+Q Best Fit', 'S','Q', 'R','location','best')
subtitle("Objective Function: " + best_fit_objective )
exportgraphics(gcf,[path '/'  ,'BestFit.png']);
saveas(gcf,[path '/' , 'BestFit.fig']);

toc;
%% FUNCTIONS%%

%cost function (keep option for both L1 and L2)
function cost = objective(p,time,ydata,flag,st_dev)
cost = 0;
[t, y] = ode23s(@(t,x) three_pop_model(t,x,p),time,[1,0,0]);

%initial conditions: fixed based on data (could fit as parameters though)
z=1:length(time);
modeldata(z)=y(z,1)+y(z,2)+y(z,3); 

if flag == 1 % L1
    cost =sum(abs(ydata-modeldata));%./st_dev);
elseif flag == 2 % L2
    cost =sum(((ydata-modeldata).^2));%./(st_dev.^2));
end
end
%plotting

%model

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

N=S+R+Q;


xp(1) = r_s*S-d*S*(1-exp(-1*gamma_1*t))-q*S*(1-exp(-1*gamma_2*t));
xp(2) = q*S*(1-exp(-1*gamma_2*t))-res*Q;
xp(3) = res*Q+r_r*R-e*d*R*(1-exp(-1*gamma_1*t));


xp=xp';
end