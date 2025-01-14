%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Fitting 3-population induced resistance model to dose-specific time   %
% course data for COLO858 cells treated Vemurafenib                     %
% Authors: Jana Gevertz and Samantha Prosperi                           %
% Updated: 11/19/2024                                                   %
%                                                                       %
% Normalized dose-specific time course data is fit to the model output  %
% S+Q+R using a multistart fmincon algorithm. The model used is:        %
%  S' = r_s*S-d*S*(1-exp(-1*gamma_1*t))-q*S*(1-exp(-1*gamma_2*t))       %
%  Q' = q*S*(1-exp(-1*gamma_2*t))-res*Q                                 %
%  R' = res*Q+r_r*R-e*d*R*(1-exp(-1*gamma_1*t))                         %
% with no pre-existing resistance. Parameters fit per dose are:         %
% 1) r_S = growth rate of sensitive (S) cells                           %
% 2) d = death rate of S                                                %
% 3) q= rate of S --> quiescient (Q)                                    %
% 4) res = rate at which Q transition to resistance (R)                 %
% 5) r_R = growth rate of R (assumed to be <= r_S)                      %
% 6) e=epsilon: where death rate of R d_R is = e*d_S (assumed to be     %
%    between 0 and 1 so R are less responsive to drug than S)           %
% The delay terms on death (gamma_1) and resistance (gamma_2) are both  %
% fixed at 0.01, and thus not fit.                                      %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load in the time and data
clear all; close all; clc; tic;
rng(1); % fixed seed for debugging
options = optimoptions(@fmincon,'MaxFunctionEvaluations', 3*10000,...
    'Display','notify'); %fmincon options

%% Read in and store experimental data
data=readtable('pcbi.1007688_reformat_ALL.xlsx');
data=table2array(data);
time=data(1:750,1); %time 0 to 100
fulltime=data(:,1);
prompt = "What treatment do you want to fit?\n Enter enter 3.2, 1, 0.32, .1, .032, or 0 for uM treatment: ";
treat_flag = input(prompt); % Stores value answered to prompt
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

%% User chooses cost function to minize: paper uses L1
prompt2 = "\nWhat cost function do you want to minimize?\n Enter 1 for L1, 2 for L2: ";
cost_flag = input(prompt2); % Stores value answered to prompt
if ((cost_flag~=1)&&(cost_flag~=2))
    fprintf('Entered cost_flag = %f but can only be 1 or 2 - exiting\n',...
        cost_flag);
    stop
end

%% Set output directory
cl = clock; clN = 0;
for ii = 2:5
    clN = floor(100*clN + cl(ii));
end
if cost_flag == 1
    path = ['Output_fmincon_', num2str(treat_flag), ...
        'uM_L1_AutoTracking_3pop_' , num2str(clN)];
elseif cost_flag == 2
    path = ['Output_fmincon_', num2str(treat_flag), ...
        'uM_L2_AutoTracking_3pop' , num2str(clN)];
end
if exist(path, 'dir') ~= 7
    mkdir(path)
end
diary_file = [path '/diary_output'];
diary(diary_file)

%% Average truncated data at each time point across replicates and normalize
num_timePts = size(A,1);
Amean = zeros(1,num_timePts);
Astd  = zeros(1,num_timePts);
for i=1:num_timePts
    Amean(i) = mean(A(i,:)/A(1,:));
    Astd(i) = std(A(i,:)/A(1,:));
end
time = time'; 

%% Plot the raw data
figure; 
plot(time,Amean,'ro'); 
xlabel("time(hrs)")
ylabel("raw cell count")
title(treat_flag+ "uM treatment of COLO858 cells with Vem")

%% Parameter fitting
% ORDER OF PARAMS: r_s, d, q, res, r_r, e 
% Parameter bounds for multi start and fmincon
%    r_s     d    q     res   r_r   e 
lb = [0.001  1e-4 1e-4  1e-4  1e-4  1e-4];
ub = [0.15   0.7  0.6   0.6   0.6   1];

% Multi start setup
numParams = length(lb);
numSobols = 1000; % how many random points to sample
n_skip = 1000; n_leap = 0; % Parameters needed by sobolset 
uniform_sobol = sobolset(numParams,'Skip',n_skip,'Leap',n_leap);
uniform_sobol = net(uniform_sobol,numSobols);
% Rescale parameters to be in range [lb,ub]
uniform_sobol_scaled = zeros(numSobols,numParams);
for i = 1:numParams
    uniform_sobol_scaled(:,i) = (ub(i) - lb(i))*uniform_sobol(:,i)+lb(i);
end

% Inequality bounds for fmincon: (r_R)-(r_S)<=0 => r_S>=r_R
A_ineq=[-1 0  0 0 1 0 ];
b      = 0;

% Objective function to minimize
fun = @(z)objective(z,time,Amean,cost_flag);

% Multi-start minimization across Mpool parallel pools
Mpools = 4; % number of parallel pools
param = zeros(numSobols,numParams);
fits = zeros(1,numSobols);
exitflag = zeros(1,numSobols);
parfor (i = 1:numSobols, Mpools)
    %call fmincon
    if mod(i,10) == 0 
        fprintf('Up to multistart fitting #%d\n',i)
    end
    [param(i,:), fits(i),exitflag(i)] = fmincon(fun,...
        uniform_sobol_scaled(i,:),A_ineq,b,[],[],lb,ub,[],options);
end

%% Sort best fit per multistart from best to worst and save and then 
%% save all with cost function within 20% of optimal cost
[fit_sorted, index] = sort(fits); 

% Best fit 
best_fit_objective_list(1) = fit_sorted(1);
best_fit_params_list(1,:)  = param(index(1),:);
perc_dif(1)=0;

% Finding all fits within 5% of the best fit
best_counter=2;
for i=2:length(fit_sorted)
    %percent difference from best objective function
    perc_dif(i)=(fit_sorted(i)-fit_sorted(1))/fit_sorted(1);

    if perc_dif(i) < 0.05 % within 5% from optimal
        best_fit_objective_list(best_counter) = fit_sorted(i);
        best_fit_params_list(best_counter,:)  = param(index(i),:);
        best_counter = best_counter+1;
    end
end
toc;


%% Plotting the best fit model
[t,y] = ode23s(@(t,x) three_pop_model(t,x,best_fit_params_list(1,:)),...
    time,[1 0 0]);
z=1:length(time);
allpops(z)=y(z,1)+y(z,2)+y(z,3);
figure;
%plot(time,Amean, 'Color', 'b', 'LineWidth',1.5); hold on;
errorbar(time,Amean,Astd/sqrt(4),"o"); hold on;
plot(time,allpops,'LineWidth',2);%, 'Color', 'r');
plot(time, y(:,1))
plot(time, y(:,2))
plot(time, y(:,3))
%plot(time, replicates, 'Color','k')
hold off;
xlabel('Time (hr)')
ylabel('Normalized Cell Count')
title("3 Population model of resistance of " +treat_flag + "uM Vemurafenib in COLO858")
legend('Data Average','S+R+Q Best Fit', 'S','Q', 'R','location','best')
subtitle("Objective Function: " + best_fit_objective_list(1) )
exportgraphics(gcf,[path '/'  ,'BestFit.png']);
saveas(gcf,[path '/' , 'BestFit.fig']);

%% Save output to path
save([path '/', 'output.mat'],"time","Amean","Astd","param",...
    "fit_sorted","index","best_fit_params_list","best_fit_objective_list"); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS                                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = objective(p,time,ydata,flag) % Cost to minimize
    % Solves model assuming all cells initiall sensitive
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