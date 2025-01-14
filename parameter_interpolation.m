%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Performs piecewise linear interpolation of best-fit parameters at     %
% dose of 0.032, 0.32 and 3.2 to predict parameter values at dose of    %
% 0.1 and 1 muM.                                                        %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

%% Data: updated EXCEPT at doses of 0.1 and 1
drug_all = [0.032 0.1 0.32 1 3.2];
drug = drug_all(1:2:end);
drug_validate = drug_all(2:2:end);
rS_all = [0.0532 0.0630 0.0722 0.0539 0.0880];
dS_all = [0.1855 0.3394 0.4012 0.2862 0.6921];
alpha_all = [0.2520 0.2154 0.1963 0.1378 0.1049];
rR_all = [0.0532 0.0519 0.0415 0.0190 0.0070];
dR_all = [0.0751 0.0819 0.0699 0.0340 0.0137];

%% Divide rs into training and validation and use linear interpolation
rS = rS_all(1:2:end);
rS_interp = interp1(drug,rS,drug_validate); % On (dose, parameter)
rS_interp_log = interp1(log10(drug),rS,log10(drug_validate)); % On (log10(dose),parameter)
figure; 
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.5, 0.75]);
subplot(1,2,1)
semilogx(drug_all,rS_all,'o-','LineWidth',4); hold on;
semilogx(drug_validate,rS_interp,'*r','MarkerSize',10);
legend('True Best-Fit Parameter','Extrapolated Parameter','FontSize',16,...
    'Location','NorthWest');
title('Piecewise Linear Interpolation','FontSize',16);
xlabel('Dose (\muM)','FontSize',14);
ylabel('r_s','FontSize',14);
hold off;
subplot(1,2,2)
semilogx(drug_all,rS_all,'o-','LineWidth',4); hold on;
semilogx(drug_validate,rS_interp_log,'*r','MarkerSize',10);
legend('True Best-Fit Parameter','Extrapolated Parameter','FontSize',16,...
    'Location','NorthWest');
title('Piecewise Linear Interpolation (Log of Dose)','FontSize',16);
xlabel('Dose (\muM)','FontSize',14);
ylabel('r_s','FontSize',14);
hold off;

%% Divide ds into training and validation and use linear interpolation
dS = dS_all(1:2:end);
dS_interp = interp1(drug,dS,drug_validate); % On (dose, parameter)
dS_interp_log = interp1(log10(drug),dS,log10(drug_validate));  % On (log10(dose),parameter)
figure; 
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.5, 0.75]);
subplot(1,2,1)
semilogx(drug_all,dS_all,'o-','LineWidth',4); hold on;
semilogx(drug_validate,dS_interp,'*r','MarkerSize',10);
legend('True Best-Fit Parameter','Extrapolated Parameter','FontSize',16);
title('Piecewise Linear Interpolation','FontSize',16);
xlabel('Dose (\muM)','FontSize',14);
ylabel('d_s','FontSize',14);
hold off;
subplot(1,2,2)
semilogx(drug_all,dS_all,'o-','LineWidth',4); hold on;
semilogx(drug_validate,dS_interp_log,'*r','MarkerSize',10);
legend('True Best-Fit Parameter','Extrapolated Parameter','FontSize',16);
title('Piecewise Linear Interpolation (Log of Dose)','FontSize',16);
xlabel('Dose (\muM)','FontSize',14);
ylabel('d_s','FontSize',14);
hold off;

%% Divide alpha into training and validation and use linear interpolation
alpha = alpha_all(1:2:end);
alpha_interp = interp1(drug,alpha,drug_validate); % On (dose, parameter)
alpha_interp_log = interp1(log10(drug),alpha,log10(drug_validate)); % On (log10(dose),parameter)
figure; 
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.5, 0.75]);
subplot(1,2,1)
semilogx(drug_all,alpha_all,'o-','LineWidth',4); hold on;
semilogx(drug_validate,alpha_interp,'*r','MarkerSize',10);
legend('True Best-Fit Parameter','Extrapolated Parameter','FontSize',16);
title('Piecewise Linear Interpolation','FontSize',16);
xlabel('Dose (\muM)','FontSize',14);
ylabel('\alpha','FontSize',14);
hold off;
subplot(1,2,2)
semilogx(drug_all,alpha_all,'o-','LineWidth',4); hold on;
semilogx(drug_validate,alpha_interp_log,'*r','MarkerSize',10);
legend('True Best-Fit Parameter','Extrapolated Parameter','FontSize',16);
title('Piecewise Linear Interpolation (Log of Dose)','FontSize',16);
xlabel('Dose (\muM)','FontSize',14);
ylabel('\alpha','FontSize',14);
hold off;

%% Divide rR into training and validation and use linear interpolation
rR = rR_all(1:2:end);
rR_interp = interp1(drug,rR,drug_validate); % On (dose, parameter)
rR_interp_log = interp1(log10(drug),rR,log10(drug_validate)); % On (log10(dose),parameter)
figure; 
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.5, 0.75]);
subplot(1,2,1)
semilogx(drug_all,rR_all,'o-','LineWidth',4); hold on;
semilogx(drug_validate,rR_interp,'*r','MarkerSize',10);
legend('True Best-Fit Parameter','Extrapolated Parameter','FontSize',16);
title('Piecewise Linear Interpolation','FontSize',16);
xlabel('Dose (\muM)','FontSize',14);
ylabel('rR','FontSize',14);
hold off;
subplot (1,2,2)
semilogx(drug_all,rR_all,'o-','LineWidth',4); hold on;
semilogx(drug_validate,rR_interp_log,'*r','MarkerSize',10);
legend('True Best-Fit Parameter','Extrapolated Parameter','FontSize',16);
title('Piecewise Linear Interpolation (Log of Dose)','FontSize',16);
xlabel('Dose (\muM)','FontSize',14);
ylabel('r_R','FontSize',14);
hold off;

%% Divide dR into training and validation and use linear interpolation
dR = dR_all(1:2:end);
dR_interp = interp1(drug,dR,drug_validate); % On (dose, parameter)
dR_interp_log = interp1(log10(drug),dR,log10(drug_validate)); % On (log10(dose),parameter)
figure; 
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.5, 0.75]);
subplot(1,2,1)
semilogx(drug_all,dR_all,'o-','LineWidth',4); hold on;
semilogx(drug_validate,dR_interp,'*r','MarkerSize',10);
legend('True Best-Fit Parameter','Extrapolated Parameter','FontSize',16);
title('Piecewise Linear Interpolation','FontSize',16);
xlabel('Dose (\muM)','FontSize',14);
ylabel('d_R','FontSize',14);
hold off;
subplot(1,2,2)
semilogx(drug_all,dR_all,'o-','LineWidth',4); hold on;
semilogx(drug_validate,dR_interp_log,'*r','MarkerSize',10);
legend('True Best-Fit Parameter','Extrapolated Parameter','FontSize',16);
title('Piecewise Linear Interpolation (Log of Dose)','FontSize',16);
xlabel('Dose (\muM)','FontSize',14);
ylabel('d_R','FontSize',14);
hold off;

eps = 0.01;
figure; 
set(groot,'defaultAxesFontSize',14) % axes font size
set(groot,'defaultAxesLabelFontSize',14) % axes label font size
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.65, 0.85]);
subplot(2,3,1)
h1 = semilogx(drug,rS,'o:','LineWidth',2','MarkerSize',6); hold on;
semilogx(drug_validate,rS_interp_log,'*r','MarkerSize',10);
set(h1, 'markerfacecolor', get(h1, 'color')); 
xlabel('Dose (\muM)','FontSize',16);
ylabel('r_S','FontSize',16);
xlim([min(drug) max(drug)])
text(drug_validate(1)+eps,rS_interp_log(1),num2str(rS_interp_log(1),4),'FontSize',14)
text(drug_validate(2)+eps,rS_interp_log(2),num2str(rS_interp_log(2),4),'FontSize',14)
hold off;

subplot(2,3,2)
h1 = semilogx(drug,dS,'o:','LineWidth',2','MarkerSize',6); hold on;
semilogx(drug_validate,dS_interp_log,'*r','MarkerSize',10);
set(h1, 'markerfacecolor', get(h1, 'color')); 
xlabel('Dose (\muM)','FontSize',16);
ylabel('d_S','FontSize',16);
xlim([min(drug) max(drug)])
text(drug_validate(1)+eps,dS_interp_log(1),num2str(dS_interp_log(1),4),'FontSize',14)
text(drug_validate(2)+eps,dS_interp_log(2),num2str(dS_interp_log(2),4),'FontSize',14)
hold off;

subplot(2,3,3)
h1 = semilogx(drug,alpha,'o:','LineWidth',2','MarkerSize',6); hold on;
semilogx(drug_validate,alpha_interp_log,'*r','MarkerSize',10);
set(h1, 'markerfacecolor', get(h1, 'color')); 
xlabel('Dose (\muM)','FontSize',16);
ylabel('\alpha','FontSize',16);
xlim([min(drug) max(drug)])
text(drug_validate(1)+eps,alpha_interp_log(1),num2str(alpha_interp_log(1),4),'FontSize',14)
text(drug_validate(2)+eps,alpha_interp_log(2),num2str(alpha_interp_log(2),4),'FontSize',14)
hold off;

subplot(2,3,4)
h1 = semilogx(drug,rR,'o:','LineWidth',2','MarkerSize',6); hold on;
semilogx(drug_validate,rR_interp_log,'*r','MarkerSize',10);
set(h1, 'markerfacecolor', get(h1, 'color')); 
xlabel('Dose (\muM)','FontSize',16);
ylabel('r_R','FontSize',16);
xlim([min(drug) max(drug)])
text(drug_validate(1)+eps,rR_interp_log(1),num2str(rR_interp_log(1),4),'FontSize',14)
text(drug_validate(2)+eps,rR_interp_log(2),num2str(rR_interp_log(2),4),'FontSize',14)
hold off;

subplot(2,3,5)
h1 = semilogx(drug,dR,'o','LineWidth',1','MarkerSize',6); hold on;
semilogx(drug_validate,dR_interp_log,'*r','MarkerSize',10);
h2 = semilogx(drug,dR,':','LineWidth',2); 
set(h1, 'markerfacecolor', get(h1, 'color')); 
set(h2, 'color', get(h1, 'color')); 
xlabel('Dose (\muM)','FontSize',16);
ylabel('d_R','FontSize',16);
xlim([min(drug) max(drug)])
legend('Best-Fit','Interpolated','FontSize',16);
text(drug_validate(1)+eps,dR_interp_log(1),num2str(dR_interp_log(1),4),'FontSize',14)
text(drug_validate(2)+eps,dR_interp_log(2),num2str(dR_interp_log(2),4),'FontSize',14)
hold off;
