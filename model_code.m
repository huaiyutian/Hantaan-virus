%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB code to fit the model described in article
%%% The main file for fig4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through smoothing bandwidths (h) and penalty weights (penalty).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Supplied dataset contains the following parameters:
% timePts: the number of data points in the dataset
% startTime: the year of the first data point
% n: the number of time intervals that comprise an annual cycle
% m: the maximum duration of time lag
% AAdensity: a column vector of length timePts with population density data of striped field mouse (Apodemus agrarius, AA)
% RNdensity: a column vector of length timePts with population density data of Norway rat (Rattus norvegicus, RN)
% RFdensity: a column vector of length timePts with population density data of buff-breasted rat (R. flavipectus, RF)
% avg_patch_size: a column vector of length timePts with average patch size data
% temp: a column vector of length timePts with temperature data
% rainfall: a column vector of length timePts with rainfall data

clear; close all; warning off;
dataset = 'dataset.mat'
%seed(123)

% Parameters for fitting the model:
nknots = 10;  % number of knots in the spline fit
h = 1;  % number of smoothing bandwidth (h) 
penalty = -5; % number of penalty weights (penalty) 

%%%%%%%%%%%%%%%%    AA    %%%%%%%%%%%%%%%%%%%
[dates, alpha1, alpha2, alpha3,alpha1hat, alpha2hat, alpha3hat, tau1, tau2, tau3,tau1hat,tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, Y,Yhat,AAdensity1, r2_logAA,  exitflag] = fit_model_AA(dataset, nknots, h, penalty, 1);% 1 indicates these results should be saved in a .mat file
     
figure;
subplot(3,3,1); hold on; ylabel('tau1'); plot(0:(length(tau1)-1), tau1, 'b'); 
subplot(3,3,2); hold on; ylabel('tau2'); plot(0:(length(tau2)-1), tau2, 'b'); 
subplot(3,3,3); hold on; ylabel('tau3'); plot(0:(length(tau3)-1), tau3, 'b'); 
subplot(3,3,4); hold on; plot(logbetarain,'b'); ylabel('betarain');
subplot(3,3,5); hold on; plot(logbetatemp,'b'); ylabel('betatemp');
subplot(3,3,6); hold on; plot(dates, betalongterm,'b'); ylabel('betalongterm'); 
subplot(3,3,7); hold on; plot(0:(length(alpha1)-1), alpha1,'b'); ylabel('alpha1'); 
subplot(3,3,8); hold on; plot(0:(length(alpha2)-1), alpha2,'b'); ylabel('alpha2');
subplot(3,3,9); hold on; plot(0:(length(alpha3)-1), alpha3,'b'); ylabel('alpha3'); 

% get log population density
figure;
subplot(2,1,1); hold on; plot(real(Y), real(Yhat),'b.');xlabel('logAA'); ylabel('logAAhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
subplot(2,1,2);hold on;plot(real(Y),'LineWidth',1);hold on;plot(real(Yhat),'LineWidth',1,'color','red'); legend('logAA','logAAhat');



%%%%%%%%%%%%%%%%    RN    %%%%%%%%%%%%%%%%%%%
[dates, alpha1, alpha2, alpha3, alpha1hat,alpha2hat, alpha3hat,  tau1, tau2, tau3, tau1hat, tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0,Y,Yhat,RNdensity1, r2_logRN, exitflag2] = fit_model_RN(dataset, nknots, h, penalty, 1);% 1 indicates these results should be saved in a .mat file

% plotting the results
figure;
subplot(3,3,1); hold on; ylabel('tau1'); plot(0:(length(tau1)-1), tau1, 'b'); 
subplot(3,3,2); hold on; ylabel('tau2'); plot(0:(length(tau2)-1), tau2, 'b'); 
subplot(3,3,3); hold on; ylabel('tau3'); plot(0:(length(tau3)-1), tau3, 'b'); 
subplot(3,3,4); hold on; plot(logbetarain,'b'); ylabel('betarain');
subplot(3,3,5); hold on; plot(logbetatemp,'b'); ylabel('betatemp');
subplot(3,3,6); hold on; plot(dates, betalongterm,'b'); ylabel('betalongterm'); 
subplot(3,3,7); hold on; plot(0:(length(alpha1)-1), alpha1,'b'); ylabel('alpha1'); 
subplot(3,3,8); hold on; plot(0:(length(alpha2)-1), alpha2,'b'); ylabel('alpha2');
subplot(3,3,9); hold on; plot(0:(length(alpha3)-1), alpha3,'b'); ylabel('alpha3'); 

% get log population density
figure;
subplot(2,1,1); hold on; plot(real(Y), real(Yhat),'b.');xlabel('logRN'); ylabel('logRNhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
subplot(2,1,2);hold on;plot(real(Y),'LineWidth',1);hold on;plot(real(Yhat),'LineWidth',1,'color','red'); legend('logRN','logRNhat');



%%%%%%%%%%%%%%%%    RF    %%%%%%%%%%%%%%%%%%%
[dates, alpha1, alpha2, alpha3, alpha1hat, alpha2hat, alpha3hat, tau1, tau2, tau3, tau1hat, tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0, Y,Yhat,RFdensity1, r2_logRF, exitflag3] = fit_model_RF(dataset, nknots, h, penalty, 1);% 1 indicates these results should be saved in a .mat file
     
% plotting the results
figure;
subplot(3,3,1); hold on; ylabel('tau1'); plot(0:(length(tau1)-1), tau1, 'b'); 
subplot(3,3,2); hold on; ylabel('tau2'); plot(0:(length(tau2)-1), tau2, 'b'); 
subplot(3,3,3); hold on; ylabel('tau3'); plot(0:(length(tau3)-1), tau3, 'b'); 
subplot(3,3,4); hold on; plot(logbetarain,'b'); ylabel('betarain');
subplot(3,3,5); hold on; plot(logbetatemp,'b'); ylabel('betatemp');
subplot(3,3,6); hold on; plot(dates, betalongterm,'b'); ylabel('betalongterm'); 
subplot(3,3,7); hold on; plot(0:(length(alpha1)-1), alpha1,'b'); ylabel('alpha1'); 
subplot(3,3,8); hold on; plot(0:(length(alpha2)-1), alpha2,'b'); ylabel('alpha2');
subplot(3,3,9); hold on; plot(0:(length(alpha3)-1), alpha3,'b'); ylabel('alpha3'); 

% get log population density
figure;
subplot(2,1,1); hold on; plot(real(Y), real(Yhat),'b.');xlabel('logRF'); ylabel('logRFhat'); v = axis; vMin = min(v(1), v(3)); vMax = max(v(2), v(4)); axis([vMin vMax vMin vMax]); line(vMin:vMax, vMin:vMax);
subplot(2,1,2);hold on;plot(real(Y),'LineWidth',1);hold on;plot(real(Yhat),'LineWidth',1,'color','red'); legend('logRF','logRFhat');



