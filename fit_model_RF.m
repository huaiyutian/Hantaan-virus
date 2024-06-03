function [dates, alpha1, alpha2, alpha3,alpha1hat, alpha2hat, alpha3hat, tau1, tau2, tau3, tau1hat, tau2hat, tau3hat, logbetarain,logbetatemp, betalongterm, Y0,Y,Yhat,RFdensity1, r2_logRF, exitflag3] = fit_model_RF(dataset, nknots, h, penalty, saveBool)
% h=1;penalty=-5;nknots = 10;   
load(dataset);

m=12;
% generating the dates of the time series
dates(1) = startTime;
for i=2:timePts
    dates(i) = dates(i-1) + 1/n;
end

% using spline inpolatation for AA
nonNaNIndices = find(~isnan(AAdensity));
AAdensity = interp1(nonNaNIndices, AAdensity(nonNaNIndices), 1:length(AAdensity), 'spline');
AAdensity = AAdensity';
% using spline inpolatation for RN
nonNaNIndices = find(~isnan(RNdensity));
RNdensity = interp1(nonNaNIndices, RNdensity(nonNaNIndices), 1:length(RNdensity), 'spline');
RNdensity = RNdensity';
% using spline inpolatation for RF
nonNaNIndices = find(~isnan(RFdensity));
RFdensity = interp1(nonNaNIndices, RFdensity(nonNaNIndices), 1:length(RFdensity), 'spline');
RFdensity = RFdensity';


% % during time points when the number of cases is negative, replace the cases with 0 cases (there should not be many of these data points)
locs1 = find(AAdensity <0); AAdensity(locs1)=0;
locs2 = find(RNdensity <0); RNdensity(locs2)=0;
locs3 = find(RFdensity <0); RFdensity(locs3)=0;

AAdensity = (AAdensity+0.1).*1.5;

% generating the delta matrix
delta = zeros(timePts, n);
for row=1:timePts
    seas = mod(row, n);
    if seas == 0
        seas = n;
    end
    delta(row,seas) = 1;
end

% taking off the last column to make the number of columns in the delta matrix only be of length n-1
delta(:,n) = [];

% generating the vector log(RF(t))
logRF = log(RFdensity);
% replace all infinite values with 0
logRF(isinf(logRF)) = 0;
    
% generating the vector log(RF(t+1))
Y = zeros(timePts,1);
for i=1:(timePts-1)
	Y(i,1) = log(RFdensity(i+1,1));
end
Y(isinf(Y)) = 0;

% generating the vector RF(t+1)
RFdensity1 = zeros(timePts,1);
for i=1:(timePts-1)
	RFdensity1(i,1) = RFdensity(i+1,1);
end

%generating the vector rainfall(t-2)
rain2 = zeros(timePts,1);
for i=1:(timePts-2)
    rain2(i+2,1) = rainfall(i,1);
end

%generating the vector temp(t-2)
temp2 = zeros(timePts,1);
for i=1:(timePts-2)
    temp2(i+2,1) = temp(i,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% generating the matrix of previous AAdensity 
Imatrix = zeros(timePts, m+1);
Imatrix(:,1) = AAdensity;
prev = 1;	
for backcounter=2:(m+1)
	for row=2:timePts
		Imatrix(row,prev+1) = Imatrix(row-1,prev);
	end
	prev = prev+1;
end

% generating the matrix of previous RNdensity 
Rmatrix = zeros(timePts, m+1);
Rmatrix(:,1) = RNdensity;
prev = 1;	
for backcounter=2:(m+1)
	for row=2:timePts
		Rmatrix(row,prev+1) = Rmatrix(row-1,prev);
	end
	prev = prev+1;
end

% generating the matrix of previous RFdensity 
Fmatrix = zeros(timePts, m+1);
Fmatrix(:,1) = RFdensity;
prev = 1;	
for backcounter=2:(m+1)
	for row=2:timePts
		Fmatrix(row,prev+1) = Fmatrix(row-1,prev);
	end
	prev = prev+1;
end

% generating the matrix of previous patch effect
Pmatrix = zeros(timePts, m+1);
Pmatrix(:,1) = avg_patch_size;
prev = 1;	
for backcounter=2:(m+1)
	for row=2:timePts
		Pmatrix(row,prev+1) = Pmatrix(row-1,prev);
	end
	prev = prev+1;
end


% cutting out the first m+1 data points and the last data point for reasons of incompleteness
dates(1:(m)) = [];
rain2([1:(m)],:) = [];
temp2([1:(m)],:) = [];
delta([1:(m)],:) = [];
logRF([1:(m)],:) = [];
Y([1:(m)],:) = [];
Imatrix([1:(m)],:) = [];
Fmatrix([1:(m)],:) = [];
Rmatrix([1:(m)],:) = [];
Pmatrix([1:(m)],:) = [];
RFdensity([1:(m)],:) = [];
RFdensity1(1:m,:) = [];


% resetting the length of the timePts
timePts = timePts-m;

% generating the long-term beta smoother
Wblt = getweightmatrix(h, timePts); 

% generating an identity matrix of the same dimensions as the long-term beta smoother
Iblt = eye(size(Wblt));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We now have all the necessary vectors and matrices that we'll keep on re-using.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepping for the first weighted least-squares regression:
% getting U
for i=1:timePts
    U(i,:) = -Pmatrix(i,:).*Imatrix(i,:);
end

% getting V
for i=1:timePts
    V(i,:) = -Pmatrix(i,:).*Rmatrix(i,:);
end

% getting W
for i=1:timePts
    W(i,:) = -Pmatrix(i,:).*Fmatrix(i,:);
end

% getting delta_rain
delta_rain=zeros(size(delta));
for i=1:timePts
    delta_rain(i,:) = delta(i,:).*rain2(i);
end

% getting delta_temp
delta_temp=zeros(size(delta));
for i=1:timePts
    delta_temp(i,:) = delta(i,:).*temp2(i);
end

X = real([delta_rain, delta_temp, -Imatrix, -Rmatrix, -Fmatrix, U, V, W]);
Y0=Y-logRF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIRST WEIGHTED-LEAST-SQUARES REGRESSION
theta = real(inv(X'*(Iblt-Wblt)*X)*X'*(Iblt-Wblt)*Y0);
logbetarain = [theta([1:n-1],1); 0];
logbetatemp = [theta([n:2*n-2],1); 0];
alpha1 = theta([2*n-1:2*n+m-1],1);
alpha2 = theta([2*n+m:2*n+2*m],1);
alpha3 = theta([2*n+2*m+1:2*n+3*m+1],1);
tau1 = theta([2*n+3*m+2:2*n+4*m+2]);
tau2 = theta([(2*n+4*m+3):(2*n+5*m+3)]);
tau3 = theta([(2*n+5*m+4):(2*n+6*m+4)]);

[tau1hat, exitflag] = getsplinefit(tau1, nknots, m, penalty);
[tau2hat, exitflag2] = getsplinefit2(tau2, nknots, m, penalty);
[tau3hat, exitflag3] = getsplinefit3(tau3, nknots, m, penalty);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBSEQUENT WEIGHTED LEAST-SQUARES REGRESSIONS: The backfitting algorithm
   
for k =1:50   
    X = [delta_rain, delta_temp,-Imatrix, -Rmatrix, -Fmatrix];
    Y1 = Y0-U*tau1hat-V*tau2hat-W*tau3hat;
    theta = real(inv(X'*(Iblt-Wblt)*X)*X'*(Iblt-Wblt)*Y1);
    logbetarain = [theta([1:(n-1)],1); 0];
    logbetatemp = [theta([n:(2*n-2)],1); 0];
    alpha1 = theta([2*n-1:2*n+m-1],1);
    alpha2 = theta([2*n+m:2*n+2*m],1);
    alpha3 = theta([2*n+2*m+1:2*n+3*m+1],1);
    
    rho = Y1-X*theta;
    Z = [U,V,W];
    
    phi = inv(Z'*(Iblt-Wblt)*Z)*Z'*(Iblt-Wblt)*rho;
    phi1 = phi([1:(m+1)],1);
    phi2 = phi([(m+2):(2*m+2)],1);
    phi3 = phi([(2*m+3):(3*m+3)],1);
    
    tau1 = tau1hat+phi1;
    tau2 = tau2hat+phi2;
    tau3 = tau3hat+phi3;
    
    [tau1hat, exitflag] = getsplinefit(real(tau1), nknots, m, penalty);
    [tau2hat, exitflag2] = getsplinefit2(real(tau2), nknots, m, penalty);
    [tau3hat, exitflag3] = getsplinefit3(real(tau3), nknots, m, penalty);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOW PROCEED WITH FINAL FIT 
tau1hat = real(tau1hat);
tau2hat = real(tau2hat);
tau3hat = real(tau3hat);
alpha1hat = real(alpha1hat);
alpha2hat = real(alpha2hat);
alpha3hat = real(alpha3hat);

% Y = f1hat + f2hat, where f1hat is the parametric part of the equation and f2hat is the nonparametric part of the equation
X = [delta_rain, delta_temp,-Imatrix, -Rmatrix, -Fmatrix];
Y1 = real(Y0-U*tau1hat-V*tau2hat-W*tau3hat);
theta = real(inv(X'*(Iblt-Wblt)*X)*X'*(Iblt-Wblt)*Y1);
f1hat = X*theta; 
logbetarain = [theta([1:(n-1)],1); 0];
logbetatemp = [theta([n:(2*n-2)],1); 0];
alpha1 = theta([2*n-1:2*n+m-1],1);
alpha2 = theta([2*n+m:2*n+2*m],1);
alpha3 = theta([2*n+2*m+1:2*n+3*m+1],1);

% f1hat is the parametric component of the fit
% f2hat is the non-parametric component of the fit

f2hat = Wblt*(Y1-f1hat); % local constant fit
betalongterm = f2hat; 

% entire fit.
% Yhat = logRF+f1hat+f2hat+U*tau1hat+V*tau2hat+W*tau3hat; 
% resids = Y-Yhat;
GRhat = f1hat+f2hat+U*tau1hat+V*tau2hat+W*tau3hat;
resids = Y0-GRhat;

% % calculate RFhat
RFhat = zeros(timePts, 1);
RFhat(16,1)=RFdensity1(16,1);
for i = 17:timePts
    RFhat(i,1)=RFhat(i-1,1)*exp(GRhat(i,1));
end

% calculate Yhat
Yhat = zeros(timePts, 1);
Yhat(1,1)=Y(1,1);
for i = 2:timePts-1
    Yhat(i,1)=Yhat(i-1)+GRhat(i,1);
end 
% RFhat = exp(Yhat);

% R^2
SSR1 = sum((Yhat-mean(Y)).^2);
SSTO1 = sum((Y-mean(Y)).^2);
r2_logRF = SSR1/SSTO1;

SSR2 = sum((RFhat-mean(RFdensity1)).^2);
SSTO2 = sum((RFdensity1-mean(RFdensity1)).^2);
r2_RF = SSR2/SSTO2;

SSR3 = sum((GRhat-mean(Y0)).^2);
SSTO3 = sum((Y0-mean(Y0)).^2);
r2_GR = SSR3/SSTO3;

if saveBool
    filename = strcat('Opt_Results_RF_', dataset);
    save(filename, 'nknots', 'm', 'n', 'penalty', 'h', 'dates', 'Wblt', 'Iblt', 'X', 'Y0','Y','Yhat','RFdensity1', 'alpha1', 'alpha1hat', 'alpha2', 'alpha2hat', 'alpha3', 'alpha3hat', 'tau1', 'tau1hat','tau2', 'tau2hat','tau3', 'tau3hat','logbetarain','logbetatemp','betalongterm', 'f2hat');
    
end