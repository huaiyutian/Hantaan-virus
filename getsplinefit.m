function [tau1hat, exitflag] = getsplinefit(tau1, nknots, m, penalty)

x = (0:m)';

% cubic spline
Xfit(:,1) = ones(length(x),1);
Xfit(:,2) = x;
Xfit(:,3) = x.^2;
Xfit(:,4) = x.^3;

% setting the knots
spacing = length(x)/(nknots + 1);  % evenly spaced knots
for i=1:nknots
    knots(i) = floor(spacing*i);
end

loc = 5;
for i=1:nknots
    Xfit(:,loc) = (max(x-knots(i),0)).^3;
    loc = loc + 1;
end

XTX = Xfit'*Xfit;

% A is the constraint matrix -> constraints here are monotonic decrease and positivity
for i=1:m
    A(i,:) = Xfit(i+1,:)-Xfit(i,:);
end
A(m+1,:) = -1*(Xfit(m+1,:));
b = zeros(m+1,1);
f = -1*(Xfit'*tau1)';

% setting the penalty
if penalty == -Inf
    H = XTX;
else
    D = zeros(size(nknots+4));
    for i=1:nknots
        D(i+4, i+4) = 10^penalty;
    end
    
    H = XTX + D;
end   

options = optimset('MaxIter',1000, 'Display', 'off');
[betahat, fval, exitflag, output] = quadprog(H, f, A, b, [], [], [], [], [], options);

tau1hat = Xfit*betahat;
