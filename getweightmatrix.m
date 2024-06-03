function [W] = getweightmatrix(h, size)
% returns a truncated Gaussian kernel weight matrix

W = zeros(size, size);
for i=1:size
	for j=1:size
        if(abs(i-j)<=h)
            W(i,j) = exp((-0.5)*(((i-j)/h)^2));
        end
    end
end

% adjusting the weights such that the sum of the weights in each row is 1
for i=1:size
	sumofRow = sum(W(i,:));
    W(i,:) = W(i,:)/sumofRow;
end
