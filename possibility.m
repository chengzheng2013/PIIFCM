function [p,sigma] = possibility(IMG,epsilon,K)

IMG = imfilter(IMG, fspecial('gaussian', 5,5), 'symmetric', 'conv');
IMG_size = size(IMG);

[up, left, right, down, upleft, upright, downleft, downright] = eigNeig(IMG); 
w = fspecial('average',3);
IMG_mean = imfilter(IMG,w,'same');

IMG_mean = IMG_mean(3:end-2,3:end-2);
mid = IMG(3:end-2,3:end-2);
var = ((upleft-IMG_mean).^2  + (up-IMG_mean).^2+...
      (upright-IMG_mean).^2  + (left-IMG_mean).^2+...
      (mid - IMG_mean).^2    + (right-IMG_mean).^2 +...
      (downleft-IMG_mean).^2 + (down-IMG_mean).^2+...
      (downright-IMG_mean).^2) / 9;

IMG_var = zeros(IMG_size);
IMG_var(3:end-2,3:end-2) = var;
%compute the mean and variance of each pixel in a 3x3 window 
mu_IMG = mean(mean(IMG_var));
sigma_IMG = std(std(IMG_var));

%set the threshold to press the influence of noise
threshold = mu_IMG + K*sigma_IMG;
IMG_1 = IMG_var > threshold;

%repeat the anomaly detection process the locate the position
max_iters = 10;
[row,col] = find(IMG_1 == 1);
X = [col,row];
for i = 1:max_iters
    [M,~] = size(X);

    mu1 = mean(X,1);
    sigma = zeros(2,2);
    for j = 1:M
        sigma = sigma + (X(j,:)-mu1)'*(X(j,:)-mu1);
    end
    sigma = sigma/M;
    p = mvnpdf(X,mu1,sigma);
    outliers = find(p < epsilon);
    if (outliers)
        X(outliers,:) = [];
    else
        break;
    end
end


sigma2 = sigma;
sigma2(1,1) = 1.2*sigma2(1,1);
sigma2(2,2) = 1.2*sigma2(2,2);
% if(sigma2(1,1)>sigma2(2,2))
%     sigma2(2,2) = 1.25*sigma2(2,2);
% else
%     sigma2(1,1) = 1.25*sigma2(1,1);
% end

[X1,X2] = meshgrid(1:1:IMG_size(2),1:1:IMG_size(1));
p = mvnpdf([X1(:),X2(:)],mu1,sigma2);
p = reshape(p,size(X1));
p = p+1e-300;

p_max = max(max(p));
p = 1/p_max*p;
p = reshape(p,IMG_size(1)*IMG_size(2),1);
end