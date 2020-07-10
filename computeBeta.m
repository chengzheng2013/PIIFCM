function [ beta ] = computeBeta( U_idx,IMG_size,gamma )
U_idx = reshape(U_idx,IMG_size);
U_idx = U_idx-1;
w = [1 1 1;1 0 1;1 1 1];
num_of_1 = imfilter(U_idx,w);

beta = U_idx.*(num_of_1*0.25-2);
beta(U_idx == 0) = -0.25*num_of_1(U_idx==0);
beta = exp(gamma*beta);

end

