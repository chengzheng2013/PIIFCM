clear ;close all ;clc

str1 = 'input_dir\';
str2 = 'output_dir\';

% initlize parameters
lambda = 2.5;
m =2;
gamma = 3.5;
tic
for flag =1:200
    fprintf('iteration %d/200...\n',flag);
    IMG = imread(strcat(str1,num2str(flag),'_bw.bmp'));
    IMG = IMG(:,:,1);
    IMG = im2double(IMG);
    IMG = PIIFCM(IMG,m,lambda,gamma);
    imwrite(IMG,strcat(str2,num2str(flag),'_bw.bmp'));
end
toc