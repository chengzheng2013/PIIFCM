function [img] = neighbour(mu,nu,pi,IMG_size)  
% Return beta 

% IMG_std = std(IMG(:));
% var = IMG_std.^2;
img = zeros(IMG_size);
mu = reshape(mu,IMG_size(1),IMG_size(2));
nu = reshape(nu,IMG_size(1),IMG_size(2));
pi = reshape(pi,IMG_size(1),IMG_size(2));

for i = 1:IMG_size(1)
    for j = 1:IMG_size(2)
        for ii = -1:1
            for jj = -1:1
                i_neig = i+ii;j_neig = j+jj;
                if (i_neig>=1 && i_neig<=IMG_size(1) &&...
					j_neig>=1 && j_neig<=IMG_size(2) &&...
                    (ii~=0 || jj~=0))
                img(i,j) = img(i,j) + ...
                        exp(-((mu(i,j)-mu(i_neig,j_neig)).^2+...
                        (nu(i,j)-nu(i_neig,j_neig)).^2+...
                        (pi(i,j)-pi(i_neig,j_neig)).^2));
                end
            end
        end
    end
end             
                
end