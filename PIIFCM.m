function X_recovered = PIIFCM(IMG,m,lambda,gamma)

IMG_size = size(IMG);

IMG1 = IMG;
% IMG = t(IMG1);

[p,sigma] = possibility(IMG,3e-6,10);

if ((sigma(1,1)/sigma(2,2))<0.5 && abs(sigma(1,2))<10)
    [p,~] = possibility(IMG1,1e-6,0);
    IMG = IMG1;
elseif ((sigma(1,1)/sigma(2,2))<1 && abs(sigma(1,2))<100)
    [p,~] = possibility(IMG1,1e-6,1);
    IMG = IMG1;
end

X = reshape(IMG,IMG_size(1)*IMG_size(2),1);%X is N-by-1 
[N,~] = size(X);

%select value for c,m,eta,epsilon,max_iters,a,b,lambda,K
c = 2;
epsilon = 1e-6;
max_iters = 100;

%initalize the fuzzy partition matrix U(0),U is N-by-c
U = rand(N,c);
U = U./repmat(sum(U,2),1,c);

%initalize the membership degree matrix mu(0),mu is N-by-1
X_max = max(X);
X_min = min(X);
mu = (X-X_min)./(X_max-X_min);

%initalize the non-membership degree matrix nu(0),nu is N-by-1
nu = (1-mu)./(1+lambda*mu);

%initalize the hesitation degree matrix pi(0),pi is N-by-1
pi = 1-mu-nu;

w1 = [0.4142,0.5,0.4142 ; 0.5,0,0.5 ; 0.4142,0.5,0.4142];
w1 = w1/sum(sum(w1));

%reshape mu,nu,pi to size(IMG)
mu_res = reshape(mu,IMG_size(1),IMG_size(2));
nu_res = reshape(nu,IMG_size(1),IMG_size(2));
pi_res = reshape(pi,IMG_size(1),IMG_size(2));

%initalize W, R, intuitionistic fuzzy distance, gamma, w, centroids,Nr
W = ones(N,c);
distance = zeros(N,c);
centroids_mu = zeros(c,1);
centroids_nu = zeros(c,1);
centroids_pi = zeros(c,1);
% J = zeros(max_iters,1);
al_be_distance = zeros(N,c);
for iter = 1:max_iters
%     W = ones(N,c);
    [~,ii] = sort(centroids_mu);
    W(:,ii(1)) = p.^(2/3);
    W(:,ii(2)) = p.^(1/6);   
%     W = ones(N,c);
    
    %compute beta for each pixel 
    [~,U_idx] = max(U,[],2);
    beta = computeBeta(U_idx,IMG_size,gamma);
    beta = beta.*reshape(p,IMG_size(1),IMG_size(2));
    %compute alpha*beta;
    alpha_beta = imfilter(beta,w1);
    al_be = reshape(alpha_beta,IMG_size(1)*IMG_size(2),1);
    %compute the alp_bet_mu,al_be_nu,and al_be_pi
    al_be_mu = imfilter(beta.*mu_res,w1);
    al_be_nu = imfilter(beta.*nu_res,w1);
    al_be_pi = imfilter(beta.*pi_res,w1);
    al_be_mu = reshape(al_be_mu,IMG_size(1)*IMG_size(2),1);
    al_be_nu = reshape(al_be_nu,IMG_size(1)*IMG_size(2),1);
    al_be_pi = reshape(al_be_pi,IMG_size(1)*IMG_size(2),1);

    %calculate the centroids_mu ,centroids_nu and centroids_pi
    centroids_mu = ((W.*(U.^m))'*(mu+al_be_mu))./...
                    (sum((W.*(U.^m).*(1+repmat(al_be,1,c))),1))';
    centroids_nu = ((W.*(U.^m))'*(nu+al_be_nu))./...
                    (sum((W.*(U.^m).*(1+repmat(al_be,1,c))),1))';
    centroids_pi = ((W.*(U.^m))'*(pi+al_be_pi))./...
                    (sum((W.*(U.^m).*(1+repmat(al_be,1,c))),1))';
    
    
    %calculate the intuitionistic fuzzy distance
    distance = (repmat(mu,1,c)-repmat(centroids_mu',N,1)).^2+...
               (repmat(nu,1,c)-repmat(centroids_nu',N,1)).^2+...
               (repmat(pi,1,c)-repmat(centroids_pi',N,1)).^2;
           
    %compute al_be_distance    
    for cen = 1:c
        dis_res = reshape(distance(:,cen),IMG_size(1),IMG_size(2));
        al_be_distance_cen = imfilter(dis_res.*beta,w1);
        al_be_distance(:,cen) = reshape(al_be_distance_cen,IMG_size(1)*IMG_size(2),1);
    end
    
    [~,ii] = sort(centroids_mu);
    W(:,ii(1)) = p.^(2/3);
    W(:,ii(2)) = p.^(1/6);
%     W = ones(N,c);
    %compute the membership matrix U
    U_new = ((W.*(distance+al_be_distance)).^(-1/(m-1)))./...
            repmat(sum(((W.*(distance+al_be_distance)).^(-1/(m-1))),2),1,c);

    %check out the termination condition
    if (max(max(abs(U_new-U))) < epsilon)
        break;
    else
        U = U_new;
    end

end

%recover the image from idx by map each pixel to the centroid value
[~,idx] = max(U,[],2);
[~,I] = sort(centroids_mu);
centroids_mu_sort = centroids_mu;
for i = 1:c
    centroids_mu_sort(I(i)) = i;
end
X_recovered = centroids_mu_sort(idx,:);
X_recovered = reshape(X_recovered,IMG_size(1),IMG_size(2));
X_recovered(X_recovered==1) = 0;
X_recovered(X_recovered==2) = 255;
X_recovered = uint8(X_recovered);

end

