function [mu, sigma, z] = gmm_mle()
%%%
% This function returns mu,sigma,val using EM(Expectation-Maximization method
% from samples collected
%%%
    
    % Returns gaussian probability
    function prob = Gaus(val,mu,sigma)
        prob = (1/sqrt(2*pi).*sigma)*exp(-0.5.*((val-mu)./sigma).^2);
    end

mu_i = [126,121,66; 146,136,62];
sigma_i = [8,20,10; 4,15,8];
w_i = [0.5,0.45,0.6; 0.5,0.55,0.4];
K = 2;
D = 3;

% Initializing variables
Samples = load('Samples.mat');
data = zeros(size(Samples.Samples));
zi = zeros(size(data,1),size(data,2),K);
zk = zeros(K,D);

mu = mu_i;
sigma = sigma_i;
z = w_i;

data(:,1) = Samples.Samples(:,1);
data(:,2) = Samples.Samples(:,2);
data(:,3) = Samples.Samples(:,3);

% Display histograms
% figure; histogram(samples_data(:,1),'EdgeColor','black');
% figure; histogram(samples_data(:,2),'EdgeColor','black');
% figure; histogram(samples_data(:,3),'EdgeColor','black');

for i=1:100
    % E-step
    for k=1:K
        for d=1:D
            z1 = z(1,d).*Gaus(data(:,d),mu(1,d),sigma(1,d)); 
            z2 = z(2,d).*Gaus(data(:,d),mu(2,d),sigma(2,d));
            zi(:,d,k) = (z(k,d).*Gaus(data(:,d),mu(k,d),sigma(k,d)))./(z1+z2);
            
            %M-step
            zk(k,d) = sum(zi(:,d,k));
            mu_updated(k,d) = sum(data(:,d).*zi(:,d,k))./(zk(k,d));
            sig_updated(k,d)= sqrt(sum(zi(:,d,k).*((data(:,d) - mu_updated(k,d)).^2))./zk(k,d));
            z_updated(k,d) = zk(k,d)./size(data,1);
        end
    end
    mu = mu_updated;
    sigma = sig_updated;
    z = z_updated;
end

save('gmm_variables','mu','sigma','z');
end
