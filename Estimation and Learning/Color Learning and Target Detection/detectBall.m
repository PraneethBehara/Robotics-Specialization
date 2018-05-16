% Robotics: Estimation and Learning 
% WEEK 1
% 
% Complete this function following the instruction. 

function [segI, loc] = detectBall(I)
% function [segI, loc] = detectBall(I)
%
% INPUT
% I       120x160x3 numerial array 
%
% OUTPUT
% segI    120x160 numeric array
% loc     1x2 or 2x1 numeric array 

    % Returns gaussian probability
    function prob = Gaus(val,mu,sigma)
        prob = (1/sqrt(2*pi).*sigma)*exp(-0.5.*((val-mu)./sigma).^2);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hard code your learned model parameters here

mu = load('gmm_variables.mat', 'mu');
sigma = load('gmm_variables.mat', 'sigma');
z = load('gmm_variables.mat', 'z');

mu = mu.mu;
sigma = sigma.sigma;
w = z.z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find ball-color pixels using your model

%Considering color space RGB
R = zeros(size(I,1),size(I,2));
G = zeros(size(I,1),size(I,2));
B = zeros(size(I,1),size(I,2));

mask = zeros(size(I,1),size(I,2));

R(:,:) = I(:,:,1);
G(:,:) = I(:,:,2);
B(:,:) = I(:,:,3);

R_probs = w(1,1).*Gaus(R,mu(1,1),sigma(1,1)) + w(2,1).*Gaus(R,mu(2,1),sigma(2,1));
G_probs = w(1,2).*Gaus(G,mu(1,2),sigma(1,2)) + w(2,2).*Gaus(G,mu(2,2),sigma(2,2));
B_probs = w(1,3).*Gaus(B,mu(1,3),sigma(1,3)) + w(2,3).*Gaus(B,mu(2,3),sigma(2,3));

thre = 0.4;

for r=1:size(I,1)
    for c=1:size(I,2)
        if (R_probs(r,c)>thre && G_probs(r,c)>thre && B_probs(r,c)>thre)
            mask(r,c) = 1;
        end
    end
end

% figure;
% imshow(mask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do more processing to segment out the right cluster of pixels.
% You may use the following functions.
%   bwconncomp
%   regionprops
% Please see example_bw.m if you need an example code.

mask = imfill(mask,'holes');
se = strel('disk',1);
mask = imerode(mask,se);

% figure;
% imshow(mask);

CC = bwconncomp(mask);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the location of the ball center
% centroid
S = regionprops(CC,'Centroid');

segI = mask;
loc = S(idx).Centroid;
% Note: In this assigment, the center of the segmented ball area will be considered for grading. 
% (You don't need to consider the whole ball shape if the ball is occluded.)

end
