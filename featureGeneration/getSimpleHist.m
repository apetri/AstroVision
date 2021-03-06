function hist = getSimpleHist( img, noise, bound, nbins, sigma )
%Description:
%   first add noise to the image, then using Gaussian filter (with input
%   sigma and 4*sigma size) to smooth the noisy image, finally calculate  
%   the histogram using input bound and nbins.
% 
%Inputs:
%--img: input image
%--noise: given noise map
%--bound: min and max boundary of the bins
%--nbins: number of the bins to calculate histogram
%--sigma: parameter of Gaussian filter
%Outputs:
%--hist: histogram of the image

%--adding noise generated by normrnd function
%img=img+normrnd(0,0.47,size(img));
%--adding given noise map 'noise'
img=img+noise;

width=ceil(4*sigma);
filter=fspecial('gaussian', width, sigma);
fimg=imfilter(img,filter);

bins=linspace(bound(1),bound(2),nbins+1);
binrange=bins(2:end-1);
binrange=binrange(:);
hist=imgHistCount(fimg,binrange);

end
