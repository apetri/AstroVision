%% Customized version of histogram counting
% binrange: obtained by cluster1D()
%
function counts = imgHistCount(img, binrange)

counts = histc(img(:), [-Inf; binrange; Inf]);
counts = counts(1:end-1);
