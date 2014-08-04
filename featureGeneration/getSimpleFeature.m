function getSimpleFeature( datasets, bound, nbins, sigma, savefile )
%Compute one level features for given datasets and parameter settings
%Inputs:
%--datasets: a cell array of the name of the datasets, e.g. {'ref','om29'}
%--bound: an array of the left and right boundary of the bins, e.g.[-0.15,0.15]
%--nbins: number of bins
%--sigma: the sigma value of the Gaussian filter
%--savefile: the suffix (before '.mat') of the file where the generated 
%features are stored

cntn=0;
%This is the folder and part of the name of the files where noise maps are stored
noisename='noise/noise_map_';
for dataset=datasets
    folder=['datasets/' dataset{1}];
    fileList=getFilesFolder(folder);
    cnt=0;
    fprintf('%s\n', folder);
    for file=fileList
        cntn=cntn+1;
        cnt=cnt+1;
        load([noisename num2str(cnt) '.mat']);
        fname=file{1};
        fprintf('%d ... ', cnt);
        try
            % might throw FITS read error
            % mat2gray to convert to [0, 1]
            img = fitsread([folder '/' fname]);
        catch err
            fprintf(['ERROR:\t' fname '\n']);
            cntn=cntn-1;
            cnt=cnt-1;
            continue;
        end 
        hist=getSimpleHist(img, noise_map, bound, nbins, sigma);
        feature(cnt,1:numel(hist))=hist;
    end
    feature=bsxfun(@rdivide,feature,sum(feature,2));
    saveVar(['Feature_' num2str(nbins) 'Bins_' savefile], dataset{1}, feature);
end

end
