function simpleChi(filename,savefile)
%Calculate the DeltaChi2 of simple features, without PCA
%Inputs:
%--filename: the name of the file where the features are stored
%--savefile: the name of the file to save the results

load(filename);

datasets={'ref','om29','si850','w12'};

for i=1:length(datasets)
    [sigma{i},miu{i}]=computeParas(eval(datasets{i}));
end

chi1=fisherFunc(sigma{1},miu{1},miu{2});
chi2=fisherFunc(sigma{1},miu{1},miu{3});
chi3=fisherFunc(sigma{1},miu{1},miu{4});

saveVar(savefile,'chi1',chi1);
saveVar(savefile,'chi2',chi2);
saveVar(savefile,'chi3',chi3);

end
