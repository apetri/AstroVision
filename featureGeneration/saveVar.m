function saveVar( fname, vname, data )
%save variable into selected file
%Inputs:
%--fname: the name of the file to store the variable
%--vname: the name with which to save the variable
%--data: the variable to save

if ~strcmp(fname(end-3:end),'.mat')
    fname=[fname '.mat'];
end

eval([vname '=data;']);

if exist(fname,'file')
    save(fname,vname,'-append');
else
    save(fname,vname);
end

end
