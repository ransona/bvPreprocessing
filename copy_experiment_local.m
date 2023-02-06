function copy_experiment_local(expID,local_repos,remote_repos)
animalID = expID(15:end);
expRootLocal = fullfile(local_repos,animalID,expID);
expRootRemote = fullfile(remote_repos,animalID,expID);
time_start = tic;
% check directory exists
if ~exist(expRootRemote)
    error('Remote directory not found');
    return;
end
if ~exist(expRootLocal)
    mkdir(expRootLocal);
end
% copy everything except the Trials Folder and tif files
allfiles = dir(expRootRemote);
excluded_folders = {'suite2p','Trials','.','..'};
excluded_extensions = {'.tif'};
for iFile = 1:length(allfiles)
    [~,~,ext] = fileparts(allfiles(iFile).name);
    if allfiles(iFile).isdir
        if ~ismember(allfiles(iFile).name,excluded_folders)
            copyfile(fullfile(expRootRemote,allfiles(iFile).name),fullfile(expRootLocal,allfiles(iFile).name));
            disp(['Copying ',allfiles(iFile).name]);
        end
    else
        if ~ismember(ext,excluded_extensions)
            copyfile(fullfile(expRootRemote,allfiles(iFile).name),fullfile(expRootLocal,allfiles(iFile).name));
            disp(['Copying ',allfiles(iFile).name]);
        end  
    end        
end
disp(['Copy completed without errors in ',num2str(toc(time_start))]);
end