function convert_stim_for_py_batch(Folder)
if ~exist('Folder','var')
  Folder = 'C:\Local_Repository\ESMT072';
end
% Folder = 'G:\.shortcut-targets-by-id\1P7g8LSE5D6vInT7OOXY1EIzJ0M4zvhos\Remote_Repository\TEST\2023-02-27_09_TEST';
FileList = dir(fullfile(Folder, '**', '*stim.mat'));
disp(['Found ',num2str(length(FileList)),' stim files']);
for iFile = 1:length(FileList)
    % load matlab stim file
    load(fullfile(FileList(iFile).folder,FileList(iFile).name));
    disp(['Starting ',expDat.expID,'(',num2str(iFile),'/',num2str(length(FileList)),')...']);
    % determine max number of feats
    max_feat = 0;
    for iStim = 1:length(expDat.stims)
        if length(expDat.stims(iStim).features)>max_feat
            max_feat = length(expDat.stims(iStim).features);
        end
    end
    % determine the total stimulus duration by finding the max start +
    % duration value
    allDurations = zeros([length(expDat.stims),1]);
    for iStim = 1:length(expDat.stims)
        for iFeat = 1:length(expDat.stims(iStim).features)
            [~,startIdx] = ismember('onset',expDat.stims(iStim).features(iFeat).params);
            [~,durationIdx] = ismember('duration',expDat.stims(iStim).features(iFeat).params);
            if str2num(expDat.stims(iStim).features(iFeat).vals{startIdx})+str2num(expDat.stims(iStim).features(iFeat).vals{durationIdx}) > allDurations(iStim)
                allDurations(iStim) =  str2num(expDat.stims(iStim).features(iFeat).vals{startIdx})+str2num(expDat.stims(iStim).features(iFeat).vals{durationIdx});
            end
        end
    end
    % make empty variable for param names
    for iFeat = 1:max_feat
        all_param_names{iFeat} = {};
    end
    % for each feature number discover all unique parameter names
    for iStim = 1:length(expDat.stims)
        for iFeat = 1:length(expDat.stims(iStim).features)
            all_param_names{iFeat} = cat(2,all_param_names{iFeat},expDat.stims(iStim).features(iFeat).params,'type');
        end
    end
    % find unique param names for each feature
    for iFeat = 1:max_feat
        all_param_names{iFeat} = unique(all_param_names{iFeat});
    end
    % build a table for csv output for each feature
    output_table = [];
    for iStim = 1:length(expDat.stims)
        for iFeat = 1:max_feat
            % check stim has enough features to probe iFeat
            if length(expDat.stims(iStim).features)>=iFeat
                for iParam = 1:length(all_param_names{iFeat})
                    % find if the param name exists in the stim/feature
                    [Lia,Locb] = ismember(all_param_names{iFeat}{iParam},expDat.stims(iStim).features(iFeat).params);
                    if Lia
                        output_table{iFeat}(iStim,iParam) = {expDat.stims(iStim).features(iFeat).vals{Locb}};
                    else
                        if strcmp(all_param_names{iFeat}{iParam},'type')
                            output_table{iFeat}(iStim,iParam) = {expDat.stims(iStim).features(iFeat).name};
                        else
                            output_table{iFeat}(iStim,iParam) = {'NaN'};
                        end
                    end
                end
            else
                % pad with nans
                output_table{iFeat}(iStim,1:length(all_param_names{iFeat})) = {'NaN'};
            end
        end
    end
    % combine the feature tables into 1 table
    combined_output_table = [];
    for iFeat = 1:max_feat
        combined_output_table = [combined_output_table,output_table{iFeat}];
    end
    % make a table where each row is a trial
    combined_output_table_all_trials = [];
    for iStim = 1:length(expDat.stimOrder)
        combined_output_table_all_trials = cat(1,combined_output_table_all_trials,combined_output_table(expDat.stimOrder(iStim),:));
    end
    % rename params to append feat number in headers
    param_header = [];
    for iFeat = 1:max_feat
        all_feat_headers = [];
        for iParam = 1:length(all_param_names{iFeat})
            all_feat_headers{iParam} = ['F',num2str(iFeat),'_',all_param_names{iFeat}{iParam}];
        end
        param_header = [param_header,all_feat_headers];
    end
    % add column for stim max length
    cellstr(num2str(allDurations))
    param_table = cell2table(cat(2,cellstr(num2str(allDurations)),combined_output_table), 'VariableNames',cat(2,{'duration'},param_header));
    all_trials_table = cell2table(cat(2,cellstr(num2str(expDat.stimOrder')),cellstr(num2str(allDurations(expDat.stimOrder))),combined_output_table_all_trials),'VariableNames',cat(2,{'stim'},{'duration'},param_header));
    % save as csv
    save_path = FileList(iFile).folder;
    expID = expDat.expID;
    writematrix(expDat.stimOrder',fullfile(save_path,[expID,'_stim_order.csv']));
    writetable(param_table,fullfile(save_path,[expID,'_stim.csv']));
    writetable(all_trials_table,fullfile(save_path,[expID,'_all_trials.csv']));
    disp(['Finished ',expDat.expID]);
end
end
