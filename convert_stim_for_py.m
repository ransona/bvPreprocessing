% determine max number of feats
max_feat = 0;
for iStim = 1:length(expDat.stims)
    if length(expDat.stims(iStim).features)>max_feat
        max_feat = length(expDat.stims(iStim).features);
    end
end
% make empty variable for param names
for iFeat = 1:max_feat
    all_param_names{iFeat} = {};
end
% for each feature number discover all unique parameter names
for iStim = 1:length(expDat.stims)
    for iFeat = 1:length(expDat.stims(iStim).features)
        all_param_names{iFeat} = cat(2,all_param_names{iFeat},expDat.stims(iStim).features(iFeat).params);
    end
end
% find unique param names for each feature
for iFeat = 1:max_feat
    all_param_names{iFeat} = unique(all_param_names{iFeat});
end
% build a table for csv output for each feature
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
                    output_table{iFeat}(iStim,iParam) = {'NaN'};
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
% rename params to append feat number in headers
param_header = [];
for iFeat = 1:max_feat
    all_feat_headers = [];
    for iParam = 1:length(all_param_names{iFeat})
        all_feat_headers{iParam} = ['F',num2str(iFeat),'_',all_param_names{iFeat}{iParam}];
    end
    param_header = [param_header,all_feat_headers];
end

param_table = cell2table(combined_output_table, 'VariableNames',param_header);
% save as csv
writetable(param_table,'G:\.shortcut-targets-by-id\18E8Ww5qCgzn27qk_LrsR3R6wMweRdJ_Z\AR_RRP\ESMT101\2022-10-19_08_ESMT101\stim.csv')

