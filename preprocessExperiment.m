function preprocessExperiment()
% updated 22.10.21 AR

dataRoot = 'D:\data';
expID = '2021-11-16_09_ESPM040';

animalID = expID(15:end);
expRoot = fullfile(dataRoot,animalID,expID);
recordingsRoot = fullfile(expRoot,'recordings');
if ~exist(recordingsRoot)
    mkdir(recordingsRoot);
end

Timeline = load(fullfile(expRoot,[expID,'_Timeline.mat']));
Timeline = Timeline.timelineSession;

%% Process Bonsai stuff
FrameEvents = readtable(fullfile(expRoot,[expID,'_FrameEvents.csv']));
FrameEvents.Properties.VariableNames = {'Frame','Timestamp','Sync','Trial'};

% find BV times when digital flips
flipTimesBV = FrameEvents.Timestamp(find((diff(FrameEvents.Sync))==-1));

% find TL times when digital flips
bvCh = find(ismember(Timeline.chNames,'Bonvision'));
tlDigThresholded = Timeline.daqData(:,bvCh)>2.5;
flipTimesTL = Timeline.time(find((diff(tlDigThresholded))==-1))';

% check NI DAQ caught as many sync pulses as BV produced
pulseDiff = length(flipTimesTL) - length(flipTimesBV);

if pulseDiff > 0
    disp([num2str(pulseDiff),' more pulses in TL']);
    flipTimesTL = flipTimesTL(1:length(flipTimesBV));
elseif pulseDiff < 0
    disp([num2str(pulseDiff*-1),' more pulses in BV']);
    error('Pulse mismatch');
else
    disp('Pulse match');
end

% make model to convert BV time to TL time
% fitlm(input,output); usage flipTimesBVPredicted = predict(mdl1,flipTimesTL);
mdl1 = fitlm(flipTimesBV,flipTimesTL,'linear');

% get trial onset times
% in BV time
trialOnsetTimesBV = [FrameEvents.Timestamp(1);FrameEvents.Timestamp(find(diff(FrameEvents.Trial)==1))];
% in TL time
trialOnsetTimesTL = predict(mdl1,trialOnsetTimesBV);

bvMeta = readtable(fullfile(expRoot,[expID,'_TrialMetadata.csv']));
% % find rows with '/gratings' commands
% stimParamRows = find(strcmp('/gratings',bvMeta.Var4));

% Find rows when stimuli start. This follows the logic that each stimulus
% is preceeded by /end and is terminated with /start (except the first
% which is not preceeded by anything
endRows   = find(strcmp('/end',bvMeta.Var4));
startRows = find(strcmp('/start',bvMeta.Var4));
stimStartRows = [1;endRows(1:end-1)+1];
stimEndRows   = startRows-1;
featsPerStim = stimEndRows - stimStartRows +1;

% do some quality controls here
% 1) number of trials same as frames file?
% 2) expected /experiment etc commands all in place?
% remove other rows
bvMeta = bvMeta(:,8);
% parse params
allTrials = [];
allTrialStruct = [];
allTrialsStr = [];
featureNumber = 0;
for iTrial = 1:size(stimStartRows,1)
    % within each trial go through each feature
    for iFeature = 0:featsPerStim(iTrial)-1
        featureNumber = featureNumber + 1;
        paramVals = string(bvMeta{stimStartRows(iTrial)+iFeature,1});
        % remove all curly and square brackets
        paramVals = erase(paramVals,["{","}","[","]"]);
        paramVals = erase(paramVals,",,");
        paramVals = strsplit(paramVals,',');
        paramVals = paramVals(2:end-1);
        allTrialsStr{featureNumber,1} = num2str(trialOnsetTimesTL(iTrial));
        allTrialsStr{featureNumber,2} = num2str(iTrial);
        allTrials(featureNumber,1) = trialOnsetTimesTL(iTrial);
        allTrials(featureNumber,2) = iTrial;
        for iParam = 1:length(paramVals)
            % +2 because the 1st/2nd columns are time and trial number
            allTrials(featureNumber,iParam+2) = str2num(paramVals(iParam));
            allTrialsStr{featureNumber,iParam+2} = char(paramVals(iParam));
        end
    end
end

% add starttime as parameter
% process the trial meta data to get stimulus parameters
paramNames_grating = {'stimnumber','featurenumber','featuretype','angle','size','x','y','contrast','opacity','phase','freq','speed','dcycle','onset','duration'};

% figure out how many unique stimulus types there are.
% this should be done later by saving the stimulus parameters in the GUI
% when experiment is run, and also saving stimulus order
% 1) turn each complete stim (i.e. all features) into 1 string
catStimParams = [];
allStimParamsCell = [];
for iStim = 1:size(stimStartRows,1)
    stimRows = find(allTrials(:,2)==iStim);
    allStimParamsCell{iStim} = allTrials(stimRows,:);
    allStimParams = allTrials(stimRows,3:end)';
    allStimParams = allStimParams(:)';
    if iStim == 1
        catStimParams = allStimParams;
    else
        catStimParams = padcat(catStimParams,allStimParams);
    end
end

uniqueStims = unique(catStimParams,'rows');

% classify each trial as one of the unique stims
for iTrial = 1:size(stimStartRows,1)
    [tf, index]=ismember(catStimParams(iTrial,:),uniqueStims,'rows');
    allTrialTypes(iTrial,1)=index;
end

% make a matrix for csv output
trialTimeMatrix = [trialOnsetTimesTL,allTrialTypes];

% store the params of each unique stim conditions in a csv
allStimTypes = [];
for iStimType = 1:size(uniqueStims,1)
    featuretype = 0; % grating - add others here later
    [~, firstInstance]=ismember(uniqueStims(iStimType,:),catStimParams,'rows');
    stimParams = allStimParamsCell{firstInstance}(:,3:end);
    stimParams = [ones([size(stimParams,1) 1])*iStimType,[1:size(stimParams,1)]',repmat(featuretype,[[size(stimParams,1) 1]]),stimParams];
    allStimTypes = [allStimTypes;stimParams];
end
    
% add running trace
Encoder = readtable(fullfile(expRoot,[expID,'_Encoder.csv']));
Encoder.Properties.VariableNames = {'Frame','Timestamp','Trial','Position'};
wheelPos = Encoder.Position;
wheelTimestamps = predict(mdl1,Encoder.Timestamp);
% resample wheel to linear timescale
wheelLinearTimescale = wheelTimestamps(1):0.01:wheelTimestamps(end);
wheelPos2 = smooth(interp1(wheelTimestamps,wheelPos,wheelLinearTimescale),50);
wheelSpeed = (([0;diff(wheelPos2)]*-1)*(62/1024))*100;

% save data
bvDataRoot = fullfile(expRoot,'bonsai');
if ~exist(bvDataRoot)
    mkdir(bvDataRoot);
end
writematrix([wheelLinearTimescale',wheelPos2],fullfile(recordingsRoot,'WheelPos.csv'));
writematrix([wheelLinearTimescale',wheelSpeed],fullfile(recordingsRoot,'WheelSpeed.csv'));
writematrix(trialTimeMatrix,fullfile(bvDataRoot,'Trials.csv'));
% info on what each stimulus number contrains - i.e. param values:
writematrix(allStimTypes,fullfile(bvDataRoot,'StimProperties.csv'));
% add more stimulus types here later:
writecell(paramNames_grating,fullfile(bvDataRoot,'FeatureParamNames_0.csv'));

%% process ca2+ imaging traces
% check suite2p folder exists to be processed
if exist(fullfile(expRoot,'suite2p'),'dir')
    doMerge = true;
    
    resampleFreq = 30;
    
    neuropilWeight = 0.7;
    
    alldF = [];
    allF = [];
    allSpikes = [];
    allDepths = [];
    allRoiPix = [];
    allRoiMaps = [];
    
    % this will be used to make all recordings 2 secs shorter than the
    % first ca trace processed to ensure all chs and depths are the same length
    expFrameLength = [];
    
    % extract the tiff header if this hasn't been done already
    %     if ~exist(fullfile(expRoot,'tifHeader.mat'))
    %         tif1 = dir(fullfile(expRoot,'*.tif'));
    %         tif1 = tif1.name;
    %         extractTiffHeader(fullfile(expRoot,tif1));
    %     end
    
    %outputTimes = Timeline.time(1):1/resampleFreq:Timeline.time(end);
    
    % check number of channels
    if exist(fullfile(expRoot,'ch2'))
        % then there are 2 functional channels
        dataPath{1} = fullfile(expRoot,'suite2p');
        dataPath{2} = fullfile(expRoot,'ch2','suite2p');
    else
        dataPath{1} = fullfile(expRoot,'suite2p');
    end
    
    % check number of depths
    depthCount = length(dir(fullfile(dataPath{1},'*plane*')));
    
    if depthCount==1
        % then we might be doing frame averaging
        %load(fullfile(expRoot,'tifHeader.mat'));
        %acqNumAveragedFrames = header.acqNumAveragedFrames;
        acqNumAveragedFrames = 1;
    else
        % then we assume no averaging
        acqNumAveragedFrames = 1;
    end
    
    % determine which channel has frame timing pulses
    neuralFramesIdx  = find(ismember(Timeline.chNames,'MicroscopeFrames'));
    neuralFramesPulses = Timeline.daqData(:,neuralFramesIdx)>1;
    % divide the frame counter by the number of depths & averaging factor.
    %Timeline.rawDAQData(:,neuralFramesIdx)=ceil(Timeline.rawDAQData(:,neuralFramesIdx)/depthCount/acqNumAveragedFrames);
    
    % determine time of each frame
    frameTimes = Timeline.time(diff(neuralFramesPulses)==1);
    framePulsesPerDepth = length(frameTimes)/length(dataPath);
    frameRate = 1/median(diff(frameTimes));
    
    % determine timeline times when we want the Ca signal of each cell
    % +1 and -1 are because we want to make sure we only include frame
    % times which are available at all depths
    outputTimes = frameTimes(1)+1:1/resampleFreq:frameTimes(end)-1;
    
    % for each channel combine all valid rois
    for iCh = 1:length(dataPath)
        alldF{iCh} = [];
        allF{iCh} = [];
        allSpikes{iCh} = [];
        allDepths{iCh} = [];
        allFOV{iCh} = [];
        for iDepth = 0:depthCount-1
            allRoiPix{iCh,iDepth+1} = [];
            allRoiMaps{iCh,iDepth+1} = [];
            % load s2p data
            Fall = load(fullfile(dataPath{iCh},['plane',num2str(iDepth)],'Fall.mat'));
            % check for mismatch between frames trigs and frames in tiff
            if abs(framePulsesPerDepth-size(Fall.F,2))/max([framePulsesPerDepth,size(Fall.F,2)])>0.01
                pcDiff = round(abs(length(frameTimes)-size(Fall.F,2))/max([length(frameTimes),size(Fall.F,2)]) * 100);
                msgbox(['There is a worrying mismatch between between frames trigs and frames in tiff - ',num2str(pcDiff),'% difference']);
                error(['There is a worrying mismatch between between frames trigs and frames in tiff - ',num2str(pcDiff),'% difference']);
            end
            % load numpy file containing cell classification
            cellValid = readNPY(fullfile(dataPath{iCh},['plane',num2str(iDepth)],'iscell.npy'));
            % overall video contamination subtraction
            % ADD BACK THESE 3 LINES
            binpath = fullfile(dataPath{iCh},['\plane',num2str(iDepth),'\data.bin']);
            meanFrameTimecourse = loadSuite2PVideoMeanFrame(binpath,size(Fall.ops.meanImg));
            meanFrameTimecourse = meanFrameTimecourse - min(meanFrameTimecourse);
            % REMOVE THIS LINE
            %meanFrameTimecourse = zeros([1,size(Fall.F,2)]);
            
            % check if any roi Fs are all zero. s2p sometimes throws these
            % up for some reason. if these are found set iscell to false
            zeroROIs = max(Fall.F,[],2)==0 & min(Fall.F,[],2)==0;
            if sum(zeroROIs)>0
                disp(['Warning: ',num2str(sum(zeroROIs)),' zero flat lined rois...']);
                % this was the old way of dealing with them
                % firstValid = find(zeroROIs==0,1);
                % Fall.F(zeroROIs,:)=repmat(Fall.F(firstValid,:),[sum(zeroROIs),1]);
                cellValid(zeroROIs,1) = 0;
                cellValid(zeroROIs,2) = 1;
            end
            
            % remove cells with iscell = 0 but keep record of original
            % suite2p output cell numbers
            Fneu = Fall.Fneu(cellValid(:,1)==1,:);
            Fa = Fall.F(cellValid(:,1)==1,:);
            Spks = Fall.spks(cellValid(:,1)==1,:);
            s2pIndices = find(cellValid(:,1)==1);
            xpix = []; ypix = [];
            validCellIDs = find(cellValid(:,1)==1);
            
            for iCell = 1:length(validCellIDs)
                currentCell = validCellIDs(iCell);
                xpix{end+1} = Fall.stat{currentCell}.xpix;
                ypix{end+1} = Fall.stat{currentCell}.ypix;
            end
            
            % remove potential stimulus artifact - i.e. mean of frame which
            % is extracted above
            Fneu = Fneu - repmat(meanFrameTimecourse,[size(Fneu,1),1]);
            Fa = Fa - repmat(meanFrameTimecourse,[size(Fa,1),1]);
            
            % neuropil subtraction
            F = Fa - (Fneu*neuropilWeight);
            
            % ensure min(corrected F) > 10;
            FMins = min(F,[],2);
            figure;
            subplot(1,2,1);
            hist(FMins);
            title({'Distribution of original','F values of ROIS'});
            if min(FMins) < 20
                disp('Frame mean and neuropil subtraction give ROIs with F < 20')
                disp(['Offsetting all F by ',num2str((min(FMins)*-1)+20)]);
                F = F + (min(FMins)*-1)+20;
            end
            FMins = min(F,[],2);
            subplot(1,2,2);
            hist(FMins);
            title({'Distribution of F values', 'of ROIS after forcing > 20'});
            drawnow
            
            % decide max experiment length in frames
            %             if isempty(expFrameLength)
            %                 expFrameLength = size(F,2)-20;
            %                 disp('CHECK THIS LINE');
            %                 frameTimes = frameTimes(1:expFrameLength*depthCount);
            %             end
            
            if doMerge
                % merge ROIS with > corrThreshold correlation
                % smoothing window for smoothing before calculating correlations
                smoothWindow = 5;
                Fsmoothed = conv2(F,ones(1,ceil(smoothWindow/frameRate))/ceil(smoothWindow/frameRate),'valid');
                corrThreshold = 0.8;
                corrMatrix = corr(Fsmoothed');
                % figure;imagesc(corrMatrix);
                % Bill's merge code
                in_cluster = []; %vector to store cells already in cluster
                clusters = {}; %cell array. Each cell is a cluster, and each element in the cell is a ROI.
                n_clusters = 0;
                n_rois = size(corrMatrix,1);
                for c = 1:n_rois
                    if ~ismember(c, in_cluster)
                        n_clusters = n_clusters + 1;
                        clusters{n_clusters} = [];
                        for c2 = c:n_rois
                            if ~ismember(c2, in_cluster)
                                if corrMatrix(c,c2) > corrThreshold
                                    clusters{n_clusters} = horzcat(clusters{n_clusters}, c2);
                                    in_cluster = [in_cluster, c2];
                                end
                            end
                        end
                    end
                end
                % plot merges
                % for iCluster = 1:20
                %     figure
                %     plot(Fsmoothed(clusters{iCluster},:)');
                %     pause(1);
                % end
                
                % do weighted averaging of rois in clusters and collect
                % together all pixels
                weightedMerge = zeros(length(clusters),size(F,2));
                weightedMergeSpks = zeros(length(clusters),size(F,2));
                % make a blank roi map
                roiMap = zeros(size(Fall.ops.meanImg));
                for iCluster = 1:length(clusters)
                    % pull out total number of pix across all rois
                    pixMerged{iCluster} = [];
                    totalPix = 0;
                    for iRoi = 1:length(clusters{iCluster})
                        roiID = clusters{iCluster}(iRoi);
                        totalPix = totalPix + length(xpix{roiID});
                        iOfPix = sub2ind(size(Fall.ops.meanImg),ypix{roiID}+1,xpix{roiID}+1);
                        pixMerged{iCluster} = [pixMerged{iCluster},iOfPix];
                    end
                    % make a weighted total adding 1 roi at a time
                    for iRoi = 1:length(clusters{iCluster})
                        roiWeight = length(xpix{clusters{iCluster}(iRoi)})/totalPix;
                        weightedMerge(iCluster,:) = weightedMerge(iCluster,:) + (F(clusters{iCluster}(iRoi),:)*roiWeight);
                        weightedMergeSpks(iCluster,:) = weightedMergeSpks(iCluster,:) + (Spks(clusters{iCluster}(iRoi),:)*roiWeight);
                    end
                end
                
                % make a roi map for the depth that can be used for longitudinal imaging etc
                for iCluster = 1:length(clusters)
                    % remove duplicate pix in the ROI
                    pixMerged{iCluster} = unique(pixMerged{iCluster});
                    % label ROI map
                    roiMap(pixMerged{iCluster}) = iCluster;
                end
                
                % Update F to merged F
                F = weightedMerge;
                Spks = weightedMergeSpks;
                roiPix = pixMerged;
                disp(['Merged ',num2str(size(F,1)),' rois --> ',num2str(length(clusters))]);
            else
                % what to do if not merging
                % make a roi map for the depth that can be used for longitudinal imaging etc
                roiPix = [];
                for iRoi = 1:size(F,1)
                    % collect pix in ROI
                    roiPix{iRoi} = sub2ind(size(Fall.ops.meanImg),ypix{iRoi}+1,xpix{iRoi}+1);
                    % label ROI map
                    roiMap(roiPix{iRoi}) = iCluster;
                end
            end
            
            % crop F down to above established max frames
            %F = F(:,1:expFrameLength);
            
            % dF/F calculation
            smoothingWindowSize = 100;
            smoothed = conv2(F,ones(1,smoothingWindowSize)/smoothingWindowSize,'same');
            % remove edge effects
            smoothed(:,1:smoothingWindowSize)=repmat(smoothed(:,smoothingWindowSize+1),[1,smoothingWindowSize]);
            smoothed(:,end-smoothingWindowSize+1:end)=repmat(smoothed(:,end-smoothingWindowSize-1),[1,smoothingWindowSize]);
            % replace nans with large values (so they don't get picked up as mins)
            smoothed(isnan(smoothed))=max(smoothed(:))*2;
            baseline = imerode(smoothed,strel('rectangle',[1 smoothingWindowSize]));
            % calculate dF/F
            dF = (F-baseline)./baseline;
            % get times of each frame
            depthFrameTimes = frameTimes(iDepth+1:depthCount:length(frameTimes));
            depthFrameTimes = depthFrameTimes(1:size(dF,2));
            % resample to get desired sampling rate
            dF = interp1(depthFrameTimes,dF',outputTimes)';
            F = interp1(depthFrameTimes,F',outputTimes)';
            Spks = interp1(depthFrameTimes,Spks',outputTimes)';
            if size(dF,2)==1
                dF = dF';
            end
            % pick out valid cells
            alldF{iCh} = [alldF{iCh};dF];
            allF{iCh} = [allF{iCh};F];
            allSpikes{iCh} = [allSpikes{iCh};Spks];
            
            allDepths{iCh} = [allDepths{iCh};repmat(iDepth,[sum(cellValid(:,1)),1])];
            allRoiPix{iCh,iDepth+1} = roiPix;
            allRoiMaps{iCh,iDepth+1} = roiMap;
            
            allFOV{iCh} = Fall.ops.meanImg;
            
        end
    end
    
    disp('Saving 2-photon data...');
    % save as CSV
    for iCh = 1:length(alldF)
        writematrix([outputTimes;alldF{iCh}]',fullfile(recordingsRoot,['dF_',num2str(iCh),'.csv']));
        writematrix([outputTimes;allF{iCh}]',fullfile(recordingsRoot,['F_',num2str(iCh),'.csv']));
        writematrix([outputTimes;allSpikes{iCh}]',fullfile(recordingsRoot,['Spikes_',num2str(iCh),'.csv']));
        writematrix([allRoiMaps{iCh}]',fullfile(recordingsRoot,['roi_',num2str(iCh),'.csv']));
        writematrix([allFOV{iCh}]',fullfile(recordingsRoot,['fov_',num2str(iCh),'.csv']));
    end
    
    % save for matlab
    s2pData.alldF = alldF;
    s2pData.allF = allF;
    s2pData.allDepths = allDepths;
    s2pData.allRoiPix = allRoiPix;
    s2pData.allRoiMaps = allRoiMaps;
    s2pData.meanFrame = Fall.ops.meanImg;
    s2pData.t = outputTimes;
    save(fullfile(recordingsRoot,'s2pData.mat'),'s2pData');
    
    % save in python format eventually...
    
end

%% process ePhys data
ePhys1Idx  = find(ismember(Timeline.chNames,'EPhys1'));
ePhys2Idx  = find(ismember(Timeline.chNames,'EPhys2'));
writematrix([Timeline.time',Timeline.daqData(:,[ePhys1Idx ePhys2Idx])],fullfile(recordingsRoot,'ephys.csv'));

%% process camera pulses
try
    load(fullfile(expRoot,[expID,'_eyeMeta1.mat']));
    camIdx  = find(ismember(Timeline.chNames,'EyeCamera'));
    camPulseTrace = Timeline.daqData(:,camIdx)>3;
    framePulseTimes = Timeline.time(find(abs(diff(camPulseTrace))==1));
    framePulseFrameNumbers = 1:100:(length(framePulseTimes))*100;
    %frameSampleTimes = interp1(framePulseTimes,framePulseFrameNumbers,
    frameSampleTimes = framePulseTimes(1):1/20:framePulseTimes(end)+median(diff(framePulseTimes));
    frameNumber = interp1(framePulseTimes',framePulseFrameNumbers',framePulseTimes(1):1/20:framePulseTimes(end)+(median(diff(framePulseTimes))*4)','linear','extrap');
    frameNumber = floor(frameNumber);
    validSamples = 1:find(frameNumber==eTrackData.frameCount,1,'last');
    frameNumber = frameNumber(validSamples);
    frameSampleTimes = frameSampleTimes(validSamples);
    frameRate = 100 / median(diff(framePulseTimes));
    disp(['Detected eye cam frame rate = ',num2str(frameRate),'Hz']);
    writematrix([frameSampleTimes',frameNumber'],fullfile(recordingsRoot,'eyeFrames.csv'));
catch
    disp('Camera data NOT processed');
end
% CSF files to create (time always in rows):
% Trial properties numeric, rows are trials (inc onset timestamps)
% Trial properties text, rows are trials (inc onset timestamps)
% Stim parameter names (first column = timestamps)
% 
% wheel position / timestamps
%
% 2P frame timestamps
% 2P ca traces
% 
% Cam frame timestamps
% Cam frame numbers
% 
% ePhys timestamps
% ePhys trace 1
% ePhys trace 2

% save all CSV files

% %% DEBUGGING 
% timeDiff = (flipTimesTL(1:end)) - (flipTimesBV(1:end));
% %timeDiff = timeDiff - mean(timeDiff);
% 
% figure;
% plot(flipTimesTL,timeDiff);
% xlabel("time elapsed in experiment according to NI DAQ (secs)")
% ylabel("diff of time of pulse edge detection in NI DAQ vs. frames file")
% 
% figure;
% diff2 = diff(timeDiff);
% histogram(diff2(diff2>.01),0:0.0167/4:0.0167*5);
% xlabel("Blip amplitude (how far out of time the 2 systems are)")
% ylabel("Occurance frequency")
% 
% figure;
% subplot(1,2,1);
% histogram(diff(flipTimesBV));
% subplot(1,2,2);
% histogram(diff(flipTimesTL));

disp('All done');
end

function frameMean = loadSuite2PVideoMeanFrame(pathToBinary,frameSize)

% fclose all
% pathToBinary = 'D:\ACC Data\2016-05-28_02_CFEB014\suite2p\plane0\data.bin';
% frameSize = [256 256];
blockSize = 1000;

finfo = dir(pathToBinary);
fsize = finfo.bytes;
fid = fopen(pathToBinary);
frameCountCalculation = fsize/frameSize(1)/frameSize(2)/2;

framesInSet = int16([256,256,frameCountCalculation]);

for iStart = 1:blockSize:frameCountCalculation
    disp(['Frame ',num2str(iStart),' of ',num2str(frameCountCalculation)]);
    lastFrame = iStart + blockSize-1;
    lastFrame = min(lastFrame,frameCountCalculation);
    framesToRead = lastFrame - iStart + 1;
    framesInSet(1:frameSize(1),1:frameSize(2),iStart:lastFrame) = reshape(fread(fid,frameSize(1)*frameSize(2)*framesToRead,'int16'),[frameSize(1),frameSize(2),framesToRead]);
end

framesInSet = mean(framesInSet,1);
framesInSet = mean(framesInSet,2);
frameMean = squeeze(framesInSet);

fclose(fid);

end

function M = padcat(original,addition)
% concatinates rows and expand original with nans if needed to accomodate addition.
% check if original or addition has more columns
oSize = size(original,2);
aSize = size(addition,2);
if oSize>aSize
    colsToAdd = oSize - aSize;
    addition = [addition,zeros([1 colsToAdd])];
elseif oSize<aSize
    colsToAdd = aSize - oSize;
    original = [original,zeros([1 colsToAdd])];
end
M = [original;addition];    
end

function data = readNPY(filename)
% Function to read NPY files into matlab.
% *** Only reads a subset of all possible NPY files, specifically N-D arrays of certain data types.
% See https://github.com/kwikteam/npy-matlab/blob/master/tests/npy.ipynb for
% more.
%

[shape, dataType, fortranOrder, littleEndian, totalHeaderLength, ~] = readNPYheader(filename);

if littleEndian
    fid = fopen(filename, 'r', 'l');
else
    fid = fopen(filename, 'r', 'b');
end

try

    [~] = fread(fid, totalHeaderLength, 'uint8');

    % read the data
    data = fread(fid, prod(shape), [dataType '=>' dataType]);

    if length(shape)>1 && ~fortranOrder
        data = reshape(data, shape(end:-1:1));
        data = permute(data, [length(shape):-1:1]);
    elseif length(shape)>1
        data = reshape(data, shape);
    end

    fclose(fid);

catch me
    fclose(fid);
    rethrow(me);
end
end

