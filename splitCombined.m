
function splitCombined()

% how it works:
% 1) take a bunch of recordings and add them to suite2p to be processed by
%    clicking 'add directory to data path'
% 2) suite2p will then collect all tiffs together and treat them as one
%    experiment
% 3) this 'combined' experiment will be processed and saved to the folder of
%    the first experiment you added
% 4) next you need to do the usual validations in suite2p:
%    a) make sure the rois that are on the left are all sensible
%    b) go to registration -> view registration metrics and make sure the
%       registration has worked well. tell me if you don't know how to
%       interpret what you see there. you can also watch this video from
%       Carsen: https://cbmm.mit.edu/video/suite2p-fast-and-accurate-pipeline-automatically-processing-functional-imaging-recordings
%       Check out at around 33-35 mins. Although the whole thing is worth
%       watching.
% 4) next open this function (splitCombined) and edit the data_root to
%    where your data is, and the first_exp_id to be the ID of the first
%    experiment you added to suite 2p.
% 5) it will transfer the relevant part of the suite2p output from each
%    experiment to each experiment's folder. you should then be able to run
%    preprocessExperiment() as usual

% Information to edit:
% data root
data_root = 'V:\Local_Repository';
first_exp_id = '2022-02-07_01_ESPM039';

% begin processing
animal_id = first_exp_id(15:end);
% the path to the first experiment in the list
exp_path = fullfile(data_root,animal_id,first_exp_id);
% check if the experiment has already been partially processed using the
% script in which case there will be a 'suite2p_combined' folder
if ~exist(fullfile(exp_path,'suite2p_combined'),'dir')
    % rename the suite2p folder in the first experiment's folder
    movefile(fullfile(exp_path,'suite2p'),fullfile(exp_path,'suite2p_combined'));
end
% set the renamed folder as the place to find the suite2p data
exp_path = fullfile(exp_path,'suite2p_combined');
% determine number of planes
planes_list = dir(fullfile(exp_path,'*plane*'));
% determine all experiment IDs that have been combined
Fall_original = load(fullfile(exp_path,['plane',num2str(0)],'Fall.mat'));
for i = 1:size(Fall_original.ops.data_path,1)
    [~,expIDs{i},~]=fileparts(Fall_original.ops.data_path(i,:));
end

% verification checks
% 1) all from same animal?
for iExp = 1:length(expIDs)
    disp(['Checking ',expIDs{iExp},'...']);
    all_animal_ids{iExp} = expIDs{iExp}(15:end);
    if iExp>1
        if strcmp(all_animal_ids{iExp},all_animal_ids{1})==0
            error('You are combining data from experiments from different animals')
        end
    end
end

disp('Verification checks passed');

all_mean_frames = [];

for iPlane = 0:length(planes_list)-1
    disp(['=====Plane ',num2str(iPlane),'=====']);
    % load data for plane
    Fall_original = load(fullfile(exp_path,['plane',num2str(iPlane)],'Fall.mat'));
    % load npy data
    npy_original.F = readNPY(fullfile(exp_path,['plane',num2str(iPlane)],'F.npy'));
    npy_original.Fneu = readNPY(fullfile(exp_path,['plane',num2str(iPlane)],'Fneu.npy'));
    npy_original.spks = readNPY(fullfile(exp_path,['plane',num2str(iPlane)],'spks.npy'));
    % split it up into seperate files
    for iExp = 1:length(expIDs)
        disp(['Processing experiment ',expIDs{iExp}]);
        % pull out frames that were extracted from tiffs in that experiment's folder
        frames_in_exp = Fall_original.ops.frames_per_folder(iExp);
        % calculate which frame in the combined data is the first from this experiment
        exp_Start_frame = sum(Fall_original.ops.frames_per_folder(1:iExp-1))+1;
        Fall = Fall_original;
        % select only frames from this experiment
        Fall.F = Fall.F(:,exp_Start_frame:exp_Start_frame+frames_in_exp-1);
        Fall.Fneu = Fall.Fneu(:,exp_Start_frame:exp_Start_frame+frames_in_exp-1);
        Fall.spks = Fall.spks(:,exp_Start_frame:exp_Start_frame+frames_in_exp-1);
        % do the same for the npy data
        npy = npy_original;
        npy.F = npy.F(:,exp_Start_frame:exp_Start_frame+frames_in_exp-1);
        npy.Fneu = npy.Fneu(:,exp_Start_frame:exp_Start_frame+frames_in_exp-1);
        npy.spks = npy.spks(:,exp_Start_frame:exp_Start_frame+frames_in_exp-1);
        % save data
        output_folder = fullfile(data_root,animal_id,expIDs{iExp},'suite2p',['plane',num2str(iPlane)]);
        [x,~] = mkdir(output_folder);
        if ~x
            error(['Error making folder to output suite2p data to for experiment ',expIDs{iExp}]);
        end
        % save matlab data
        F = Fall.F;
        Fneu = Fall.Fneu;
        spks = Fall.spks;
        ops = Fall.ops;
        stat = Fall.stat;
        exp_frames = exp_Start_frame:exp_Start_frame+frames_in_exp-1;
        save(fullfile(output_folder,'Fall.mat'),'F','Fneu','spks','ops','stat','exp_frames');
        % save npy data
        writeNPY(Fall.F,fullfile(output_folder,'F.npy'));
        writeNPY(npy.Fneu,fullfile(output_folder,'Fneu.npy'));
        writeNPY(npy.spks,fullfile(output_folder,'spks.npy'));
        % copy npy cell classification file
        copyfile(fullfile(exp_path,['plane',num2str(iPlane)],'iscell.npy'),fullfile(output_folder,'iscell.npy'));
        % copy frames from registered video bin file to split folder
        path_to_source_bin = fullfile(exp_path,['plane',num2str(iPlane)],'data.bin');
        path_to_dest_bin = fullfile(output_folder,'data.bin');
        frameSize = size(ops.meanImg);
        frames_to_copy = exp_Start_frame:exp_Start_frame+frames_in_exp-1;
        all_mean_frames{iPlane+1,iExp} = split_s2p_vid(path_to_source_bin,path_to_dest_bin,frameSize,frames_to_copy);
    end 
end

figure;
subplot(length(expIDs),length(planes_list),1);
current_plot = 0;
for iPlane = 0:length(planes_list)-1
    for iExp = 1:length(expIDs)
        current_plot = current_plot + 1;
        subplot(length(expIDs),length(planes_list),current_plot);
        imagesc(all_mean_frames{iPlane+1,iExp});
        colormap gray
    end
end
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

function writeNPY(var, filename)
% function writeNPY(var, filename)
%
% Only writes little endian, fortran (column-major) ordering; only writes
% with NPY version number 1.0.
%
% Always outputs a shape according to matlab's convention, e.g. (10, 1)
% rather than (10,).


shape = size(var);
dataType = class(var);

header = constructNPYheader(dataType, shape);

fid = fopen(filename, 'w');
fwrite(fid, header, 'uint8');
fwrite(fid, var, dataType);
fclose(fid);
end

function header = constructNPYheader(dataType, shape, varargin)

	if ~isempty(varargin)
		fortranOrder = varargin{1}; % must be true/false
		littleEndian = varargin{2}; % must be true/false
	else
		fortranOrder = true;
		littleEndian = true;
	end

    dtypesMatlab = {'uint8','uint16','uint32','uint64','int8','int16','int32','int64','single','double', 'logical'};
    dtypesNPY = {'u1', 'u2', 'u4', 'u8', 'i1', 'i2', 'i4', 'i8', 'f4', 'f8', 'b1'};

    magicString = uint8([147 78 85 77 80 89]); %x93NUMPY
    
    majorVersion = uint8(1);
    minorVersion = uint8(0);

    % build the dict specifying data type, array order, endianness, and
    % shape
    dictString = '{''descr'': ''';
    
    if littleEndian
        dictString = [dictString '<'];
    else
        dictString = [dictString '>'];
    end
    
    dictString = [dictString dtypesNPY{strcmp(dtypesMatlab,dataType)} ''', '];
    
    dictString = [dictString '''fortran_order'': '];
    
    if fortranOrder
        dictString = [dictString 'True, '];
    else
        dictString = [dictString 'False, '];
    end
    
    dictString = [dictString '''shape'': ('];
    
%     if length(shape)==1 && shape==1
%         
%     else
%         for s = 1:length(shape)
%             if s==length(shape) && shape(s)==1
%                 
%             else
%                 dictString = [dictString num2str(shape(s))];
%                 if length(shape)>1 && s+1==length(shape) && shape(s+1)==1
%                     dictString = [dictString ','];
%                 elseif length(shape)>1 && s<length(shape)
%                     dictString = [dictString ', '];
%                 end            
%             end
%         end
%         if length(shape)==1
%             dictString = [dictString ','];
%         end
%     end

    for s = 1:length(shape)
        dictString = [dictString num2str(shape(s))];
        if s<length(shape)
            dictString = [dictString ', '];
        end
    end
    
    dictString = [dictString '), '];
    
    dictString = [dictString '}'];
    
    totalHeaderLength = length(dictString)+10; % 10 is length of magicString, version, and headerLength
    
    headerLengthPadded = ceil(double(totalHeaderLength+1)/16)*16; % the whole thing should be a multiple of 16
                                                                  % I add 1 to the length in order to allow for the newline character

	% format specification is that headerlen is little endian. I believe it comes out so using this command...
    headerLength = typecast(int16(headerLengthPadded-10), 'uint8');
	
    zeroPad = zeros(1,headerLengthPadded-totalHeaderLength, 'uint8')+uint8(32); % +32 so they are spaces
    zeroPad(end) = uint8(10); % newline character
    
    header = uint8([magicString majorVersion minorVersion headerLength dictString zeroPad]);

end

function combined_data = split_s2p_vid(path_to_source_bin,path_to_dest_bin,frameSize,frames_to_copy,expected_total_frames)
frames_to_copy = double(frames_to_copy);
% fclose all
% pathToBinary = 'D:\ACC Data\2016-05-28_02_CFEB014\suite2p\plane0\data.bin';
% frameSize = [256 256];
blockSize = 1000;

finfo = dir(path_to_source_bin);
fsize = finfo.bytes;
fid = fopen(path_to_source_bin);
fid2 = fopen(path_to_dest_bin,'w');
frameCountCalculation = fsize/frameSize(1)/frameSize(2)/2;
total_frames_to_write = length(frames_to_copy);

frame_mean = [];
framesInSet = [];

% jump forward in file to start of current experiment
start_bytes = 2 * frameSize(1) * frameSize(2) * (frames_to_copy(1)-1);
fseek(fid,start_bytes,'bof');
read_data = int16(1);
for iStart = 1:blockSize:total_frames_to_write
    lastFrame = iStart + blockSize-1;
    lastFrame = min(lastFrame,total_frames_to_write);
    framesToRead = lastFrame - iStart + 1;
    disp(['Frame ',num2str(iStart+frames_to_copy(1)-1),'-',num2str(lastFrame+frames_to_copy(1)-1),' of ',num2str(frameCountCalculation)]);
    % read block of frames
    disp('Reading...');
    read_data = int16(fread(fid,frameSize(1)*frameSize(2)*framesToRead,'int16'));
    % write to other file
    disp('Writing...');
    fwrite(fid2,read_data,'int16');
    % debug
    if iStart == 1
        combined_data = reshape(read_data,[frameSize(1),frameSize(2),framesToRead]);
    %combined_data = reshape(combined_data,[frameSize(1),frameSize(2),framesToRead]);
    %figure; imagesc(mean(combined_data,3))
    end
end
fclose(fid);
fclose(fid2);
combined_data = squeeze(mean(combined_data,3));
end

