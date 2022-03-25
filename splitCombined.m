
function splitCombined()

% Information to edit:
% data root
data_root = 'D:\Data';
first_exp_id = '2022-01-21_01_ESPM039';

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
        save(fullfile(output_folder,'Fall.mat'),'Fall');
        % save npy data
        writeNPY(Fall.F,fullfile(output_folder,'F.npy'));
        writeNPY(npy.Fneu,fullfile(output_folder,'Fneu.npy'));
        writeNPY(npy.spks,fullfile(output_folder,'spks.npy'));
        % copy npy cell classification file
        copyfile(fullfile(exp_path,['plane',num2str(iPlane)],'iscell.npy'),fullfile(output_folder,'iscell.npy'));
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
