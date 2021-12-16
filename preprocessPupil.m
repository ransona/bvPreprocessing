function preprocessPupil(expID)
expID = '2021-09-30_02_TEST';
animalID = expID(15:end);
dataRoot = 'D:\data';
expRoot = fullfile(dataRoot,animalID,expID);

displayOn = true;
displayInterval = 100;

if displayOn
    f = figure;
end

disp(['Starting ',expID]);

dlc_filenames = {[expID,'_eye1_leftDLC_resnet50_Trial_newMay19shuffle1_1030000.csv'],...
    [expID,'_eye1_rightDLC_resnet50_Trial_newMay19shuffle1_1030000.csv']};

vid_filenames = {[expID,'_eye1_left.avi'],...
    [expID,'_eye1_right.avi']};

if exist(fullfile(expRoot,'dlcEye.mat'),'file')
    if ~strcmp(questdlg('DLC data already processed - reprocess?'),'Yes')
        return
    end
end

for iVid = 1:length(dlc_filenames)
    try
        videoPath = fullfile(expRoot,vid_filenames{iVid});
        v = VideoReader(videoPath);
    catch
        disp('Error: Eye video file not found');
        return
    end

    % read the csv deeplabcut output file
    dlc_data = readmatrix(fullfile(expRoot,dlc_filenames{iVid}));
    % remove first column
    dlc_data = dlc_data(:,2:end);
    eyeX = dlc_data(:,25:3:end);
    eyeY = dlc_data(:,26:3:end);
    pupilX = dlc_data(:,1:3:24);
    pupilY = dlc_data(:,2:3:24);
    % get minimum of eye x and eye y confidence from dlc
    % eye x and eye y are always needed as a minimum to process a frame so
    % we ensure below that these coordinates all have confidence > 0.8
    eyeMinConfid = min(dlc_data(:,27:3:end),[],2);
    
    firstFrame = readVideoIndex(v,1);
    frameSize = size(squeeze(firstFrame(:,:,1)));
    
    % choose approximate eye area - points outside this will be considered
    % invalid
    figure(f);
    imagesc(firstFrame);
    roiLeft = median(eyeX(:,3));
    roiTop = median(eyeY(:,2));
    roiWidth = median(eyeX(:,1))-roiLeft;
    roiHeight = median(eyeY(:,4))-roiTop;
    padding = roiWidth * 0.75;
    eyeROI = images.roi.Rectangle(gca,'Position',[roiLeft-padding,roiTop-padding,roiWidth+padding*2,roiHeight+padding*2],'StripeColor','r');
    eyeROI.wait;
    validRegionMask = createMask(eyeROI,firstFrame);
    eyeROI.delete;
    % calc some average values for eye to be used for QC later
    eyeWidth = (median(eyeX(:,1))-median(eyeX(:,3)));
    % clip coordinates to frame size
    eyeX(eyeX>frameSize(2))=frameSize(2);
    eyeY(eyeY>frameSize(1))=frameSize(1);
    pupilX(pupilX>frameSize(2))=frameSize(2);
    pupilY(pupilY>frameSize(1))=frameSize(1);
    eyeX(eyeX<1)=1;
    eyeY(eyeY<1)=1;
    pupilX(pupilX<1)=1;
    pupilY(pupilY<1)=1;
    
    eyeDat = [];
    lastFrame = tic;
    
    for iFrame = 1:size(dlc_data,1)
        % do QC to make sure the eye corners have been well detected
        pointsValid = validRegionMask(sub2ind(size(validRegionMask),ceil(eyeY(iFrame,:)),ceil(eyeX(iFrame,:))));
        % check spacing of left and right corners is about right compared to
        % median of whole recording
        cornerDistanceDiff = abs((eyeX(iFrame,1)-eyeX(iFrame,3))-eyeWidth)/eyeWidth;
        % check top and bottom lid mid points are around halfway between eye
        % corners in x direction
        min_corner_middle_distance = min(abs([eyeX(iFrame,3)-eyeX(iFrame,2),...
            eyeX(iFrame,3)-eyeX(iFrame,4),...
            eyeX(iFrame,1)-eyeX(iFrame,2),...
            eyeX(iFrame,1)-eyeX(iFrame,4)]));
        
        if min(pointsValid)==1 && (min_corner_middle_distance / eyeWidth)>0.33
            
            % fit two parabolas - one for each eye lid
            % points:
            % 1 = lateral / 2 = sup / 3 = medial / 4 = inf
            topLid = fit(eyeX(iFrame,[1 2 3])',eyeY(iFrame,[1 2 3])','poly2');
            botLid = fit(eyeX(iFrame,[1 4 3])',eyeY(iFrame,[1 4 3])','poly2');
            
            % generate points
            xVals = linspace(eyeX(iFrame,1),eyeX(iFrame,3));
            % for upper lid
            yVals = topLid.p1*xVals.^2 + topLid.p2*xVals + topLid.p3;
            % for lower lid
            yVals = [yVals,botLid.p1*fliplr(xVals).^2 + botLid.p2*fliplr(xVals) + botLid.p3];
            xVals = [xVals,fliplr(xVals)];
            % check if y values
            % make a poly mask using the points
            eyeMask = (poly2mask(xVals,yVals,frameSize(1),frameSize(2)));
            % check if each pupil point is in the eye mask and exclude it if not
            pupilIdx = sub2ind(frameSize,round(pupilY(iFrame,:)),round(pupilX(iFrame,:)));
            inEye = eyeMask(pupilIdx);
            
            % fit a circle to those pupil points within the eye
            if sum(inEye) > 2
                [xCenter, yCenter, radius, ~] = circlefit(pupilX(iFrame,inEye),pupilY(iFrame,inEye));
            else
                % not enough points to fit circle
                xCenter = nan;
                yCenter = nan;
                radius = nan;
            end
            
            eyeDat.x(iFrame) = xCenter;
            eyeDat.y(iFrame) = yCenter;
            
            if sum(inEye) < 2
                eyeDat.qc(iFrame) = 2; % indicates QC failed due to pupil fit
            else
                eyeDat.qc(iFrame) = 0; % indicates QC passed
            end
            
            eyeDat.radius(iFrame) = radius;
            eyeDat.topLid{iFrame} = topLid;
            eyeDat.botLid{iFrame} = botLid;
            eyeDat.inEye{iFrame} = inEye;
            
            if mod(iFrame,displayInterval)==0
                disp('#####################');
                disp([num2str(iFrame),'/',num2str(size(dlc_data,1)),' - ',num2str(iFrame/size(dlc_data,1)*100),'% complete']);
                disp(['Frame rate = ',num2str(1/toc(lastFrame))]);
                disp('#####################');
            end
            
        else
            % frame has failed QC
            eyeDat.x(iFrame) = nan;
            eyeDat.y(iFrame) = nan;
            eyeDat.radius(iFrame) = nan;
            eyeDat.topLid{iFrame} = nan;
            eyeDat.botLid{iFrame} = nan;
            eyeDat.inEye{iFrame} = nan;
            eyeDat.qc(iFrame) = 2; % indicates all eye corners not well detected
        end
        
        lastFrame = tic;
        
        if displayOn
            if mod(iFrame,displayInterval)==0
                if eyeDat.qc(iFrame) == 0
                    % quality control passed
                    set(groot,'CurrentFigure',f);
                    currentFrame = readVideoIndex(v,iFrame);
                    currentFrame = squeeze(currentFrame(:,:,1));
                    % all plotting
                    hold off
                    imagesc((currentFrame));
                    colormap gray
                    % plot all coordinates
                    hold on
                    % plot eye
                    plot(xVals,yVals,'y');
                    % plot valid pupil points
                    scatter(pupilX(iFrame,inEye),pupilY(iFrame,inEye),'g');
                    % invalid
                    scatter(pupilX(iFrame,~inEye),pupilY(iFrame,~inEye),'r');
                    % draw pupil circle
                    viscircles([xCenter yCenter],radius);
                    drawnow
                else
                    % quality control NOT passed
                    set(groot,'CurrentFigure',f);
                    currentFrame = readVideoIndex(v,iFrame);
                    currentFrame = squeeze(currentFrame(:,:,1));
                    % all plotting
                    hold off
                    imagesc((currentFrame));
                    colormap gray
                    % plot all coordinates
                    hold on
                    scatter(pupilX(iFrame,:),pupilY(iFrame,:),'g');
                    % plot valid pupil points
                    scatter(pupilX(iFrame,:),pupilY(iFrame,:),'g');
                    % eye points
                    scatter(eyeX(iFrame,:),eyeY(iFrame,:));
                    drawnow
                end
            end
        end
    end
    
    % do some further quality control
    % remove points where pupil looks dodgy
%     validFrames = ones(1,length(videoTiming.pupil.x));
%     validFrames(isnan(videoTiming.pupil.x))=0;
%     validFrames(videoTiming.pupil.radius > 100)=0;
%     validFrames(videoTiming.pupil.x > median(videoTiming.pupil.x(validFrames==1)+100))=0;
%     validFrames(videoTiming.pupil.x < median(videoTiming.pupil.x(validFrames==1)-100))=0;
%     validFrames(videoTiming.pupil.y > median(videoTiming.pupil.y(validFrames==1)+100))=0;
%     validFrames(videoTiming.pupil.y < median(videoTiming.pupil.y(validFrames==1)-100))=0;
%     videoTiming.pupil.x(validFrames==0) = nan;
%     videoTiming.pupil.y(validFrames==0) = nan;
%     videoTiming.pupil.radius(validFrames==0) = nan;
%     
     % calculate pupil velocity
     xdiffs = diff(eyeDat.x);
     ydiffs = diff(eyeDat.x);
     eucla_diff = sqrt(xdiffs.^2 + ydiffs.^2);
     eyeDat.velocity = conv(eucla_diff,ones(1,10),'same');
    eyeDat.velocity = [eyeDat.velocity,eyeDat.velocity(end)];
    if iVid == 1
        save(fullfile(expRoot,'dlcEyeLeft.mat'),'eyeDat');
    else
        save(fullfile(expRoot,'dlcEyeRight.mat'),'eyeDat');
    end

end
end

function [xCenter, yCenter, radius, a] = circlefit(x, y)
% circlefit(): Fits a circle through a set of points in the x - y plane.
% USAGE :
% [xCenter, yCenter, radius, a] = circlefit(X, Y)
% The output is the center point (xCenter, yCenter) and the radius of the fitted circle.
% "a" is an optional output vector describing the coefficients in the circle's equation:
%     x ^ 2 + y ^ 2 + a(1) * x + a(2) * y + a(3) = 0
% by Bucher Izhak 25 - Oct - 1991

numPoints = numel(x);
xx = x .* x;
yy = y .* y;
xy = x .* y;
A = [sum(x),  sum(y),  numPoints;
    sum(xy), sum(yy), sum(y);
    sum(xx), sum(xy), sum(x)];
B = [-sum(xx + yy) ;
    -sum(xx .* y + yy .* y);
    -sum(xx .* x + xy .* y)];
a = A \ B;
xCenter = -.5 * a(1);
yCenter = -.5 * a(2);
radius  =  sqrt((a(1) ^ 2 + a(2) ^ 2) / 4 - a(3));
end