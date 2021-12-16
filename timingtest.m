global eye_check_data;
currentFrame = round(handles.slider2.Value);
handles.frame_num_text.String = ['Frame number: ',num2str(currentFrame)];
leftImage = eye_check_data.v_left.read(currentFrame);
rightImage = eye_check_data.v_right.read(currentFrame);
% DO LEFT EYE
axes(handles.axes1);
imagesc(squeeze(leftImage(:,:,1)));
axis off
colormap gray
eyeX = eye_check_data.dlc0.left(:,25:3:end);
eyeY = eye_check_data.dlc0.left(:,26:3:end);
pupilX = eye_check_data.dlc0.left(:,1:3:24);
pupilY = eye_check_data.dlc0.left(:,2:3:24);
% generate points for lids
xVals = linspace(eyeX(currentFrame,1),eyeX(currentFrame,3));
topLid = eye_check_data.dlc1.left.eyeDat.topLid{currentFrame};
botLid = eye_check_data.dlc1.left.eyeDat.botLid{currentFrame};
% for upper lid
yVals = topLid.p1*xVals.^2 + topLid.p2*xVals + topLid.p3;
% for lower lid
yVals = [yVals,botLid.p1*fliplr(xVals).^2 + botLid.p2*fliplr(xVals) + botLid.p3];
xVals = [xVals,fliplr(xVals)];
% draw the eye outline on
hold on;
plot(xVals,yVals,'y');
% plot points
% plot valid pupil points
inEye = eye_check_data.dlc1.left.eyeDat.inEye{currentFrame};
scatter(pupilX(currentFrame,inEye),pupilY(currentFrame,inEye),'g');
% invalid
scatter(pupilX(currentFrame,~inEye),pupilY(currentFrame,~inEye),'r');
% draw pupil
viscircles([eye_check_data.dlc1.left.eyeDat.x(currentFrame) eye_check_data.dlc1.left.eyeDat.y(currentFrame)],eye_check_data.dlc1.left.eyeDat.radius(currentFrame));
% DO RIGHT EYE
axes(handles.axes2);
imagesc(squeeze(rightImage(:,:,1)));
axis off
colormap gray
eyeX = eye_check_data.dlc0.right(:,25:3:end);
eyeY = eye_check_data.dlc0.right(:,26:3:end);
pupilX = eye_check_data.dlc0.right(:,1:3:24);
pupilY = eye_check_data.dlc0.right(:,2:3:24);
% generate points for lids
xVals = linspace(eyeX(currentFrame,1),eyeX(currentFrame,3));
topLid = eye_check_data.dlc1.right.eyeDat.topLid{currentFrame};
botLid = eye_check_data.dlc1.right.eyeDat.botLid{currentFrame};
% for upper lid
yVals = topLid.p1*xVals.^2 + topLid.p2*xVals + topLid.p3;
% for lower lid
yVals = [yVals,botLid.p1*fliplr(xVals).^2 + botLid.p2*fliplr(xVals) + botLid.p3];
xVals = [xVals,fliplr(xVals)];
% draw the eye outline on
hold on;
plot(xVals,yVals,'y');
% plot points
% plot valid pupil points
inEye = eye_check_data.dlc1.right.eyeDat.inEye{currentFrame};
scatter(pupilX(currentFrame,inEye),pupilY(currentFrame,inEye),'g');
% invalid
scatter(pupilX(currentFrame,~inEye),pupilY(currentFrame,~inEye),'r');
% draw pupil
viscircles([eye_check_data.dlc1.right.eyeDat.x(currentFrame) eye_check_data.dlc1.right.eyeDat.y(currentFrame)],eye_check_data.dlc1.right.eyeDat.radius(currentFrame));