%% plot OS tuning numbers

%Table of contents
%
% 1-initial loading
% 2-waveforms




%% 1-initial loading

load osCategorizedTable.mat
neuronTable = osCategorizedTable;

% load neuronTable.mat

% % Make sure these indices are correct:
mainDir = 1;
orthoDir = 2;

neuronTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal") | strcmp(neuronTable.location,"ventroNasal"),:);
neuronTable = neuronTable(strcmp(neuronTable.calciumSensor,"calDye")  | strcmp(neuronTable.calciumSensor,"calBryte"),:);

% neuronTable.OSrealTheta = neuronTable.OSrealTheta+deg2rad(neuronTable.degCorr);

neuronTable = neuronTable(strcmp(neuronTable.location,"ventroNasal"),:);
% neuronTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal"),:);

NRtable = neuronTable(neuronTable.age > 29 & strcmp(neuronTable.condition,"NR"),:);
DRtable = neuronTable(neuronTable.age > 29 & strcmp(neuronTable.condition,"DR"),:);
B2KOtable = neuronTable(neuronTable.age > 29 & strcmp(neuronTable.condition,"B2KO"),:);





% neuronTable = neuronTable(neuronTable.distanceFromON < 1000,:);

% figure, histogram(neuronTable.OSrealTheta, 100);

%% 2-heatmaps and waveforms

neuronTable = NRtable;

% Set OS thresholds??
OSIsigThresh = 0.95;
VARthresh = 1000; %Set to 10 for no Thresh; Set to 0.2 for reasonable thresh
dFoFthresh = 0;
OSIthresh = 0.1;
neuronTable = neuronTable(neuronTable.OSIcircVar > OSIthresh & neuronTable.OSIsigCircVar > OSIsigThresh &...
    neuronTable.varSum < VARthresh & max(neuronTable.meanRespToBars,[],2) > dFoFthresh,:);

concatResp = nan(length(neuronTable.DSI),20,8);
concatResp(:,:,1) = neuronTable.wvfRespToBars(:,1:20);
concatResp(:,:,2) = neuronTable.wvfRespToBars(:,21:40);
concatResp(:,:,3) = neuronTable.wvfRespToBars(:,41:60);
concatResp(:,:,4) = neuronTable.wvfRespToBars(:,61:80);
concatResp(:,:,5) = neuronTable.wvfRespToBars(:,81:100);
concatResp(:,:,6) = neuronTable.wvfRespToBars(:,101:120);
concatResp(:,:,7) = neuronTable.wvfRespToBars(:,121:140);
concatResp(:,:,8) = neuronTable.wvfRespToBars(:,141:160);

meanRespAllDir = mean(concatResp,3,'omitnan');

figure, hold
plot(mean(meanRespAllDir),'k')
plot(mean(meanRespAllDir(strcmp(neuronTable.ooIDX,"sON"),:)),'r')
plot(mean(meanRespAllDir(strcmp(neuronTable.ooIDX,"tON"),:)),'b')
plot(mean(meanRespAllDir(strcmp(neuronTable.ooIDX,"tOFF"),:)),'g')
ylim([-0.1 0.3])


figure, hold
plot(mean(neuronTable.wvfRespToBars(strcmp(neuronTable.ooIDX,"sON") & neuronTable.idxOSdir == 1,:)),'g')
plot(mean(neuronTable.wvfRespToBars(strcmp(neuronTable.ooIDX,"sON") & neuronTable.idxOSdir == 2,:)),'m')


figure, hold
plot(mean(neuronTable.wvfRespToBars(strcmp(neuronTable.ooIDX,"tON") & neuronTable.idxOSdir == 1,:)),'g')
plot(mean(neuronTable.wvfRespToBars(strcmp(neuronTable.ooIDX,"tON") & neuronTable.idxOSdir == 2,:)),'m')


figure, hold
plot(mean(neuronTable.wvfRespToBars(strcmp(neuronTable.ooIDX,"tOFF") & neuronTable.idxOSdir == 1,:)),'g')
plot(mean(neuronTable.wvfRespToBars(strcmp(neuronTable.ooIDX,"tOFF") & neuronTable.idxOSdir == 2,:)),'m')


%% Polar plots of tuning curves

%for OFFt
offPeaks = [16:20:156, 16];

offtwvfsMain = neuronTable.wvfRespToBars(strcmp(neuronTable.ooIDX,"tOFF") & neuronTable.idxOSdir == mainDir,:);
offtwvfsOrtho = neuronTable.wvfRespToBars(strcmp(neuronTable.ooIDX,"tOFF") & neuronTable.idxOSdir == orthoDir,:);
offtwvfsMainTuning = neuronTable.OSIcircVar(strcmp(neuronTable.ooIDX,"tOFF") & neuronTable.idxOSdir == mainDir,:);
offtwvfsOrthoTuning = neuronTable.OSIcircVar(strcmp(neuronTable.ooIDX,"tOFF") & neuronTable.idxOSdir == orthoDir,:);

mean_offtwvfsMain = mean(offtwvfsMain(:,offPeaks));
std_offtwvfsMain = std(offtwvfsMain(:,offPeaks));

mean_offtwvfsOrtho = mean(offtwvfsOrtho(:,offPeaks));

figure
polarplot(mean_offtwvfsMain)
rlim([0 0.3])

figure
polarplot(mean_offtwvfsOrtho)
rlim([0 0.3])




%for ONs
onPeaks = [9:20:149, 9];
onswvfsMain = neuronTable.wvfRespToBars(strcmp(neuronTable.ooIDX,"sON") & neuronTable.idxOSdir == mainDir,:);
onswvfsOrtho = neuronTable.wvfRespToBars(strcmp(neuronTable.ooIDX,"sON") & neuronTable.idxOSdir == orthoDir,:);
onswvfsMainTuning = neuronTable.OSIcircVar(strcmp(neuronTable.ooIDX,"sON") & neuronTable.idxOSdir == mainDir,:);
onswvfsOrthoTuning = neuronTable.OSIcircVar(strcmp(neuronTable.ooIDX,"sON") & neuronTable.idxOSdir == orthoDir,:);

mean_onswvfsMain = mean(onswvfsMain(:,onPeaks));
std_onswvfsMain = std(onswvfsOrtho(:,onPeaks));

mean_onswvfsOrtho = mean(onswvfsOrtho(:,onPeaks));

figure
polarplot(mean_onswvfsMain)
rlim([0 0.3])

figure
polarplot(mean_onswvfsOrtho)
rlim([0 0.3])



%for ONt
onPeaks = [9:20:149, 9];
ontwvfsMain = neuronTable.wvfRespToBars(strcmp(neuronTable.ooIDX,"tON") & neuronTable.idxOSdir == mainDir,:);
ontwvfsOrtho = neuronTable.wvfRespToBars(strcmp(neuronTable.ooIDX,"tON") & neuronTable.idxOSdir == orthoDir,:);
ontwvfsMainTuning = neuronTable.OSIcircVar(strcmp(neuronTable.ooIDX,"tON") & neuronTable.idxOSdir == mainDir,:);
ontwvfsOrthoTuning = neuronTable.OSIcircVar(strcmp(neuronTable.ooIDX,"tON") & neuronTable.idxOSdir == orthoDir,:);

mean_ontwvfsMain = mean(ontwvfsMain(:,onPeaks));
std_ontwvfsMain = std(ontwvfsOrtho(:,onPeaks));

mean_ontwvfsOrtho = mean(ontwvfsOrtho(:,onPeaks));

figure
polarplot(mean_ontwvfsMain)
rlim([0 0.3])

figure
polarplot(mean_ontwvfsOrtho)
rlim([0 0.3])



