load osCategorizedTable.mat

% NRvnTable = osCategorizedTable(strcmp(osCategorizedTable.condition,"NR"),:);
NRvnTable = osCategorizedTable(strcmp(osCategorizedTable.condition,"NR") | strcmp(osCategorizedTable.condition,"DR"),:);
dataTable = NRvnTable;

DistFromMainAxis = 0;
DistAlongMainAxis = 0;

%Only mainDir OS
dataTable = dataTable(dataTable.idxOSdir == 1,:);
% dataTable = dataTable(strcmp(dataTable.ooIDX,"tON"),:);
% dataTable = dataTable(dataTable.distanceFromON > 1000,:);
% dataTable = dataTable(dataTable.OSIsigCircVar > 0.99,:);
% dataTable = dataTable(dataTable.OSIcircVar > 0.1,:);

%% Load the table
data_guide_name = 'NRvnLocations.xlsx';

% Detect import options for the data guide spreadsheet
opts = detectImportOptions(data_guide_name);

% Ensure that the 'data_name' column of the table is a string array
opts = setvartype(opts,{'fileName', 'dirVentral','Location'},'string');

expTable = readtable(data_guide_name, opts');

dataTable.OSrealTheta(strcmp(dataTable.location,"ventroTemporal")) = dataTable.OSrealTheta(strcmp(dataTable.location,"ventroTemporal")).*-1;


%% go through FOVs

allAngles = [];
allMeanAngles = [];
toVentralAngles = [];
toVentralAnglesAll = [];
awayVentralAngles = [];
awayVentralAnglesAll = [];

VNangles = [];
VTangles = [];

listFovs = unique(dataTable.fileName);
for i = 1:length(listFovs)
    tempTable = dataTable(strcmp(dataTable.fileName,listFovs(i)),:);
    meanTheta = mean(tempTable.OSrealTheta);


    
    if expTable.Location(i) == "VN"
        VNangles = [VNangles;meanTheta];
        tempTable.OSrealTheta = tempTable.OSrealTheta-(-0.8338);
        meanTheta = meanTheta -(-0.8338);
    elseif expTable.Location(i) == "VT"
        VTangles = [VTangles;meanTheta];
        tempTable.OSrealTheta = tempTable.OSrealTheta-(-0.7969);
        meanTheta = meanTheta -(-0.7969);
    end


    allAngles = [allAngles;tempTable.OSrealTheta];
    allMeanAngles = [allMeanAngles;meanTheta];

    if expTable.distVT(i)> DistFromMainAxis && expTable.distVN(i) > DistAlongMainAxis
        toVentralAngles = [toVentralAngles;meanTheta];
        toVentralAnglesAll = [toVentralAnglesAll;tempTable.OSrealTheta];
    elseif expTable.distVT(i) < -DistFromMainAxis && expTable.distVN(i) > DistAlongMainAxis
        awayVentralAngles = [awayVentralAngles;meanTheta];
        awayVentralAnglesAll = [awayVentralAnglesAll;tempTable.OSrealTheta];
    end


end

BindEdges = [-50:1:75];
% BindEdges = [-50:1:0];
% BindEdges = [0:1:75];

figure,
subplot(2,1,1)
histogram(rad2deg(toVentralAnglesAll),'BinEdges',BindEdges)
subplot(2,1,2)
histogram(rad2deg(awayVentralAnglesAll),'BinEdges',BindEdges)

toVentHistCounts = histcounts(rad2deg(toVentralAnglesAll),'BinEdges',BindEdges);
awayVentHistCounts = histcounts(rad2deg(awayVentralAnglesAll),'BinEdges',BindEdges);

figure, hold
plot(toVentHistCounts./max(toVentHistCounts),'k');
plot(awayVentHistCounts./max(awayVentHistCounts),'r');

figure, hold
cdfplot(toVentHistCounts./max(toVentHistCounts));
cdfplot(awayVentHistCounts./max(awayVentHistCounts));

% figure, violin({rad2deg(toVentralAnglesAll), rad2deg(awayVentralAnglesAll)})


[h,p,ks2stat] = kstest2(toVentralAnglesAll, awayVentralAnglesAll);

lengthToVent = length(toVentralAnglesAll);
lengthAwayVent = length(awayVentralAnglesAll);
allAngles = [toVentralAnglesAll;awayVentralAnglesAll];
lengthAngles = length(allAngles);

pRand = nan(1000,1);
kstatRand = nan(1000,1);

for i = 1:1000
    tempIdx = randperm(lengthAngles);
    
    [dummy,pRand(i),kstatRand(i)] = kstest2(allAngles(tempIdx(1:lengthToVent)), allAngles(tempIdx(lengthToVent+1:lengthAngles)));

end

% cummulative sum plots
BindEdges = [-50:1:75];

POStoVentHistCounts = histcounts(rad2deg(toVentralAnglesAll),'BinEdges',BindEdges);
POSawayVentHistCounts = histcounts(rad2deg(awayVentralAnglesAll),'BinEdges',BindEdges);

figure, hold
plot((cumsum((POStoVentHistCounts))/max(cumsum(POStoVentHistCounts))),'r');
plot((cumsum((POSawayVentHistCounts))/max(cumsum(POSawayVentHistCounts))),'k');






% figure, hold
% cdfplot(fliplr(toVentHistCounts./max(toVentHistCounts)));
% cdfplot(fliplr(awayVentHistCounts./max(awayVentHistCounts)));


