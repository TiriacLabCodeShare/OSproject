%% plot OS tuning numbers

%Table of contents
%
% 1-initial loading
% 2-Tuning: just NR, all cells
% 3-Tuning: EO vs NR vs DR vs B2KO, avgs per FOV
% 4-Proportion of subtypes per FOV and across conditions
% 5-Proportions: Specifically just Directions.
% 6-Proportion: ON/OFF subtype prop comparisons against Gcamps




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

% neuronTable = neuronTable(strcmp(neuronTable.location,"ventroNasal"),:);
% neuronTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal"),:);

NRtable = neuronTable(neuronTable.age > 29 & strcmp(neuronTable.condition,"NR"),:);
DRtable = neuronTable(neuronTable.age > 29 & strcmp(neuronTable.condition,"DR"),:);
B2KOtable = neuronTable(neuronTable.age > 29 & strcmp(neuronTable.condition,"B2KO"),:);
EOtable = neuronTable(neuronTable.age < 15 & strcmp(neuronTable.condition,"EO"),:);



% neuronTable = neuronTable(neuronTable.distanceFromON < 1000,:);

% figure, histogram(neuronTable.OSrealTheta, 100);


tempWvfs = [NRtable.wvfRespToBars(strcmp(NRtable.location,"ventroNasal") & NRtable.idxOSdir == 1 & strcmp(NRtable.ooIDX,"sON"),:);...
    NRtable.wvfRespToBars(strcmp(NRtable.location,"ventroNasal") & NRtable.idxOSdir == 2 & strcmp(NRtable.ooIDX,"sON"),:);...
    NRtable.wvfRespToBars(strcmp(NRtable.location,"ventroNasal") & NRtable.idxOSdir == 1 & strcmp(NRtable.ooIDX,"tON"),:);...
    NRtable.wvfRespToBars(strcmp(NRtable.location,"ventroNasal") & NRtable.idxOSdir == 2 & strcmp(NRtable.ooIDX,"tON"),:);...
    NRtable.wvfRespToBars(strcmp(NRtable.location,"ventroNasal") & NRtable.idxOSdir == 1 & strcmp(NRtable.ooIDX,"tOFF"),:);...
    NRtable.wvfRespToBars(strcmp(NRtable.location,"ventroNasal") & NRtable.idxOSdir == 2 & strcmp(NRtable.ooIDX,"tOFF"),:);...
    ];

figure, 
imagesc(tempWvfs)
colormap copper
caxis([0 0.3])

tempWvfs = [NRtable.wvfRespToBars(strcmp(NRtable.location,"ventroTemporal") & NRtable.idxOSdir == 1 & strcmp(NRtable.ooIDX,"sON"),:);...
    NRtable.wvfRespToBars(strcmp(NRtable.location,"ventroTemporal") & NRtable.idxOSdir == 2 & strcmp(NRtable.ooIDX,"sON"),:);...
    NRtable.wvfRespToBars(strcmp(NRtable.location,"ventroTemporal") & NRtable.idxOSdir == 1 & strcmp(NRtable.ooIDX,"tON"),:);...
    NRtable.wvfRespToBars(strcmp(NRtable.location,"ventroTemporal") & NRtable.idxOSdir == 2 & strcmp(NRtable.ooIDX,"tON"),:);...
    NRtable.wvfRespToBars(strcmp(NRtable.location,"ventroTemporal") & NRtable.idxOSdir == 1 & strcmp(NRtable.ooIDX,"tOFF"),:);...
    NRtable.wvfRespToBars(strcmp(NRtable.location,"ventroTemporal") & NRtable.idxOSdir == 2 & strcmp(NRtable.ooIDX,"tOFF"),:);...
    ];

figure, 
imagesc(tempWvfs)
colormap copper
caxis([0 0.3])





%% 2- all NR OS tuning


% Make the groups to plot
x0 = NRtable.OSIcircVar;
x1 = NRtable.OSIcircVar(NRtable.idxOSdir == mainDir);
x2 = NRtable.OSIcircVar(NRtable.idxOSdir == orthoDir);
x3 = NRtable.OSIcircVar(strcmp(NRtable.ooIDX, "tON"));
x4 = NRtable.OSIcircVar(strcmp(NRtable.ooIDX, "sON"));
x5 = NRtable.OSIcircVar(strcmp(NRtable.ooIDX, "tOFF"));
x6 = NRtable.OSIcircVar(NRtable.idxOSdir == mainDir & strcmp(NRtable.ooIDX, "tON"));
x7 = NRtable.OSIcircVar(NRtable.idxOSdir == orthoDir & strcmp(NRtable.ooIDX, "tON"));
x8 = NRtable.OSIcircVar(NRtable.idxOSdir == mainDir & strcmp(NRtable.ooIDX, "sON"));
x9 = NRtable.OSIcircVar(NRtable.idxOSdir == orthoDir & strcmp(NRtable.ooIDX, "sON"));
x10 = NRtable.OSIcircVar(NRtable.idxOSdir == mainDir & strcmp(NRtable.ooIDX, "tOFF"));
x11 = NRtable.OSIcircVar(NRtable.idxOSdir == orthoDir & strcmp(NRtable.ooIDX, "tOFF"));
x = [x0; x1; x2; x3; x4; x5; x6 ; x7; x8; x9; x10; x11];

g0 = repmat({'all OS'},length(x0),1);
g1 = repmat({'mainDir'},length(x1),1);
g2 = repmat({'orthoDir'},length(x2),1);
g3 = repmat({'tON'},length(x3),1);
g4 = repmat({'sON'},length(x4),1);
g5 = repmat({'tOFF'},length(x5),1);
g6 = repmat({'tON_mainDir'},length(x6),1);
g7 = repmat({'tON_orthoDir'},length(x7),1);
g8 = repmat({'sON_mainDir'},length(x8),1);
g9 = repmat({'sON_orthoDir'},length(x9),1);
g10 = repmat({'tOFF_mainDir'},length(x10),1);
g11 = repmat({'tOFF_orthoDir'},length(x11),1);
g = [g0; g1; g2; g3; g4; g5; g6; g7; g8; g9; g10; g11];

f1 = figure('Name','NR proportions, boxplot');
boxplot(x,g,'Notch','on','Symbol', ".", 'OutlierSize',4)
ylim([0 0.5])



% figure
% swarmchart(NRtable.idxOSdir,NRtable.OSIcircVar,'.');
% ylim([0 0.5])
% 
% figure
% swarmchart(NRtable.idxOS,NRtable.OSIcircVar,'.');
% ylim([0 0.5])




%% 3 - PER FOV tuning and across exp conditions

[tuningPerFOV_NRtable, NR_listFOVs] = tuningPerFOV(NRtable);
[tuningPerFOV_DRtable, DR_listFOVs] = tuningPerFOV(DRtable);
[tuningPerFOV_B2KOtable, B2KO_listFOVs] = tuningPerFOV(B2KOtable);
[tuningPerFOV_EOtable, EO_listFOVs] = tuningPerFOV(EOtable);

NR_meanMat = [tuningPerFOV_NRtable.tONhor;tuningPerFOV_NRtable.tONver;tuningPerFOV_NRtable.sONhor;tuningPerFOV_NRtable.sONver;tuningPerFOV_NRtable.tOFFhor;tuningPerFOV_NRtable.tOFFver];
NR_categoriesMat = [ones(length(NR_listFOVs),1);2*ones(length(NR_listFOVs),1);3*ones(length(NR_listFOVs),1);4*ones(length(NR_listFOVs),1);5*ones(length(NR_listFOVs),1);6*ones(length(NR_listFOVs),1)];
DR_meanMat = [tuningPerFOV_DRtable.tONhor;tuningPerFOV_DRtable.tONver;tuningPerFOV_DRtable.sONhor;tuningPerFOV_DRtable.sONver;tuningPerFOV_DRtable.tOFFhor;tuningPerFOV_DRtable.tOFFver];
DR_categoriesMat = [ones(length(DR_listFOVs),1);2*ones(length(DR_listFOVs),1);3*ones(length(DR_listFOVs),1);4*ones(length(DR_listFOVs),1);5*ones(length(DR_listFOVs),1);6*ones(length(DR_listFOVs),1)];
B2KO_meanMat = [tuningPerFOV_B2KOtable.tONhor;tuningPerFOV_B2KOtable.tONver;tuningPerFOV_B2KOtable.sONhor;tuningPerFOV_B2KOtable.sONver;tuningPerFOV_B2KOtable.tOFFhor;tuningPerFOV_B2KOtable.tOFFver];
B2KO_categoriesMat = [ones(length(B2KO_listFOVs),1);2*ones(length(B2KO_listFOVs),1);3*ones(length(B2KO_listFOVs),1);4*ones(length(B2KO_listFOVs),1);5*ones(length(B2KO_listFOVs),1);6*ones(length(B2KO_listFOVs),1)];
EO_meanMat = [tuningPerFOV_EOtable.tONhor;tuningPerFOV_EOtable.tONver;tuningPerFOV_EOtable.sONhor;tuningPerFOV_EOtable.sONver;tuningPerFOV_EOtable.tOFFhor;tuningPerFOV_EOtable.tOFFver];
EO_categoriesMat = [ones(length(EO_listFOVs),1);2*ones(length(EO_listFOVs),1);3*ones(length(EO_listFOVs),1);4*ones(length(EO_listFOVs),1);5*ones(length(EO_listFOVs),1);6*ones(length(EO_listFOVs),1)];

neuronTable = NRtable;

%NR
tuningAllCells = [(tuningPerFOV_NRtable.tONhor+tuningPerFOV_NRtable.tOFFhor+tuningPerFOV_NRtable.sONhor)./3;...
    (tuningPerFOV_NRtable.tONver+tuningPerFOV_NRtable.tOFFver+tuningPerFOV_NRtable.sONver)./3;...
    (tuningPerFOV_NRtable.tONhor+tuningPerFOV_NRtable.tONver)./2;...
    (tuningPerFOV_NRtable.sONhor+tuningPerFOV_NRtable.sONver)./2;...
    (tuningPerFOV_NRtable.tOFFhor+tuningPerFOV_NRtable.tOFFver)./2;...
    NR_meanMat];
tuningAllCellsCat = [ones(length(NR_listFOVs),1);2*ones(length(NR_listFOVs),1);3*ones(length(NR_listFOVs),1);4*ones(length(NR_listFOVs),1);5*ones(length(NR_listFOVs),1);6*ones(length(NR_listFOVs),1);...
    7*ones(length(NR_listFOVs),1);8*ones(length(NR_listFOVs),1);9*ones(length(NR_listFOVs),1);10*ones(length(NR_listFOVs),1);11*ones(length(NR_listFOVs),1)];
tuningAllCellsViolin = [(tuningPerFOV_NRtable.tONhor+tuningPerFOV_NRtable.tOFFhor+tuningPerFOV_NRtable.sONhor)./3,...
    (tuningPerFOV_NRtable.tONver+tuningPerFOV_NRtable.tOFFver+tuningPerFOV_NRtable.sONver)./3,...
    (tuningPerFOV_NRtable.tONhor+tuningPerFOV_NRtable.tONver)./2,...
    (tuningPerFOV_NRtable.sONhor+tuningPerFOV_NRtable.sONver)./2,...
    (tuningPerFOV_NRtable.tOFFhor+tuningPerFOV_NRtable.tOFFver)./2,...
    tuningPerFOV_NRtable.tONhor,...
    tuningPerFOV_NRtable.tONver,...
    tuningPerFOV_NRtable.sONhor,...
    tuningPerFOV_NRtable.sONver,...
    tuningPerFOV_NRtable.tOFFhor,...
    tuningPerFOV_NRtable.tOFFver];

figTuningAllCells = figure('Name','tuning for all cells');
% plotSpread(tuningMat,'distributionIdx',tuningMatCat)
swarmchart(tuningAllCellsCat,tuningAllCells,'.')
ylim([0 0.3])

ftuningAllCellsViolin = figure('Name','NR tuning all cells_violin');
violin(tuningAllCellsViolin)
ylim([0 0.3])

%StatsForNRtuning
tempTable = tuningPerFOV_NRtable;
lengthTable = length(tempTable.tOFFhor);
statsTableTuningNR = table;
statsTableTuningNR.data = [tempTable.tOFFhor;tempTable.tOFFver;tempTable.tONhor;tempTable.tONver;tempTable.sONhor;tempTable.sONver];
statsTableTuningNR.dir = [repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1]);repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1]);repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1])];
statsTableTuningNR.subtype = [repmat("tOFF",[2*lengthTable 1]);repmat("tON",[2*lengthTable 1]);repmat("sON",[2*lengthTable 1])];
% [p,tbl,stats] = anovan(statsTableTuningNR.data,{statsTableTuningNR.dir,statsTableTuningNR.subtype},'model','interaction','varnames',{'dir','subtype'});
statsTableTuningNR.condition = repmat("NR",[length(statsTableTuningNR.data),1]);

%DR
tuningAllCells = [(tuningPerFOV_DRtable.tONhor+tuningPerFOV_DRtable.tOFFhor+tuningPerFOV_DRtable.sONhor)./3;...
    (tuningPerFOV_DRtable.tONver+tuningPerFOV_DRtable.tOFFver+tuningPerFOV_DRtable.sONver)./3;...
    (tuningPerFOV_DRtable.tONhor+tuningPerFOV_DRtable.tONver)./2;...
    (tuningPerFOV_DRtable.sONhor+tuningPerFOV_DRtable.sONver)./2;...
    (tuningPerFOV_DRtable.tOFFhor+tuningPerFOV_DRtable.tOFFver)./2;...
    DR_meanMat];
tuningAllCellsCat = [ones(length(DR_listFOVs),1);2*ones(length(DR_listFOVs),1);3*ones(length(DR_listFOVs),1);4*ones(length(DR_listFOVs),1);5*ones(length(DR_listFOVs),1);6*ones(length(DR_listFOVs),1);...
    7*ones(length(DR_listFOVs),1);8*ones(length(DR_listFOVs),1);9*ones(length(DR_listFOVs),1);10*ones(length(DR_listFOVs),1);11*ones(length(DR_listFOVs),1)];
tuningAllCellsViolin = [(tuningPerFOV_DRtable.tONhor+tuningPerFOV_DRtable.tOFFhor+tuningPerFOV_DRtable.sONhor)./3,...
    (tuningPerFOV_DRtable.tONver+tuningPerFOV_DRtable.tOFFver+tuningPerFOV_DRtable.sONver)./3,...
    (tuningPerFOV_DRtable.tONhor+tuningPerFOV_DRtable.tONver)./2,...
    (tuningPerFOV_DRtable.sONhor+tuningPerFOV_DRtable.sONver)./2,...
    (tuningPerFOV_DRtable.tOFFhor+tuningPerFOV_DRtable.tOFFver)./2,...
    tuningPerFOV_DRtable.tONhor,...
    tuningPerFOV_DRtable.tONver,...
    tuningPerFOV_DRtable.sONhor,...
    tuningPerFOV_DRtable.sONver,...
    tuningPerFOV_DRtable.tOFFhor,...
    tuningPerFOV_DRtable.tOFFver];

figTuningAllCells = figure('Name','DR tuning for all cells');
% plotSpread(tuningMat,'distributionIdx',tuningMatCat)
swarmchart(tuningAllCellsCat,tuningAllCells,'.')
ylim([0 0.3])

ftuningAllCellsViolin = figure('Name','DR tuning all cells_violin');
violin(tuningAllCellsViolin)
ylim([0 0.3])

%StatsForDRtuning
tempTable = tuningPerFOV_DRtable;
lengthTable = length(tempTable.tOFFhor);
statsTableTuningDR = table;
statsTableTuningDR.data = [tempTable.tOFFhor;tempTable.tOFFver;tempTable.tONhor;tempTable.tONver;tempTable.sONhor;tempTable.sONver];
statsTableTuningDR.dir = [repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1]);repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1]);repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1])];
statsTableTuningDR.subtype = [repmat("tOFF",[2*lengthTable 1]);repmat("tON",[2*lengthTable 1]);repmat("sON",[2*lengthTable 1])];
% [p,tbl,stats] = anovan(statsTableTuningDR.data,{statsTableTuningDR.dir,statsTableTuningDR.subtype},'model','interaction','varnames',{'dir','subtype'});
statsTableTuningDR.condition = repmat("DR",[length(statsTableTuningDR.data),1]);

%EO
tuningAllCells = [(tuningPerFOV_EOtable.tONhor+tuningPerFOV_EOtable.tOFFhor+tuningPerFOV_EOtable.sONhor)./3;...
    (tuningPerFOV_EOtable.tONver+tuningPerFOV_EOtable.tOFFver+tuningPerFOV_EOtable.sONver)./3;...
    (tuningPerFOV_EOtable.tONhor+tuningPerFOV_EOtable.tONver)./2;...
    (tuningPerFOV_EOtable.sONhor+tuningPerFOV_EOtable.sONver)./2;...
    (tuningPerFOV_EOtable.tOFFhor+tuningPerFOV_EOtable.tOFFver)./2;...
    EO_meanMat];
tuningAllCellsCat = [ones(length(EO_listFOVs),1);2*ones(length(EO_listFOVs),1);3*ones(length(EO_listFOVs),1);4*ones(length(EO_listFOVs),1);5*ones(length(EO_listFOVs),1);6*ones(length(EO_listFOVs),1);...
    7*ones(length(EO_listFOVs),1);8*ones(length(EO_listFOVs),1);9*ones(length(EO_listFOVs),1);10*ones(length(EO_listFOVs),1);11*ones(length(EO_listFOVs),1)];
tuningAllCellsViolin = [(tuningPerFOV_EOtable.tONhor+tuningPerFOV_EOtable.tOFFhor+tuningPerFOV_EOtable.sONhor)./3,...
    (tuningPerFOV_EOtable.tONver+tuningPerFOV_EOtable.tOFFver+tuningPerFOV_EOtable.sONver)./3,...
    (tuningPerFOV_EOtable.tONhor+tuningPerFOV_EOtable.tONver)./2,...
    (tuningPerFOV_EOtable.sONhor+tuningPerFOV_EOtable.sONver)./2,...
    (tuningPerFOV_EOtable.tOFFhor+tuningPerFOV_EOtable.tOFFver)./2,...
    tuningPerFOV_EOtable.tONhor,...
    tuningPerFOV_EOtable.tONver,...
    tuningPerFOV_EOtable.sONhor,...
    tuningPerFOV_EOtable.sONver,...
    tuningPerFOV_EOtable.tOFFhor,...
    tuningPerFOV_EOtable.tOFFver];

figTuningAllCells = figure('Name','EO tuning for all cells');
% plotSpread(tuningMat,'distributionIdx',tuningMatCat)
swarmchart(tuningAllCellsCat,tuningAllCells,'.')
ylim([0 0.3])

ftuningAllCellsViolin = figure('Name','EO tuning all cells_violin');
violin(tuningAllCellsViolin)
ylim([0 0.3])

%B2KO
tuningAllCells = [(tuningPerFOV_B2KOtable.tONhor+tuningPerFOV_B2KOtable.tOFFhor+tuningPerFOV_B2KOtable.sONhor)./3;...
    (tuningPerFOV_B2KOtable.tONver+tuningPerFOV_B2KOtable.tOFFver+tuningPerFOV_B2KOtable.sONver)./3;...
    (tuningPerFOV_B2KOtable.tONhor+tuningPerFOV_B2KOtable.tONver)./2;...
    (tuningPerFOV_B2KOtable.sONhor+tuningPerFOV_B2KOtable.sONver)./2;...
    (tuningPerFOV_B2KOtable.tOFFhor+tuningPerFOV_B2KOtable.tOFFver)./2;...
    B2KO_meanMat];
tuningAllCellsCat = [ones(length(B2KO_listFOVs),1);2*ones(length(B2KO_listFOVs),1);3*ones(length(B2KO_listFOVs),1);4*ones(length(B2KO_listFOVs),1);5*ones(length(B2KO_listFOVs),1);6*ones(length(B2KO_listFOVs),1);...
    7*ones(length(B2KO_listFOVs),1);8*ones(length(B2KO_listFOVs),1);9*ones(length(B2KO_listFOVs),1);10*ones(length(B2KO_listFOVs),1);11*ones(length(B2KO_listFOVs),1)];
tuningAllCellsViolin = [(tuningPerFOV_B2KOtable.tONhor+tuningPerFOV_B2KOtable.tOFFhor+tuningPerFOV_B2KOtable.sONhor)./3,...
    (tuningPerFOV_B2KOtable.tONver+tuningPerFOV_B2KOtable.tOFFver+tuningPerFOV_B2KOtable.sONver)./3,...
    (tuningPerFOV_B2KOtable.tONhor+tuningPerFOV_B2KOtable.tONver)./2,...
    (tuningPerFOV_B2KOtable.sONhor+tuningPerFOV_B2KOtable.sONver)./2,...
    (tuningPerFOV_B2KOtable.tOFFhor+tuningPerFOV_B2KOtable.tOFFver)./2,...
    tuningPerFOV_B2KOtable.tONhor,...
    tuningPerFOV_B2KOtable.tONver,...
    tuningPerFOV_B2KOtable.sONhor,...
    tuningPerFOV_B2KOtable.sONver,...
    tuningPerFOV_B2KOtable.tOFFhor,...
    tuningPerFOV_B2KOtable.tOFFver];

figTuningAllCells = figure('Name','B2KO tuning for all cells');
% plotSpread(tuningMat,'distributionIdx',tuningMatCat)
swarmchart(tuningAllCellsCat,tuningAllCells,'.')
ylim([0 0.3])

ftuningAllCellsViolin = figure('Name','B2KO tuning all cells_violin');
violin(tuningAllCellsViolin)
ylim([0 0.3])

%StatsForB2KOtuning
tempTable = tuningPerFOV_B2KOtable;
lengthTable = length(tempTable.tOFFhor);
statsTableTuningB2KO = table;
statsTableTuningB2KO.data = [tempTable.tOFFhor;tempTable.tOFFver;tempTable.tONhor;tempTable.tONver;tempTable.sONhor;tempTable.sONver];
statsTableTuningB2KO.dir = [repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1]);repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1]);repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1])];
statsTableTuningB2KO.subtype = [repmat("tOFF",[2*lengthTable 1]);repmat("tON",[2*lengthTable 1]);repmat("sON",[2*lengthTable 1])];
% [p,tbl,stats] = anovan(statsTableTuningB2KO.data,{statsTableTuningB2KO.dir,statsTableTuningDR.subtype},'model','interaction','varnames',{'dir','subtype'});
statsTableTuningB2KO.condition = repmat("B2KO",[length(statsTableTuningB2KO.data),1]);

f2 = figure('Name','NR tunings');
swarmchart(NR_categoriesMat,NR_meanMat);
hold
plot(1:6,[nanmean(tuningPerFOV_NRtable.tONhor) nanmean(tuningPerFOV_NRtable.tONver) nanmean(tuningPerFOV_NRtable.sONhor) nanmean(tuningPerFOV_NRtable.sONver) nanmean(tuningPerFOV_NRtable.tOFFhor) nanmean(tuningPerFOV_NRtable.tOFFver)], 'dg')
xlabel('OS subtype');
xticklabels({'','tONmainDir','tONorthDir','sONmainDir','sONorthDir','sOFFmainDir','sOFForthDir',})
ylabel('Circular variance');
ylim([0 0.5])

f3 = figure('Name','DR tunings');
swarmchart(DR_categoriesMat,DR_meanMat);
hold
plot(1:6,[nanmean(tuningPerFOV_DRtable.tONhor) nanmean(tuningPerFOV_DRtable.tONver) nanmean(tuningPerFOV_DRtable.sONhor) nanmean(tuningPerFOV_DRtable.sONver) nanmean(tuningPerFOV_DRtable.tOFFhor) nanmean(tuningPerFOV_DRtable.tOFFver)], 'dg')
xlabel('OS subtype');
xticklabels({'','tONmainDir','tONorthDir','sONmainDir','sONorthDir','sOFFmainDir','sOFForthDir',})
ylabel('Circular variance');
ylim([0 0.5])

f4 = figure('Name','B2KO tunings');
swarmchart(B2KO_categoriesMat,B2KO_meanMat);
hold
plot(1:6,[nanmean(tuningPerFOV_B2KOtable.tONhor) nanmean(tuningPerFOV_B2KOtable.tONver) nanmean(tuningPerFOV_B2KOtable.sONhor) nanmean(tuningPerFOV_B2KOtable.sONver) nanmean(tuningPerFOV_B2KOtable.tOFFhor) nanmean(tuningPerFOV_B2KOtable.tOFFver)], 'dg')
xlabel('OS subtype');
xticklabels({'','tONmainDir','tONorthDir','sONmainDir','sONorthDir','sOFFmainDir','sOFForthDir',})
ylabel('Circular variance');
ylim([0 0.5])

f5 = figure('Name','EO tunings');
swarmchart(EO_categoriesMat,EO_meanMat);
hold
plot(1:6,[nanmean(tuningPerFOV_EOtable.tONhor) nanmean(tuningPerFOV_EOtable.tONver) nanmean(tuningPerFOV_EOtable.sONhor) nanmean(tuningPerFOV_EOtable.sONver) nanmean(tuningPerFOV_EOtable.tOFFhor) nanmean(tuningPerFOV_EOtable.tOFFver)], 'dg')
xlabel('OS subtype');
xticklabels({'','tONmainDir','tONorthDir','sONmainDir','sONorthDir','sOFFmainDir','sOFForthDir',})
ylabel('Circular variance');
ylim([0 0.5])

%Stats for tuning, all conditions
statsTableTuning = [statsTableTuningNR;statsTableTuningDR;statsTableTuningB2KO];
% [p,tbl,stats] = anovan(statsTableTuning.data,{statsTableTuning.dir,statsTableTuning.subtype,statsTableTuning.condition},'model','interaction','varnames',{'dir','subtype','direction'});
% [results,~,~,gnames] = multcompare(stats,"Dimension",[2]);



%% 4 - Proportion of subtypes per FOV and across conditions

[props_NRtable, NR_listFOVs] = cellsPerFOV(NRtable);
[props_DRtable, DR_listFOVs] = cellsPerFOV(DRtable);
[props_B2KOtable, B2KO_listFOVs] = cellsPerFOV(B2KOtable);
[props_EOtable, EO_listFOVs] = cellsPerFOV(EOtable);

NR_props = [props_NRtable.tONhor./props_NRtable.numOScells;props_NRtable.tONver./props_NRtable.numOScells;props_NRtable.sONhor./props_NRtable.numOScells;props_NRtable.sONver./props_NRtable.numOScells;props_NRtable.tOFFhor./props_NRtable.numOScells;props_NRtable.tOFFver./props_NRtable.numOScells;];
DR_props = [props_DRtable.tONhor./props_DRtable.numOScells;props_DRtable.tONver./props_DRtable.numOScells;props_DRtable.sONhor./props_DRtable.numOScells;props_DRtable.sONver./props_DRtable.numOScells;props_DRtable.tOFFhor./props_DRtable.numOScells;props_DRtable.tOFFver./props_DRtable.numOScells;];
B2KO_props = [props_B2KOtable.tONhor./props_B2KOtable.numOScells;props_B2KOtable.tONver./props_B2KOtable.numOScells;props_B2KOtable.sONhor./props_B2KOtable.numOScells;props_B2KOtable.sONver./props_B2KOtable.numOScells;props_B2KOtable.tOFFhor./props_B2KOtable.numOScells;props_B2KOtable.tOFFver./props_B2KOtable.numOScells;];
EO_props = [props_EOtable.tONhor./props_EOtable.numOScells;props_EOtable.tONver./props_EOtable.numOScells;props_EOtable.sONhor./props_EOtable.numOScells;props_EOtable.sONver./props_EOtable.numOScells;props_EOtable.tOFFhor./props_EOtable.numOScells;props_EOtable.tOFFver./props_EOtable.numOScells;];

%NR
propsAllCells = [(props_NRtable.tONhor+props_NRtable.tOFFhor+props_NRtable.sONhor)./props_NRtable.numOScells;...
    (props_NRtable.tONver+props_NRtable.tOFFver+props_NRtable.sONver)./props_NRtable.numOScells;...
    (props_NRtable.tONhor+props_NRtable.tONver)./props_NRtable.numOScells;...
    (props_NRtable.sONhor+props_NRtable.sONver)./props_NRtable.numOScells;...
    (props_NRtable.tOFFhor+props_NRtable.tOFFver)./props_NRtable.numOScells;...
    NR_props];
propsAllCellsViolin = [(props_NRtable.tONhor+props_NRtable.tOFFhor+props_NRtable.sONhor)./props_NRtable.numOScells,...
    (props_NRtable.tONver+props_NRtable.tOFFver+props_NRtable.sONver)./props_NRtable.numOScells,...
    (props_NRtable.tONhor+props_NRtable.tONver)./props_NRtable.numOScells,...
    (props_NRtable.sONhor+props_NRtable.sONver)./props_NRtable.numOScells,...
    (props_NRtable.tOFFhor+props_NRtable.tOFFver)./props_NRtable.numOScells,...
    (props_NRtable.tONhor)./props_NRtable.numOScells,...
    (props_NRtable.tONver)./props_NRtable.numOScells,...
    (props_NRtable.sONhor)./props_NRtable.numOScells,...
    (props_NRtable.sONver)./props_NRtable.numOScells,...
    (props_NRtable.tOFFhor)./props_NRtable.numOScells,...
    (props_NRtable.tOFFver)./props_NRtable.numOScells];
propsAllCellsCat = [ones(length(NR_listFOVs),1);2*ones(length(NR_listFOVs),1);3*ones(length(NR_listFOVs),1);4*ones(length(NR_listFOVs),1);5*ones(length(NR_listFOVs),1);6*ones(length(NR_listFOVs),1);...
    7*ones(length(NR_listFOVs),1);8*ones(length(NR_listFOVs),1);9*ones(length(NR_listFOVs),1);10*ones(length(NR_listFOVs),1);11*ones(length(NR_listFOVs),1)];

fpropsAllCells = figure('Name','NR proportions all cells');
plotSpread(propsAllCells,'distributionIdx',propsAllCellsCat)
ylim([0 1])

fpropsAllCellsViolin = figure('Name','NR proportions all cells_violin');
violin(propsAllCellsViolin)
ylim([0 1])

%StatsForNRprops
tempTable = props_NRtable;
tempTable.tOFFhor = tempTable.tOFFhor./tempTable.numOScells;
tempTable.tOFFver = tempTable.tOFFver./tempTable.numOScells;
tempTable.tONhor = tempTable.tONhor./tempTable.numOScells;
tempTable.tONver = tempTable.tONver./tempTable.numOScells;
tempTable.sONhor = tempTable.sONhor./tempTable.numOScells;
tempTable.sONver = tempTable.sONver./tempTable.numOScells;
lengthTable = length(tempTable.numOScells);
statsTableNR = table;
statsTableNR.data = [tempTable.tOFFhor;tempTable.tOFFver;tempTable.tONhor;tempTable.tONver;tempTable.sONhor;tempTable.sONver];
statsTableNR.dir = [repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1]);repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1]);repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1])];
statsTableNR.subtype = [repmat("tOFF",[2*lengthTable 1]);repmat("tON",[2*lengthTable 1]);repmat("sON",[2*lengthTable 1])];
statsTableNR.condition = repmat("NR",[length(statsTableNR.data) 1]);
% [p,tbl,stats] = anovan(statsTableNR.data,{statsTableNR.dir,statsTableNR.subtype},'model','interaction','varnames',{'dir','subtype'});

%DR
propsAllCells = [(props_DRtable.tONhor+props_DRtable.tOFFhor+props_DRtable.sONhor)./props_DRtable.numOScells;...
    (props_DRtable.tONver+props_DRtable.tOFFver+props_DRtable.sONver)./props_DRtable.numOScells;...
    (props_DRtable.tONhor+props_DRtable.tONver)./props_DRtable.numOScells;...
    (props_DRtable.sONhor+props_DRtable.sONver)./props_DRtable.numOScells;...
    (props_DRtable.tOFFhor+props_DRtable.tOFFver)./props_DRtable.numOScells;...
    DR_props];
propsAllCellsViolin = [(props_DRtable.tONhor+props_DRtable.tOFFhor+props_DRtable.sONhor)./props_DRtable.numOScells,...
    (props_DRtable.tONver+props_DRtable.tOFFver+props_DRtable.sONver)./props_DRtable.numOScells,...
    (props_DRtable.tONhor+props_DRtable.tONver)./props_DRtable.numOScells,...
    (props_DRtable.sONhor+props_DRtable.sONver)./props_DRtable.numOScells,...
    (props_DRtable.tOFFhor+props_DRtable.tOFFver)./props_DRtable.numOScells,...
    (props_DRtable.tONhor)./props_DRtable.numOScells,...
    (props_DRtable.tONver)./props_DRtable.numOScells,...
    (props_DRtable.sONhor)./props_DRtable.numOScells,...
    (props_DRtable.sONver)./props_DRtable.numOScells,...
    (props_DRtable.tOFFhor)./props_DRtable.numOScells,...
    (props_DRtable.tOFFver)./props_DRtable.numOScells];
propsAllCellsCat = [ones(length(DR_listFOVs),1);2*ones(length(DR_listFOVs),1);3*ones(length(DR_listFOVs),1);4*ones(length(DR_listFOVs),1);5*ones(length(DR_listFOVs),1);6*ones(length(DR_listFOVs),1);...
    7*ones(length(DR_listFOVs),1);8*ones(length(DR_listFOVs),1);9*ones(length(DR_listFOVs),1);10*ones(length(DR_listFOVs),1);11*ones(length(DR_listFOVs),1)];

fpropsAllCells = figure('Name','DR proportions all cells');
plotSpread(propsAllCells,'distributionIdx',propsAllCellsCat)
ylim([0 1])

fpropsAllCellsViolin = figure('Name','DR proportions all cells_violin');
violin(propsAllCellsViolin)
ylim([0 1])

%StatsForDRprops
tempTable = props_DRtable;
tempTable.tOFFhor = tempTable.tOFFhor./tempTable.numOScells;
tempTable.tOFFver = tempTable.tOFFver./tempTable.numOScells;
tempTable.tONhor = tempTable.tONhor./tempTable.numOScells;
tempTable.tONver = tempTable.tONver./tempTable.numOScells;
tempTable.sONhor = tempTable.sONhor./tempTable.numOScells;
tempTable.sONver = tempTable.sONver./tempTable.numOScells;
lengthTable = length(tempTable.numOScells);
statsTableDR = table;
statsTableDR.data = [tempTable.tOFFhor;tempTable.tOFFver;tempTable.tONhor;tempTable.tONver;tempTable.sONhor;tempTable.sONver];
statsTableDR.dir = [repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1]);repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1]);repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1])];
statsTableDR.subtype = [repmat("tOFF",[2*lengthTable 1]);repmat("tON",[2*lengthTable 1]);repmat("sON",[2*lengthTable 1])];
statsTableDR.condition = repmat("DR",[length(statsTableDR.data) 1]);

%EO
propsAllCells = [(props_EOtable.tONhor+props_EOtable.tOFFhor+props_EOtable.sONhor)./props_EOtable.numOScells;...
    (props_EOtable.tONver+props_EOtable.tOFFver+props_EOtable.sONver)./props_EOtable.numOScells;...
    (props_EOtable.tONhor+props_EOtable.tONver)./props_EOtable.numOScells;...
    (props_EOtable.sONhor+props_EOtable.sONver)./props_EOtable.numOScells;...
    (props_EOtable.tOFFhor+props_EOtable.tOFFver)./props_EOtable.numOScells;...
    EO_props];
propsAllCellsViolin = [(props_EOtable.tONhor+props_EOtable.tOFFhor+props_EOtable.sONhor)./props_EOtable.numOScells,...
    (props_EOtable.tONver+props_EOtable.tOFFver+props_EOtable.sONver)./props_EOtable.numOScells,...
    (props_EOtable.tONhor+props_EOtable.tONver)./props_EOtable.numOScells,...
    (props_EOtable.sONhor+props_EOtable.sONver)./props_EOtable.numOScells,...
    (props_EOtable.tOFFhor+props_EOtable.tOFFver)./props_EOtable.numOScells,...
    (props_EOtable.tONhor)./props_EOtable.numOScells,...
    (props_EOtable.tONver)./props_EOtable.numOScells,...
    (props_EOtable.sONhor)./props_EOtable.numOScells,...
    (props_EOtable.sONver)./props_EOtable.numOScells,...
    (props_EOtable.tOFFhor)./props_EOtable.numOScells,...
    (props_EOtable.tOFFver)./props_EOtable.numOScells];
propsAllCellsCat = [ones(length(EO_listFOVs),1);2*ones(length(EO_listFOVs),1);3*ones(length(EO_listFOVs),1);4*ones(length(EO_listFOVs),1);5*ones(length(EO_listFOVs),1);6*ones(length(EO_listFOVs),1);...
    7*ones(length(EO_listFOVs),1);8*ones(length(EO_listFOVs),1);9*ones(length(EO_listFOVs),1);10*ones(length(EO_listFOVs),1);11*ones(length(EO_listFOVs),1)];

fpropsAllCells = figure('Name','EO proportions all cells');
plotSpread(propsAllCells,'distributionIdx',propsAllCellsCat)
ylim([0 1])

fpropsAllCellsViolin = figure('Name','EO proportions all cells_violin');
violin(propsAllCellsViolin)
ylim([0 1])

%B2KO
propsAllCells = [(props_B2KOtable.tONhor+props_B2KOtable.tOFFhor+props_B2KOtable.sONhor)./props_B2KOtable.numOScells;...
    (props_B2KOtable.tONver+props_B2KOtable.tOFFver+props_B2KOtable.sONver)./props_B2KOtable.numOScells;...
    (props_B2KOtable.tONhor+props_B2KOtable.tONver)./props_B2KOtable.numOScells;...
    (props_B2KOtable.sONhor+props_B2KOtable.sONver)./props_B2KOtable.numOScells;...
    (props_B2KOtable.tOFFhor+props_B2KOtable.tOFFver)./props_B2KOtable.numOScells;...
    B2KO_props];
propsAllCellsViolin = [(props_B2KOtable.tONhor+props_B2KOtable.tOFFhor+props_B2KOtable.sONhor)./props_B2KOtable.numOScells,...
    (props_B2KOtable.tONver+props_B2KOtable.tOFFver+props_B2KOtable.sONver)./props_B2KOtable.numOScells,...
    (props_B2KOtable.tONhor+props_B2KOtable.tONver)./props_B2KOtable.numOScells,...
    (props_B2KOtable.sONhor+props_B2KOtable.sONver)./props_B2KOtable.numOScells,...
    (props_B2KOtable.tOFFhor+props_B2KOtable.tOFFver)./props_B2KOtable.numOScells,...
    (props_B2KOtable.tONhor)./props_B2KOtable.numOScells,...
    (props_B2KOtable.tONver)./props_B2KOtable.numOScells,...
    (props_B2KOtable.sONhor)./props_B2KOtable.numOScells,...
    (props_B2KOtable.sONver)./props_B2KOtable.numOScells,...
    (props_B2KOtable.tOFFhor)./props_B2KOtable.numOScells,...
    (props_B2KOtable.tOFFver)./props_B2KOtable.numOScells];
propsAllCellsCat = [ones(length(B2KO_listFOVs),1);2*ones(length(B2KO_listFOVs),1);3*ones(length(B2KO_listFOVs),1);4*ones(length(B2KO_listFOVs),1);5*ones(length(B2KO_listFOVs),1);6*ones(length(B2KO_listFOVs),1);...
    7*ones(length(B2KO_listFOVs),1);8*ones(length(B2KO_listFOVs),1);9*ones(length(B2KO_listFOVs),1);10*ones(length(B2KO_listFOVs),1);11*ones(length(B2KO_listFOVs),1)];

fpropsAllCells = figure('Name','B2KO proportions all cells');
plotSpread(propsAllCells,'distributionIdx',propsAllCellsCat)
ylim([0 1])

fpropsAllCellsViolin = figure('Name','B2KO proportions all cells_violin');
violin(propsAllCellsViolin)
ylim([0 1])

%StatsForB2KOprops
tempTable = props_B2KOtable;
tempTable.tOFFhor = tempTable.tOFFhor./tempTable.numOScells;
tempTable.tOFFver = tempTable.tOFFver./tempTable.numOScells;
tempTable.tONhor = tempTable.tONhor./tempTable.numOScells;
tempTable.tONver = tempTable.tONver./tempTable.numOScells;
tempTable.sONhor = tempTable.sONhor./tempTable.numOScells;
tempTable.sONver = tempTable.sONver./tempTable.numOScells;
lengthTable = length(tempTable.numOScells);
statsTableB2KO = table;
statsTableB2KO.data = [tempTable.tOFFhor;tempTable.tOFFver;tempTable.tONhor;tempTable.tONver;tempTable.sONhor;tempTable.sONver];
statsTableB2KO.dir = [repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1]);repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1]);repmat("horz",[lengthTable 1]);repmat("vert",[lengthTable 1])];
statsTableB2KO.subtype = [repmat("tOFF",[2*lengthTable 1]);repmat("tON",[2*lengthTable 1]);repmat("sON",[2*lengthTable 1])];
statsTableB2KO.condition = repmat("B2KO",[length(statsTableB2KO.data) 1]);

% swarmchart(propsAllCellsCat,propsAllCells);
% hold
% % plot(1:6,[nanmean(props_NRtable.tONhor./props_NRtable.numOScells) nanmean(props_NRtable.tONver./props_NRtable.numOScells) nanmean(props_NRtable.sONhor./props_NRtable.numOScells) nanmean(props_NRtable.sONver./props_NRtable.numOScells) nanmean(props_NRtable.tOFFhor./props_NRtable.numOScells) nanmean(props_NRtable.tOFFver./props_NRtable.numOScells)], 'dg')
% xlabel('OS subtype');
% xticklabels({'','main','orth','tON','sON','tOFF','tONmainDir','tONorthDir','sONmainDir','sONorthDir','sOFFmainDir','sOFForthDir',})
% ylabel('Proportion of OS cells');
% ylim([0 1])


f6 = figure('Name','NR proportions');
swarmchart(NR_categoriesMat,NR_props);
hold
plot(1:6,[nanmean(props_NRtable.tONhor./props_NRtable.numOScells) nanmean(props_NRtable.tONver./props_NRtable.numOScells) nanmean(props_NRtable.sONhor./props_NRtable.numOScells) nanmean(props_NRtable.sONver./props_NRtable.numOScells) nanmean(props_NRtable.tOFFhor./props_NRtable.numOScells) nanmean(props_NRtable.tOFFver./props_NRtable.numOScells)], 'dg')
xlabel('OS subtype');
xticklabels({'','tONmainDir','tONorthDir','sONmainDir','sONorthDir','sOFFmainDir','sOFForthDir',})
ylabel('Proportion of OS cells');
ylim([0 1])

f7 = figure('Name','DR proportions');
swarmchart(DR_categoriesMat,DR_props);
hold
plot(1:6,[nanmean(props_DRtable.tONhor./props_DRtable.numOScells) nanmean(props_DRtable.tONver./props_DRtable.numOScells) nanmean(props_DRtable.sONhor./props_DRtable.numOScells) nanmean(props_DRtable.sONver./props_DRtable.numOScells) nanmean(props_DRtable.tOFFhor./props_DRtable.numOScells) nanmean(props_DRtable.tOFFver./props_DRtable.numOScells)], 'dg')
xlabel('OS subtype');
xticklabels({'','tONmainDir','tONorthDir','sONmainDir','sONorthDir','sOFFmainDir','sOFForthDir',})
ylabel('Proportion of OS cells');
ylim([0 1])

f8 = figure('Name','B2KO proportions');
swarmchart(B2KO_categoriesMat,B2KO_props);
hold
plot(1:6,[nanmean(props_B2KOtable.tONhor./props_B2KOtable.numOScells) nanmean(props_B2KOtable.tONver./props_B2KOtable.numOScells) nanmean(props_B2KOtable.sONhor./props_B2KOtable.numOScells) nanmean(props_B2KOtable.sONver./props_B2KOtable.numOScells) nanmean(props_B2KOtable.tOFFhor./props_B2KOtable.numOScells) nanmean(props_B2KOtable.tOFFver./props_B2KOtable.numOScells)], 'dg')
xlabel('OS subtype');
xticklabels({'','tONmainDir','tONorthDir','sONmainDir','sONorthDir','sOFFmainDir','sOFForthDir',})
ylabel('Proportion of OS cells');
ylim([0 1])

f9 = figure('Name','EO proportions');
swarmchart(EO_categoriesMat,EO_props);
hold
plot(1:6,[nanmean(props_EOtable.tONhor./props_EOtable.numOScells) nanmean(props_EOtable.tONver./props_EOtable.numOScells) nanmean(props_EOtable.sONhor./props_EOtable.numOScells) nanmean(props_EOtable.sONver./props_EOtable.numOScells) nanmean(props_EOtable.tOFFhor./props_EOtable.numOScells) nanmean(props_EOtable.tOFFver./props_EOtable.numOScells)], 'dg')
xlabel('OS subtype');
xticklabels({'','tONmainDir','tONorthDir','sONmainDir','sONorthDir','sOFFmainDir','sOFForthDir',})
ylabel('Proportion of OS cells');
ylim([0 1])

%Stats for all proportions
statsAllProps = [statsTableNR;statsTableDR;statsTableB2KO];
[p,tbl,stats] = anovan(statsAllProps.data,{statsAllProps.dir,statsAllProps.subtype,statsAllProps.condition},'model','interaction','varnames',{'dir','subtype','condition'});
[results,~,~,gnames] = multcompare(stats,"Dimension",[2]);


return

%% 5 - Proportions: Focusing on just directions

props_NRtable = props_EOtable;
NR_listFOVs = EO_listFOVs;

props_NRtable.horCells = props_NRtable.tOFFhor + props_NRtable.tONhor + props_NRtable.sONhor;
props_NRtable.verCells = props_NRtable.tOFFver + props_NRtable.tONver + props_NRtable.sONver;

NR_dirAnalysis = [props_NRtable.horCells./props_NRtable.numOScells; props_NRtable.verCells./props_NRtable.numOScells;...
    props_NRtable.tONhor./(props_NRtable.tONhor+props_NRtable.tONver);props_NRtable.tONver./(props_NRtable.tONhor+props_NRtable.tONver);...
    props_NRtable.sONhor./(props_NRtable.sONhor+props_NRtable.sONver);props_NRtable.sONver./(props_NRtable.sONhor+props_NRtable.sONver);...
    props_NRtable.tOFFhor./(props_NRtable.tOFFhor+props_NRtable.tOFFver);props_NRtable.tOFFver./(props_NRtable.tOFFhor+props_NRtable.tOFFver)];
NR_dirCategories = [ones(length(NR_listFOVs),1); 2*ones(length(NR_listFOVs),1);...
    3*ones(length(NR_listFOVs),1); 4*ones(length(NR_listFOVs),1);...
    5*ones(length(NR_listFOVs),1); 6*ones(length(NR_listFOVs),1);...
    7*ones(length(NR_listFOVs),1); 8*ones(length(NR_listFOVs),1)];

figure('Name','EO Direction proportion analysis');
swarmchart(NR_dirCategories,NR_dirAnalysis);
hold
plot(1:8,[nanmean(props_NRtable.horCells./props_NRtable.numOScells) nanmean(props_NRtable.verCells./props_NRtable.numOScells)...
    nanmean(props_NRtable.tONhor./(props_NRtable.tONhor+props_NRtable.tONver)) nanmean(props_NRtable.tONver./(props_NRtable.tONhor+props_NRtable.tONver))...
    nanmean(props_NRtable.sONhor./(props_NRtable.sONhor+props_NRtable.sONver)) nanmean(props_NRtable.sONver./(props_NRtable.sONhor+props_NRtable.sONver))...
    nanmean(props_NRtable.tOFFhor./(props_NRtable.tOFFhor+props_NRtable.tOFFver)) nanmean(props_NRtable.tOFFver./(props_NRtable.tOFFhor+props_NRtable.tOFFver))], 'dr')
xlabel('OS subtype');
xticklabels({'','horz','vert','tONhor','tONver','sONhor','sONver','tOFFhor','tOFFver'})
ylabel('Proportion of OS cells');
ylim([0 1])

return


%% 6 - ON/OFF proportion comparisons against gcamps

load gcampOScategorizedTable.mat

[propsOOonly_NR, NR_listFOVs] = cellsPerFOV_ooONLY(NRtable);
[propsOOonly_G6, G6_listFOVs] = cellsPerFOV_ooONLY(gcampOScategorizedTable);

propData = [propsOOonly_NR.tON./propsOOonly_NR.numOScells; propsOOonly_G6.tON./propsOOonly_G6.numOScells; propsOOonly_NR.sON./propsOOonly_NR.numOScells;propsOOonly_G6.sON./propsOOonly_G6.numOScells; propsOOonly_NR.tOFF./propsOOonly_NR.numOScells;propsOOonly_G6.tOFF./propsOOonly_G6.numOScells];
propCategories = [ones(length(NR_listFOVs),1);2*ones(length(G6_listFOVs),1); 3*ones(length(NR_listFOVs),1);4*ones(length(G6_listFOVs),1); 5*ones(length(NR_listFOVs),1);6*ones(length(G6_listFOVs),1)];

figure('Name','NR vs GCAMP proportions out of OS cells');
swarmchart(propCategories,propData);
hold
plot(1:6,[nanmean(propsOOonly_NR.tON./propsOOonly_NR.numOScells) nanmean(propsOOonly_G6.tON./propsOOonly_G6.numOScells) nanmean(propsOOonly_NR.sON./propsOOonly_NR.numOScells) nanmean(propsOOonly_G6.sON./propsOOonly_G6.numOScells) nanmean(propsOOonly_NR.tOFF./propsOOonly_NR.numOScells) nanmean(propsOOonly_G6.tOFF./propsOOonly_G6.numOScells)], 'dg')
xlabel('OS subtype');
xticklabels({'','NR tON','G6 tON', 'NR sON','G6 sON', 'NR tOFF','G6 tOFF',})
ylabel('Proportion of OS cells');
ylim([0 1])

% propDataAllCells = [propsOOonly_NR.tON./propsOOonly_NR.numCells; propsOOonly_G6.tON./propsOOonly_G6.numCells; propsOOonly_NR.sON./propsOOonly_NR.numCells;propsOOonly_G6.sON./propsOOonly_G6.numCells; propsOOonly_NR.tOFF./propsOOonly_NR.numCells;propsOOonly_G6.tOFF./propsOOonly_G6.numCells];
propDataAllCells = [propsOOonly_NR.tON; propsOOonly_G6.tON; propsOOonly_NR.sON;propsOOonly_G6.sON; propsOOonly_NR.tOFF;propsOOonly_G6.tOFF];

figure('Name','NR vs GCAMP proportions out of all Cells');
swarmchart(propCategories,propDataAllCells);
hold
% plot(1:6,[nanmean(propsOOonly_NR.tON./propsOOonly_NR.numCells) nanmean(propsOOonly_G6.tON./propsOOonly_G6.numCells) nanmean(propsOOonly_NR.sON./propsOOonly_NR.numCells) nanmean(propsOOonly_G6.sON./propsOOonly_G6.numCells) nanmean(propsOOonly_NR.tOFF./propsOOonly_NR.numCells) nanmean(propsOOonly_G6.tOFF./propsOOonly_G6.numCells)], 'dg')
xlabel('OS subtype');
xticklabels({'','NR tON','G6 tON', 'NR sON','G6 sON', 'NR tOFF','G6 tOFF',})
ylabel('Proportion of OS cells');
ylim([0 1])

%% 7 - Waveforms

orderedWvfs = [NRtable.wvfRespToBars(strcmp(NRtable.ooIDX,"sON") & NRtable.idxOSdir == 1,:);...
    NRtable.wvfRespToBars(strcmp(NRtable.ooIDX,"sON") & NRtable.idxOSdir == 2,:);...
    NRtable.wvfRespToBars(strcmp(NRtable.ooIDX,"tON") & NRtable.idxOSdir == 1,:);...
    NRtable.wvfRespToBars(strcmp(NRtable.ooIDX,"tON") & NRtable.idxOSdir == 2,:);...
    NRtable.wvfRespToBars(strcmp(NRtable.ooIDX,"tOFF") & NRtable.idxOSdir == 1,:);...
    NRtable.wvfRespToBars(strcmp(NRtable.ooIDX,"tOFF") & NRtable.idxOSdir == 2,:)];

figure, imagesc(orderedWvfs)
caxis([0 1])






% Stats

% Copy pasted example, follow this format:
% y = [52.7 57.5 45.9 44.5 53.0 57.0 45.9 44.0]';
% g1 = [1 2 1 2 1 2 1 2]; 
% g2 = {'hi';'hi';'lo';'lo';'hi';'hi';'lo';'lo'}; 
% g3 = {'NR';'DR';'B2KO'};
% 
% 
% p = anovan(y,{g1,g2,g3})



%% Whole population analysis

% tempTable = osCategorizedTable(osCategorizedTable.idxOS == 3,:);
% 
% figure,
% boxplot(tempTable.OSIcircVar,tempTable.condition)




%% Functions

function [resultsTable, listFOVs] = tuningPerFOV(inputTable)

    resultsTable = table;

    listFOVs = unique(inputTable.fileName);
    resultsTable.tOFFhor = nan(length(listFOVs),1);
    resultsTable.tOFFver = nan(length(listFOVs),1);
    resultsTable.tONhor = nan(length(listFOVs),1);
    resultsTable.tONver = nan(length(listFOVs),1);
    resultsTable.sONhor = nan(length(listFOVs),1);
    resultsTable.sONver = nan(length(listFOVs),1);
    
    
    for i = 1:length(listFOVs)
        tempTable = inputTable(strcmp(inputTable.fileName,listFOVs(i)),:);
    
        resultsTable.tOFFhor(i) = mean(tempTable.OSIcircVar(tempTable.idxOSdir == 1 & strcmp(tempTable.ooIDX, "tOFF")),'omitnan');
        resultsTable.tOFFver(i) = mean(tempTable.OSIcircVar(tempTable.idxOSdir == 2 & strcmp(tempTable.ooIDX, "tOFF")),'omitnan');
        resultsTable.tONhor(i) = mean(tempTable.OSIcircVar(tempTable.idxOSdir == 1 & strcmp(tempTable.ooIDX, "tON")),'omitnan');    
        resultsTable.tONver(i) = mean(tempTable.OSIcircVar(tempTable.idxOSdir == 2 & strcmp(tempTable.ooIDX, "tON")),'omitnan');
        resultsTable.sONhor(i) = mean(tempTable.OSIcircVar(tempTable.idxOSdir == 1 & strcmp(tempTable.ooIDX, "sON")),'omitnan');
        resultsTable.sONver(i) = mean(tempTable.OSIcircVar(tempTable.idxOSdir == 2 & strcmp(tempTable.ooIDX, "sON")),'omitnan');
    
    end

end

function [resultsTable, listFOVs] = cellsPerFOV(inputTable)

    resultsTable = table;

    listFOVs = unique(inputTable.fileName);
    resultsTable.tOFFhor = nan(length(listFOVs),1);
    resultsTable.tOFFver = nan(length(listFOVs),1);
    resultsTable.tONhor = nan(length(listFOVs),1);
    resultsTable.tONver = nan(length(listFOVs),1);
    resultsTable.sONhor = nan(length(listFOVs),1);
    resultsTable.sONver = nan(length(listFOVs),1);
    resultsTable.numOScells = nan(length(listFOVs),1);
    resultsTable.numCells = nan(length(listFOVs),1);
    
    
    for i = 1:length(listFOVs)
        tempTable = inputTable(strcmp(inputTable.fileName,listFOVs(i)),:);

        resultsTable.numOScells(i) = length(tempTable.DSI);
        resultsTable.numCells(i) = max(tempTable.neuronNumWithinFOV);
    
        resultsTable.tOFFhor(i) = sum(tempTable.idxOSdir == 1 & strcmp(tempTable.ooIDX, "tOFF"));
        resultsTable.tOFFver(i) = sum(tempTable.idxOSdir == 2 & strcmp(tempTable.ooIDX, "tOFF"));
        resultsTable.tONhor(i) = sum(tempTable.idxOSdir == 1 & strcmp(tempTable.ooIDX, "tON"));    
        resultsTable.tONver(i) = sum(tempTable.idxOSdir == 2 & strcmp(tempTable.ooIDX, "tON"));
        resultsTable.sONhor(i) = sum(tempTable.idxOSdir == 1 & strcmp(tempTable.ooIDX, "sON"));
        resultsTable.sONver(i) = sum(tempTable.idxOSdir == 2 & strcmp(tempTable.ooIDX, "sON"));
    
    end

    resultsTable.numOSperAllCells = resultsTable.numOScells./resultsTable.numCells;

end

function [resultsTable, listFOVs] = cellsPerFOV_ooONLY(inputTable)

    resultsTable = table;

    listFOVs = unique(inputTable.fileName);
    resultsTable.tOFF = nan(length(listFOVs),1);
    resultsTable.tON = nan(length(listFOVs),1);
    resultsTable.sON = nan(length(listFOVs),1);
    resultsTable.numOScells = nan(length(listFOVs),1);
    resultsTable.numCells = nan(length(listFOVs),1);
    
    
    for i = 1:length(listFOVs)
        tempTable = inputTable(strcmp(inputTable.fileName,listFOVs(i)),:);

        resultsTable.numOScells(i) = length(tempTable.DSI);
        resultsTable.numCells(i) = max(tempTable.neuronNumWithinFOV);
    
        resultsTable.tOFF(i) = sum(strcmp(tempTable.ooIDX, "tOFF"));
        resultsTable.tON(i) = sum(strcmp(tempTable.ooIDX, "tON"));    
        resultsTable.sON(i) = sum(strcmp(tempTable.ooIDX, "sON"));
    
    end

    resultsTable.numOSperAllCells = resultsTable.numOScells./resultsTable.numCells;

end




