%% Build the neuron table
% this code will build a table where every row is a neuron
% Created by Alex T on 10/11/2022
%
% This code needs an excel spreadsheet with metadata in it.
% Right now, it will read 'OSProjectMetaData.xlsx'
% Every row in that excell spreadsheet is a field of view.
% The buildNeuronTable code will run all the field of views
%
% The end result of this code is that it saves 'neuronTable.mat'
% This is the table where each row is a neuron
% Other codes will load this table and run analysis on it

%% Load the table
data_guide_name = 'OSProjectMetaData2.xlsx';

% Detect import options for the data guide spreadsheet
opts = detectImportOptions(data_guide_name);

% Ensure that the 'data_name' column of the table is a string array
opts = setvartype(opts,{'ExperimentDate', 'FileName','LightCondition','GFPLabel','eye', 'location', 'CalciumSensor','orientationTransform'},'string');

expTable = readtable(data_guide_name, opts');

[num_files dummy] = size(expTable);

%% Intialize the neuron table (called neuronTable)

totalNeurons = nansum(expTable.numberOfNeurons);
neuronCounter = 1; %This value will keep track of which neuron I am on

neuronTable = table('Size', [totalNeurons 33], 'VariableTypes', {'double','double','string','string','string',...
    'double','string','string','string','string','string','double','double','double','string','double','double',...
    'string','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double',...
    });

neuronTable.Properties.VariableNames = {'neuronNum', 'neuronNumWithinFOV', 'animalID','fileName', 'sex',...
    'age', 'condition', 'calciumSensor','GFPretina','GFPid','eye','X','Y','distanceFromON','location','degCorr','frameRate',...
    'cellID','DSI','DSIsig','vecSum','prefDir','prefDirCorr','varSum','OSI','OSTheta','OSIsig','OSrealTheta','OSrealThetaCorr','OSIcircVar','OSIsigCircVar','OSrealVS','QI',...
    };

neuronTable.QIdir = NaN*ones([totalNeurons,8]);
neuronTable.wvfRespToBars = NaN*ones([totalNeurons,160]);
neuronTable.wvfRespToBarsAllReps = NaN*ones([totalNeurons,160,3]);
neuronTable.meanRespToBars = NaN*ones([totalNeurons,9]);
neuronTable.meanRespToBarsNormalized = NaN*ones([totalNeurons,9]);
neuronTable.cellLoc = NaN*ones([totalNeurons,2]);
% neuronTable.shuffledDSI = NaN*ones([totalNeurons,1000]);
% neuronTable.shuffledOSI = NaN*ones([totalNeurons,1000]);

neuronTable.GFPid = repmat("none",[totalNeurons 1]);


%% Here is the for loop that goes through every file


for i = 1:num_files   % this will run through all files
% for i = 1:1 % For testing code

    %load file(i) and initialize some variables
    load(expTable.FileName(i)) 
    animalID = char(expTable.FileName(i));
    animalID = animalID(1:6);
    
    %Transform directions using key
    realOrientations = str2num(expTable.orientationTransform(i));

    newInd = zeros(24,1);
    replaceDir = realOrientations';
    temp = unique(textFileArray(:,1));
    
    for tt = 1:numel(temp)
        indices = textFileArray(:,1) == temp(tt);
        newInd(indices) = replaceDir(tt);
    end
    textFileArrayOld = textFileArray;
    textFileArray(:,1) = newInd; 

    % Call function calcDS to calculate DS
    [DSI, vecSum, vecTheta, OSI, OSTheta, OSrealTheta, OSIcircVar, OSrealVS, wvf_resp_mean,wvf_resp_allReps, sumVar, rhos_all, rhos_norm_all,QI, QIdir] = calcDS(roiInt, textFileArray, expTable.FrameRate(i));
   
    % Determine if DSI is significant via bootstrapping
    shuffledDSI = nan(length(DSI),1000);
    shuffledOSI = nan(length(OSI),1000);
    shuffledOSIcircVar = nan(length(OSI),1000);
    parfor j = 1:1000
        shuffledText = textFileArray;
        shuffledText(:,1) = shuffleTrialDirections(textFileArray(:,1),8, 3);

        [shuffledDSI(:,j), dummy, dummy, shuffledOSI(:,j), dummy, dummy, shuffledOSIcircVar(:,j), dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy] = calcDS(roiInt, shuffledText, expTable.FrameRate(i));
    end

    %Determine DSI significance
    [val,idx]=min(abs(DSI-sort(shuffledDSI,2)),[],2);
    DSIsigTemp = idx/1000;

    %Determine OSI significance
    [val,idx]=min(abs(OSI-sort(shuffledOSI,2)),[],2);
    OSIsigTemp = idx/1000;

    [val,idx]=min(abs(OSIcircVar-sort(shuffledOSIcircVar,2)),[],2);
    OSIsigCircVar = idx/1000;
    
%     % Call function to make the CellID string array
%     [cellID, numID] = determineID(str2num(expTable.On_OffCells(i)), str2num(expTable.OnCells(i)), str2num(expTable.BadCells(i)),expTable.numberOfNeurons(i));
    
%     % Call function to add labels to GFPid so we can identify GFP+ cells
%     if expTable.GFPLabel_(i) ~= "none"
%         [GFPid] = determineGFPid(str2num(expTable.GFPcellID(i)),expTable.numberOfNeurons(i),expTable.GFPLabel_(i));
%         neuronTable.GFPid(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = GFPid;
%     end
    

    % Load neuron info in the neuronTable
    neuronTable.neuronNum(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = [neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1]';
    neuronTable.neuronNumWithinFOV(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = [1:expTable.numberOfNeurons(i)]';
    neuronTable.animalID(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(animalID,[expTable.numberOfNeurons(i),1]);
    neuronTable.fileName(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.FileName(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.sex(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.Sex(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.age(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.Age(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.condition(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.LightCondition(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.calciumSensor(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.CalciumSensor(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.GFPretina(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.GFPLabel(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.eye(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.eye(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.X(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.x(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.Y(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.y(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.distanceFromON(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.DistanceFromONinUm(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.location(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.location(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.degCorr(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.degCorrection(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.frameRate(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.FrameRate(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.DSI(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = DSI;
    neuronTable.DSIsig(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = DSIsigTemp;
%     neuronTable.shuffledDSI(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1,:) = shuffledDSI;
    neuronTable.vecSum(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = vecSum;
    neuronTable.prefDir(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = vecTheta;
    neuronTable.varSum(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = sumVar;
    neuronTable.OSI(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = OSI;
    neuronTable.OSIsig(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = OSIsigTemp;
    neuronTable.OSIsigCircVar(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = OSIsigCircVar;
%     neuronTable.shuffledOSI(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1,:) = shuffledOSI;
    neuronTable.OSTheta(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = OSTheta;
    neuronTable.OSrealTheta(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = OSrealTheta;
    neuronTable.OSIcircVar(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = OSIcircVar;
    neuronTable.OSrealVS(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = OSrealVS;
    neuronTable.QI(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = QI;
    neuronTable.QIdir(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1,:) = QIdir;
    neuronTable.wvfRespToBars(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1,:) = wvf_resp_mean;
    neuronTable.wvfRespToBarsAllReps(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1,:,:) = wvf_resp_allReps;
    neuronTable.meanRespToBars(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1,:) = rhos_all;
    neuronTable.meanRespToBarsNormalized(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1,:) = rhos_norm_all;

    if exist('cellCenters','var') == 1
        neuronTable.cellLoc(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1,:) = cellCenters;
        clear cellCenters
    end
    
    %Update how many neurons we have analyzed so far (+1)
    neuronCounter = neuronCounter+expTable.numberOfNeurons(i);
    
    if floor(i/10)==i/10
        i
    end
end

neuronTable.prefDirCorr = neuronTable.prefDir + deg2rad(neuronTable.degCorr);
for i = 1:totalNeurons
    %Fix prefDirCorr
    if neuronTable.prefDirCorr(i) > 2*pi
        neuronTable.prefDirCorr(i) = neuronTable.prefDirCorr(i) - 2*pi;
    elseif neuronTable.prefDirCorr(i) < 0
        neuronTable.prefDirCorr(i) = neuronTable.prefDirCorr(i) + 2*pi;
    end
    
%     %Determine DSI significance
%     [val,idx]=min(abs(neuronTable.DSI(i)-sort(neuronTable.shuffledDSI(i,:))));
%     neuronTable.DSIsig(i) = idx/1000;
% 
%     %Determine OSI significance
%     [val,idx]=min(abs(neuronTable.OSI(i)-sort(neuronTable.shuffledOSI(i,:))));
%     neuronTable.OSIsig(i) = idx/1000;
end

neuronTable.OSrealThetaCorr = neuronTable.OSrealThetaCorr + deg2rad(neuronTable.degCorr);

% %classify all non-BAD hb9 and drd4 cells as ON-OFF cells
% neuronTable.cellID(strcmp(neuronTable.GFPid, 'hb9')==1 & strcmp(neuronTable.cellID,'BAD')==0) = 'ON-OFF';
% neuronTable.cellID(strcmp(neuronTable.GFPid, 'drd4')==1 & strcmp(neuronTable.cellID,'BAD')==0) = 'ON-OFF';


% save('neuronTable.mat','neuronTable');
% 
% 







%% Functions

%DS function
function [DSI,vecSum, vecTheta, OSI, OSTheta, OSrealTheta, OSIcircVar, OSrealVS, wvf_resp_mean, wvf_resp_allReps, sumVar, rhos_all, rhos_norm_all,QI,QIdir] = calcDS(roiInt, textFileArray, frameRate)

    [num_trials dummy] = size(textFileArray);
    [num_cells num_frames] = size(roiInt);

    % Find max intensity at every bar presentation (not discriminating between ON and OFF)

    timeOfStim = mean(textFileArray(:,8)-textFileArray(:,7));
    timeOfStim_frames = ceil(timeOfStim * frameRate);
    winWvf = 3;

    barResp = zeros(num_cells,num_trials);
    wvf_resp = zeros(num_cells, timeOfStim_frames+2*winWvf ,num_trials); % waveforms of all responses
    wvf_resp_allComb = [];

    for i = 1:num_trials

        stimFrameStart = floor(textFileArray(i,7)*frameRate);
        stimFrameEnd = ceil(textFileArray(i,8)*frameRate);
        barResp(:,i) = max(roiInt(:,stimFrameStart:stimFrameEnd),[],2);

        wvf_resp(:,:,i) = roiInt(:,stimFrameStart-winWvf:stimFrameStart+timeOfStim_frames+winWvf-1);
        wvf_resp_allComb = [wvf_resp_allComb; roiInt(:,stimFrameStart-winWvf:stimFrameStart+timeOfStim_frames+winWvf-1)];

    end

    ind0 = find(textFileArray(:,1) == 0);
    ind45 = find(textFileArray(:,1) == 45);
    ind90 = find(textFileArray(:,1) == 90);
    ind135 = find(textFileArray(:,1) == 135);
    ind180 = find(textFileArray(:,1) == 180);
    ind225 = find(textFileArray(:,1) == 225);
    ind270 = find(textFileArray(:,1) == 270);
    ind315 = find(textFileArray(:,1) == 315);

    %Zeros start
    OSI = zeros(num_cells,1);
    OSrealTheta = zeros(num_cells,1);
    OScircVar = zeros(num_cells,1);
    OSrealVS = zeros(num_cells,1);

    DSI = zeros(num_cells,1);
    OSTheta = zeros(num_cells,1);
    vecSum = zeros(num_cells,1);
    vecTheta = zeros(num_cells,1);
    maxDF = zeros(num_cells,1);
    sumVar = zeros(num_cells,1);

    rhos_all = zeros(num_cells, 9);
    rhos1_all = zeros(num_cells, 9);
    rhos2_all = zeros(num_cells, 9);
    rhos3_all = zeros(num_cells, 9);
    rhos_var = zeros(num_cells, 9);
    rhos_norm_all = zeros(num_cells, 9);
    %Zeroes end

    thetas = [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4, 0]; %Theta values for 0 45 90 135 180 225 270 315

    for i = 1:num_cells
        
        %Calculate OSI
        rhos = [mean([barResp(i,ind0) barResp(i,ind180)]), mean([barResp(i,ind45) barResp(i,ind225)]), mean([barResp(i,ind90) barResp(i,ind270)]),mean([barResp(i,ind135) barResp(i,ind315)]),mean([barResp(i,ind0) barResp(i,ind180)])];
        rhosAllDir = [mean(barResp(i,ind0)), mean(barResp(i,ind45)), mean(barResp(i,ind90)), mean(barResp(i,ind135)), mean(barResp(i,ind180)), mean(barResp(i,ind225)), mean(barResp(i,ind270)), mean(barResp(i,ind315)), mean(barResp(i,ind0))];
        rhosNorm = rhos/sum(rhos(1:4));

        [temp, prefInd] = max(rhosNorm(1:4));

        prefOrientation = thetas(prefInd);
        nullInd = prefInd - 2;
        if nullInd < 1
            nullInd = 4+nullInd;
        end
        nullOrientation = thetas(nullInd);
    
        OSI_temp = (rhos(prefInd)-rhos(nullInd))/(rhos(prefInd)+rhos(nullInd));
    
        OSI(i,1) = OSI_temp;
        OSTheta(i,1) = prefOrientation; % this is not the real OS dir
        
        %Code to figure out real OS direction
        OSnum = nan(1,8);
        OSdenom = nan(1,8);
        
        for jj = 1:8
            OSnum(jj) = rhosAllDir(jj)*exp(2i*thetas(jj));
            OSdenom(jj) = rhosAllDir(jj);
        end
        
        circVar = sum(OSnum)/sum(OSdenom);
        
        OSIcircVar(i,1) = abs(circVar);
        OSrealTheta(i,1) = (angle(circVar)/2);


%         if prefInd == 3
%             rhosForOStheta = [rhosNorm(1:4), rhosNorm(1)];
%             thetasForOStheta = [thetas(1:5)];
%         elseif prefInd == 4
%             rhosForOStheta = [rhosNorm(2:4), rhosNorm(1:2)];
%             thetasForOStheta = [thetas(2:6)];
%         elseif prefInd == 1
%             rhosForOStheta = [rhosNorm(3:4), rhosNorm(1:3)];
%             thetasForOStheta = [thetas(7:8), thetas(1:3)];
%         elseif prefInd == 2
%             rhosForOStheta = [rhosNorm(4), rhosNorm(1:4)];
%             thetasForOStheta = [thetas(8), thetas(1:4)];
%         end
%         [x,y] = pol2cart(thetasForOStheta,rhosForOStheta);
%         vectSumX = sum(x);
%         vectSumY = sum(y);
%         [OSrealTheta(i,1), OSrealVS(i,1)] = cart2pol(vectSumX,vectSumY);



        %Calculate DSI
        rhos = [mean(barResp(i,ind0)), mean(barResp(i,ind45)), mean(barResp(i,ind90)), mean(barResp(i,ind135)), mean(barResp(i,ind180)), mean(barResp(i,ind225)), mean(barResp(i,ind270)), mean(barResp(i,ind315)), mean(barResp(i,ind0))];
        rhos1 = [barResp(i,ind0(1)), barResp(i,ind45(1)), barResp(i,ind90(1)), barResp(i,ind135(1)), barResp(i,ind180(1)), barResp(i,ind225(1)), barResp(i,ind270(1)), barResp(i,ind315(1)), barResp(i,ind0(1))];
        rhos2 = [barResp(i,ind0(2)), barResp(i,ind45(2)), barResp(i,ind90(2)), barResp(i,ind135(2)), barResp(i,ind180(2)), barResp(i,ind225(2)), barResp(i,ind270(2)), barResp(i,ind315(2)), barResp(i,ind0(2))];
        rhos3 = [barResp(i,ind0(3)), barResp(i,ind45(3)), barResp(i,ind90(3)), barResp(i,ind135(3)), barResp(i,ind180(3)), barResp(i,ind225(3)), barResp(i,ind270(3)), barResp(i,ind315(3)), barResp(i,ind0(3))];
        rhos_var(i,:) = [var(barResp(i,ind0)), var(barResp(i,ind45)), var(barResp(i,ind90)), var(barResp(i,ind135)), var(barResp(i,ind180)), var(barResp(i,ind225)), var(barResp(i,ind270)), var(barResp(i,ind315)), var(barResp(i,ind0))];
        sumVar(i) = sum(rhos_var(i,:));
        
        
        rhosNorm = rhos/sum(rhos(1:8));
        [x,y] = pol2cart(thetas,rhosNorm);
        vectSumX = sum(x(1:8));
        vectSumY = sum(y(1:8));
        [the, rho] = cart2pol(vectSumX,vectSumY);
        if the <0
            the = 2*pi+the;
        end

        %what is the pref dir?

        if the > 5.8905
            prefInd = 1;
        else
            [temp,prefInd] = min(abs(thetas-the));
        end
        prefDir = thetas(prefInd);
        nullInd = prefInd - 4;
        if nullInd < 1
            nullInd = 8+nullInd;
        end

        DSI_temp = (rhos(prefInd)-rhos(nullInd))/(rhos(prefInd)+rhos(nullInd));


        %Store for cells
        DSI(i,1) = DSI_temp;
        vecSum(i,1) = rho;
        vecTheta(i,1) = the;
        maxDF(i,1) = max(rhos); %What was the max DF

        rhos_all(i,:) = rhos;
        rhos1_all(i,:) = rhos1;
        rhos2_all(i,:) = rhos2;
        rhos3_all(i,:) = rhos3;
        rhos_norm_all(i,:) = rhosNorm;
        %End store for cells


    end
    
    
    reordered_ind = [ind0 ; ind45 ; ind90 ; ind135 ; ind180 ; ind225 ; ind270 ; ind315];
    wvf_resp_reordered = wvf_resp(:,:,reordered_ind);
    wvf_resp_mean = [mean(wvf_resp_reordered(:,:,1:3),3),mean(wvf_resp_reordered(:,:,4:6),3),mean(wvf_resp_reordered(:,:,7:9),3),mean(wvf_resp_reordered(:,:,10:12),3),mean(wvf_resp_reordered(:,:,13:15),3),mean(wvf_resp_reordered(:,:,16:18),3),mean(wvf_resp_reordered(:,:,19:21),3),mean(wvf_resp_reordered(:,:,22:24),3)];
    wvf_resp_allReps = [wvf_resp_reordered(:,:,1:3),wvf_resp_reordered(:,:,4:6),wvf_resp_reordered(:,:,7:9),wvf_resp_reordered(:,:,10:12),wvf_resp_reordered(:,:,13:15),wvf_resp_reordered(:,:,16:18),wvf_resp_reordered(:,:,19:21),wvf_resp_reordered(:,:,22:24)];
    
    % QI calculation
    QI = nan(num_cells,1);
    QIdir = nan(num_cells,8);
    for i = 1:num_cells
        currWvf = reshape(wvf_resp_allReps(i,:,:),[160,3])';
        QI(i) = var(mean(currWvf,1))/mean(var(currWvf,0,2));

        currWvf = reshape(wvf_resp_reordered(i,:,1:3),[20,3])';
        QIdir(i,1) = var(mean(currWvf,1))/mean(var(currWvf,0,2));

        currWvf = reshape(wvf_resp_reordered(i,:,4:6),[20,3])';
        QIdir(i,2) = var(mean(currWvf,1))/mean(var(currWvf,0,2));

        currWvf = reshape(wvf_resp_reordered(i,:,7:9),[20,3])';
        QIdir(i,3) = var(mean(currWvf,1))/mean(var(currWvf,0,2));

        currWvf = reshape(wvf_resp_reordered(i,:,10:12),[20,3])';
        QIdir(i,4) = var(mean(currWvf,1))/mean(var(currWvf,0,2));

        currWvf = reshape(wvf_resp_reordered(i,:,13:15),[20,3])';
        QIdir(i,5) = var(mean(currWvf,1))/mean(var(currWvf,0,2));

        currWvf = reshape(wvf_resp_reordered(i,:,16:18),[20,3])';
        QIdir(i,6) = var(mean(currWvf,1))/mean(var(currWvf,0,2));

        currWvf = reshape(wvf_resp_reordered(i,:,19:21),[20,3])';
        QIdir(i,7) = var(mean(currWvf,1))/mean(var(currWvf,0,2));

        currWvf = reshape(wvf_resp_reordered(i,:,22:24),[20,3])';
        QIdir(i,8) = var(mean(currWvf,1))/mean(var(currWvf,0,2));
    end


end

