% These MATLAB scripts are prepared by A.M.E for the following paper,
% Ahmet M. Elbir, "CNN-based Precoder and Combiner Design in mmWave MIMO Systems", IEEE Communications Letters, in press.
% please cite the above work if you use this file. For any comments and
% questions please email: ahmetmelbir@gmail.com

% clearvars -except fileNameLabelsTrain
% fileNameLabelsTrain = 'scenarioID_52_52'; % Train data.
% load(fileNameLabelsTrain)
%% CNN AS classification type
dataAntennaSelection = XAS;
labelsAntennaSelection = categorical(Y);
sizeInputAntennaSelection = size(dataAntennaSelection);
Qf = numel(categories(labelsAntennaSelection));

%% CNN RF for regression type
dataFRFChainSelection = XRF;
labelsRFChainSelection = YFRFr;
sizeInputFRFChainSelection = size(dataFRFChainSelection);
sizeOutputFRFChainSelection = size(labelsRFChainSelection);

dataWRFChainSelection = XRF;
labelsWRFChainSelection = YWRFr;
sizeInputWRFChainSelection = size(dataWRFChainSelection);
sizeOutputWRFChainSelection = size(labelsWRFChainSelection);
%% val. for regression.

idx = randperm(size(dataAntennaSelection,4),floor(.3*sizeInputAntennaSelection(end)));
valDataAntennaSelection = dataAntennaSelection(:,:,:,idx);
valLabelsAntennaSelection = labelsAntennaSelection(:,idx);
dataAntennaSelection(:,:,:,idx) = [];
labelsAntennaSelection(:,idx) = [];

idx = randperm(size(dataFRFChainSelection,4),floor(.3*sizeInputFRFChainSelection(end)));
valDataRFChainSelection = dataFRFChainSelection(:,:,:,idx);
valLabelsFRFChainSelection = labelsRFChainSelection(idx,:);
dataFRFChainSelection(:,:,:,idx) = [];
labelsRFChainSelection(idx,:) = [];

idx = randperm(size(dataFRFChainSelection,4),floor(.3*sizeInputFRFChainSelection(end)));
valDataWRFChainSelection = dataWRFChainSelection(:,:,:,idx);
valLabelsWRFChainSelection = labelsWRFChainSelection(idx,:);
dataWRFChainSelection(:,:,:,idx) = [];
labelsWRFChainSelection(idx,:) = [];
%% settings.
miniBatchSize = 200;
numValidationsPerEpoch = 5000;
validationFrequency = 50;
%% DNN Layers.
layersAntennaSelection = [imageInputLayer([sizeInputAntennaSelection(1:3)],'Normalization', 'zerocenter');
    convolution2dLayer(2,2^6);
    batchNormalizationLayer
    reluLayer();
%                         maxPooling2dLayer(2,'Stride',2);
    convolution2dLayer(2,2^6);
    batchNormalizationLayer
    reluLayer();
%                         maxPooling2dLayer(2,'Stride',2);
    convolution2dLayer(2,2^6);
    batchNormalizationLayer
    reluLayer();
%                     maxPooling2dLayer(2,'Stride',2);
%     convolution2dLayer(2,2^6);
%     batchNormalizationLayer
%     reluLayer();
    fullyConnectedLayer(2^7);
    reluLayer();
    dropoutLayer();
    fullyConnectedLayer(2^7);
    reluLayer();
    dropoutLayer();
    fullyConnectedLayer(Qf);
    softmaxLayer();
    classificationLayer();
    ]

    optsAntennaSelection = trainingOptions('sgdm',...
    'Momentum', 0.9,...
    'InitialLearnRate',0.01,... % The default value is 0.01.
    'MaxEpochs',5000,...
    'MiniBatchSize',miniBatchSize,... % The default is 128.
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.9,...
    'LearnRateDropPeriod',30,...
    'L2Regularization',0.0001,... % The default is 0.0001.
    'ExecutionEnvironment', 'auto',...
    'ValidationData',{valDataAntennaSelection,valLabelsAntennaSelection},...
    'ValidationFrequency',validationFrequency,...
    'ValidationPatience', 200,...
    'Plots','none',...
    'Shuffle','every-epoch',...
    'OutputFcn',@(info)stopIfAccuracyNotImproving(info,3));

fprintf(2,['Train CNN For Antenna Selection\n'])
%%
warning off parallel:gpu:device:DeviceLibsNeedsRecompiling
try
    gpuArray.eye(2)^2;
catch ME
end
try
    nnet.internal.cnngpu.reluForward(1);
catch ME
end
%%
convnetAntennaSelection = trainNetwork(dataAntennaSelection, labelsAntennaSelection, layersAntennaSelection, optsAntennaSelection);
%% DNN for HB.
layersFRFChainSelection = [imageInputLayer([sizeInputFRFChainSelection(1:3)],'Normalization', 'zerocenter');
%     convolution2dLayer(2^1,2^3);
%     batchNormalizationLayer
%     reluLayer();
%     maxPooling2dLayer([2 2],'Stride',2);
%     convolution2dLayer(2^1,2^3);
%     batchNormalizationLayer
%     reluLayer();
%     convolution2dLayer(2,2^6);
%     batchNormalizationLayer
%     reluLayer();
%     fullyConnectedLayer(2^10);
    fullyConnectedLayer(2^13);
%     reluLayer();
%     dropoutLayer();
    fullyConnectedLayer(2^8);
%     fullyConnectedLayer(2^7);
%     reluLayer();
%     dropoutLayer();
%     fullyConnectedLayer(QFRF);
%     softmaxLayer();
%     classificationLayer();
    fullyConnectedLayer(sizeOutputFRFChainSelection(2))
    regressionLayer()
    ];
optsFRFSelection = trainingOptions('sgdm',...
    'Momentum', 0.9,...
    'InitialLearnRate',0.0005,... % The default value is 0.01.
    'MaxEpochs',5000,...
    'MiniBatchSize',miniBatchSize,... % The default is 128.
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',.9,...
    'LearnRateDropPeriod',30,...
    'L2Regularization',0.0001,... % The default is 0.0001.
    'ExecutionEnvironment', 'auto',...
    'ValidationData',{valDataRFChainSelection,valLabelsFRFChainSelection},...
    'ValidationFrequency',validationFrequency,...
    'ValidationPatience', 200,...
    'Plots','none',...
    'Shuffle','every-epoch',...
    'OutputFcn',@(info)stopIfAccuracyNotImproving(info,3));
fprintf(2,['Train CNN For Frf \n'])
timeCNNRF = tic;
convnetFRFSelection = trainNetwork(dataFRFChainSelection, labelsRFChainSelection, layersFRFChainSelection, optsFRFSelection);
% stopp


layersWRFChainSelection = [layersFRFChainSelection(1:end-2);
    fullyConnectedLayer(sizeOutputWRFChainSelection(2))
    regressionLayer()
    ];
optsWRFSelection = trainingOptions('sgdm',...
    'Momentum', 0.9,...
    'InitialLearnRate',0.0005,... % The default value is 0.01.
    'MaxEpochs',5000,...
    'MiniBatchSize',miniBatchSize,... % The default is 128.
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',.9,...
    'LearnRateDropPeriod',30,...
    'L2Regularization',0.0001,... % The default is 0.0001.
    'ExecutionEnvironment', 'auto',...
    'ValidationData',{valDataWRFChainSelection,valLabelsWRFChainSelection},...
    'ValidationFrequency',validationFrequency,...
    'ValidationPatience', 200,...
    'Plots','none',...
    'Shuffle','every-epoch',...
    'OutputFcn',@(info)stopIfAccuracyNotImproving(info,3));
fprintf(2,['Train CNN For Wrf \n'])
convnetWRFSelection = trainNetwork(dataWRFChainSelection, labelsWRFChainSelection, layersWRFChainSelection, optsWRFSelection);
timeCNNRF = toc(timeCNNRF)
% stopp
% xLabels(:,1) = valLabelsAntennaSelection;
% xLabels(:,2) = classify(convnetAntennaSelection,valDataAntennaSelection);
% accuracySourceVal = [mean(xLabels(:,2) == xLabels(:,1))];
%% Save.
% stopp
% load('scenarioIndex')
% iExp = iExp + 1;
% fileNameSavedNetwork = ['scenarioID_' num2str(fileNameLabelsTrain(12:end)), '_Exp' num2str(iExp)];
% fileNameLabelsTrainNet = varToSave(iLabels,iTrain,iExp,fileNameSavedNetwork,NtRF,Nr,NrRF,Nrs,...
%     bestAntennas,subSet,XAS,XRF,Y,YFRFr,YWRFr,Z,opts,convnetAntennaSelection,convnetFRFSelection,convnetWRFSelection,idx);
% save('scenarioIndex','iLabels','iTrain','iExp')
% fileNameLabelsTrainNet
% %% Run performance test
% beep 
run main03_TestNetwork.m