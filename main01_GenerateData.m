% These MATLAB scripts are prepared by A.M.E for the following paper,
% Ahmet M. Elbir, "CNN-based Precoder and Combiner Design in mmWave MIMO Systems", IEEE Communications Letters, in press.
% please cite the above work if you use this file. For any comments and
% questions please email: ahmetmelbir@gmail.com

clear all;
addpath('./AltMin/Narrowband');
addpath(genpath('./AltMin'));
opts.Nt_param = 16; % number of Tx antennas.
NtRF = 4; % number of PS at Tx.
NrRF = 4; % number of PS at Rx.
Nr = 16; % number of Rx antennas.
Nrs = 16; % number of selected antennas out of Nr.
%%
% opts.qBits = 6;
opts.Nray_param = 4; % number of rays for each user.
opts.Ncl_param = 4; % number of clusters.
opts.angspread = 5; % angle spread of users.
%% Generate Input data for CNN
opts.selectOutputAsPhases = 1;
opts.snr_param = [0]; % SNR dB.
opts.Ns_param = [2]; % number of data streams.
opts.Nreal = 100; % number of realizations.
opts.Nchannels = 10; % number of channel matrices for the input data, choose small values for fast results in training.
opts.fixedUsers = 0;
opts.fixedChannelGain = 0;
% opts.generateNoisyChannels = 1;
opts.noiseLevelHdB = [15 20 25]; % dB.
opts.inputH = 1;
opts.inputRy = 0;
timeGenerate = tic;
[XAS,XRF,Y,YFRFr,YWRFr,bestAntennas,subSet,Z,opts] = generateMIMO(Nr,Nrs,NrRF,NtRF,opts);
timeGenerate = toc(timeGenerate);
% stopp
%% Save the training data
% iLabels = 1;
% iTrain = 1;
% iExp = 1;
% save('scenarioIndex','iLabels','iTrain','iExp')
% load('scenarioIndex') % where the scenario indices are saved. 
% iLabels = iLabels + 1;
% iTrain = iTrain + 1;
% iLabels = str2num(fileNameLabels(12:14)); % run for the fist time to set
% label number.
% fileNameLabelsTrain = varToSave(iLabels,iTrain,NtRF,Nr,NrRF,Nrs,...
%     bestAntennas,subSet,XAS,XRF,Y,YFRFr,YWRFr,Z,opts);
% save('scenarioIndex','iLabels','iTrain','iExp')
% display([ num2str(iLabels), '_' num2str(iTrain),'_' num2str(iExp) ', Ns=' num2str(opts.Ns_param) ', ' num2str(Nr) '/' num2str(Nrs) ])
% beep
% fileNameLabelsTrain
run main02_TrainNetwork.m
