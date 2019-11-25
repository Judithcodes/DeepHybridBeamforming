% These MATLAB scripts are prepared by A.M.E for the following paper,
% Ahmet M. Elbir, "CNN-based Precoder and Combiner Design in mmWave MIMO Systems", IEEE Communications Letters, in press.
% please cite the above work if you use this file. For any comments and
% questions please email: ahmetmelbir@gmail.com

addpath('./AltMin/Narrowband');
addpath(genpath('./AltMin'));
% profile on
%%
SNR_index = -20:5:20; % noise added to the received signal.
SNR_TEST_DATA = 10;%-10:5:30; % synthetic noise added to the test data.
xIndex = SNR_index;
% xIndex = SNR_TEST_DATA;
Ntrial = 10;
%% Select Methods.
t = num2cell(1:11);
[iFullOPT, iFullSOMP, iDASOPT, iRASOPT, iDBest, iDASHB, iRASDHB, iDASAltMin, iRASAltMin, iDASSOMP, iRASSOMP] = deal(t{:});
% tempSet = [iDASHB iDASSOMP iDBest];
if Nr == Nrs %% Full Array Test
    selectedMethods = [iDASOPT, iDBest, iDASHB, iDASAltMin, iDASSOMP];
    LegendsR{iDASOPT} = 'OPT'; MarkerSet{iDASOPT} = 'none';LineStyleSet{iDASOPT} = '-'; % optimum BF. subarray
    LegendsR{iDASHB} = 'HBDL' ;MarkerSet{iDASHB} = 'pentagram';LineStyleSet{iDASHB} = '-'; % DAS + DHB
    LegendsR{iDBest} = 'Best' ;MarkerSet{iDBest} = 'pentagram';LineStyleSet{iDBest} = '-.'; % DAS + DHB
    LegendsR{iDASAltMin} = 'PE Alt-Min';MarkerSet{iDASAltMin} = 'v';LineStyleSet{iDASAltMin} = '-'; % DAS + MOAltMin
    LegendsR{iDASSOMP}  = 'SOMP'; MarkerSet{iDASSOMP}  = '<';LineStyleSet{iDASSOMP}  = '-.';
else % Antenna Selection Test.
    selectedMethods = [iDASOPT, iRASOPT, iDASHB, iRASDHB, iDASAltMin, iRASAltMin, iDASSOMP, iRASSOMP];
    % selectedMethods = [iFullOPT, iFullSOMP, iDASOPT, iRASOPT, iDASHB, iRASDHB, iDASAltMin, iRASAltMin, iDASSOMP, iRASSOMP];
    LegendsR{iFullOPT} = 'Full Array, OPT' ;MarkerSet{iFullOPT} = 'o';LineStyleSet{iFullOPT} = '-';% optimum BF. Full array
    % LegendsR{iFullDHB} = 'Full array, DHB' ;MarkerSet{iFullDHB} = 'pentagram';LineStyleSet{iFullDHB} = '-'; % DAS + DHB
    % LegendsR{iFullAltMin} = 'Subarray, DAS + MO Alt-Min';MarkerSet{iFullAltMin} = 'v';LineStyleSet{iFullAltMin} = '-'; % DAS + MOAltMin
    LegendsR{iFullSOMP} = 'Full Array, SOMP' ; MarkerSet{iFullSOMP} = 'none' ;LineStyleSet{iFullSOMP} = '-' ;% full array + SOMP
    LegendsR{iDASOPT} = 'DAS + OPT'; MarkerSet{iDASOPT} = 'none';LineStyleSet{iDASOPT} = '-'; % optimum BF. subarray
    
    LegendsR{iDASHB} = 'DAS + DHB' ;MarkerSet{iDASHB} = 'pentagram';LineStyleSet{iDASHB} = '-'; % DAS + DHB
    LegendsR{iRASDHB} = 'RAS + DHB' ;MarkerSet{iRASDHB} = '*';LineStyleSet{iRASDHB} = '-.'; % RAS + DHB
    LegendsR{iDASAltMin} = 'DAS + MO Alt-Min';MarkerSet{iDASAltMin} = 'v';LineStyleSet{iDASAltMin} = '-'; % DAS + MOAltMin
    
    LegendsR{iRASAltMin} = 'RAS + MO Alt-Min'; MarkerSet{iRASAltMin} = '^';LineStyleSet{iRASAltMin} = '-'; % RAS + MOAltMin
    LegendsR{iRASSOMP}  = 'RAS + SOMP'; MarkerSet{iRASSOMP}  = '<';LineStyleSet{iRASSOMP}  = '-';
    LegendsR{iDASSOMP}  = 'DAS + SOMP'; MarkerSet{iDASSOMP}  = '<';LineStyleSet{iDASSOMP}  = '-.';
    
    LegendsR{iRASOPT} = 'RAS + OPT'; MarkerSet{iRASOPT} = 'none';LineStyleSet{iRASOPT} = '-.';
end

sA = cell(numel(selectedMethods),1); % selected antennas.
sRF = cell(numel(selectedMethods),1); % selected RF Chains.
%% Test Data.
dataAntennaSelection = XAS;
labelsAntennaSelection = categorical(Y);
Nray = opts.Nray_param;
Ncl = opts.Ncl_param;
Nscatter = Nray*Ncl;
%%
if length(idx) > Ntrial
    trialIndex = idx(1:Ntrial); % Validation data indices.
else
    trialIndex = idx;%(1:Ntrial); % Validation data indices.
end
%% Generate Array Positions.
rng(4096);
c = 3e8;
fc = 28e9;
lambda = c/fc;
Nt = opts.Nt_param(1);
txarray = phased.PartitionedArray(...
    'Array',phased.URA([sqrt(Nt) sqrt(Nt)],lambda/2),...
    'SubarraySelection',ones(NtRF,Nt),'SubarraySteering','Custom');
rxarray = phased.PartitionedArray(...
    'Array',phased.URA([sqrt(Nr) sqrt(Nr)],lambda/2),...
    'SubarraySelection',ones(NrRF,Nr),'SubarraySteering','Custom');
txpos = getElementPosition(txarray)/lambda;
rxpos = getElementPosition(rxarray)/lambda;
%%
Nx = length(xIndex);
Ropt0 = zeros(Nx,Ntrial,3); Ropt = zeros(Nx,numel(selectedMethods));
Rhyb0 = zeros(Nx,Ntrial,3); Rhyb = zeros(Nx,numel(selectedMethods));
for i = 1:length(xIndex)
    %% SNR TEST
        snr = db2pow(xIndex(i));
        snrTEST = SNR_TEST_DATA;
        Ns =  opts.Ns_param;
        Ncl = opts.Ncl_param;
    xLabel = 'SNR, [dB]';
    %% SNRTESTDATA TEST
%     snr = db2pow(SNR_index);
%     snrTEST = SNR_TEST_DATA(i);
%     Ns =  opts.Ns_param;
%     Ncl = opts.Ncl_param;
%     xLabel = 'SNR_{TEST}, [dB]';
    %%
    for iTrial = trialIndex
        %% Generate Data.
        %         H = Z(1,iTrial).H; % No noise.
        H = awgn(Z(1,iTrial).H,snrTEST,'Measured');
        Ar = Z(1,iTrial).Ar;
        At = Z(1,iTrial).At;
        Arb = Z(1,iTrial).Arb;
        Atb = Z(1,iTrial).Atb;
        txang = Z(1,iTrial).txang;
        rxang = Z(1,iTrial).rxang;
        if opts.selectOutputAsPhases == 1
            FrfSelected = Z(1,iTrial).FrfSelected;
            WrfSelected = Z(1,iTrial).WrfSelected;
            FbbSelected = Z(1,iTrial).FbbSelected;
            WbbSelected = Z(1,iTrial).WbbSelected;
        else
            DOASelected = Z(1,iTrial).DOASelected;
            DODSelected = Z(1,iTrial).DODSelected;
        end
        randomSelection = subSet(randi([size(subSet,1)],1),:);
        for iM = selectedMethods
            switch iM
                case iFullOPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% Full Array, baseband beamformer. Optimum.
                    sA{iM} = 1:Nr;
                    %%
                    [Fopt,Wopt] = helperOptimalHybridWeights(H(sA{iM},:),Ns,1/snr);
                    Ropt0(i,iTrial,iM) = helperComputeSpectralEfficiency(H(sA{iM},:),Fopt,Wopt,Ns,snr);
                    Rhyb0(i,iTrial,iM) = Ropt0(i,iTrial,iM);
                case iRASOPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% RAS, baseband beamformer. Optimum.
                    sA{iM} = randomSelection;
                    %%
                    [Fopt,Wopt] = helperOptimalHybridWeights(H(sA{iM},:),Ns,1/snr);
                    Ropt0(i,iTrial,iM) = helperComputeSpectralEfficiency(H(sA{iM},:),Fopt,Wopt,Ns,snr);
                    Rhyb0(i,iTrial,iM) = Ropt0(i,iTrial,iM);
                case iDASOPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% DAS, Deep Antenna Selection.
                    XTest(:,:,1) = abs(H);
                    XTest(:,:,2) = real(H);
                    XTest(:,:,3) = imag(H);
                    YTest = Z(1,iTrial).Y;
                    timeCNNTemp = tic;
                    [YPred,~] = classify(convnetAntennaSelection,XTest); % find the best array index.
                    timeCNN0(iTrial) = toc(timeCNNTemp);
                    YPred = str2num(char(YPred));
                    accuracy0(i,iTrial,1) = mean(YPred == (YTest).'); % compute accuracy.
                    sA{iM} = subSet(YPred,:);
                    %%
                    [Fopt,Wopt] = helperOptimalHybridWeights(H(sA{iM},:),Ns,1/snr);
                    Ropt0(i,iTrial,iM) = helperComputeSpectralEfficiency(H(sA{iM},:),Fopt,Wopt,Ns,snr);
                    Rhyb0(i,iTrial,iM) = Ropt0(i,iTrial,iM);
                case iDBest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    sA{iM} = sA{iDASOPT};
                    %%
                    [Fopt,Wopt] = helperOptimalHybridWeights(H(sA{iM},:),Ns,1/snr);
                    Frf = FrfSelected;
                    Fbb = FbbSelected;
%                     Fbb = (Frf'*Frf)\Frf'*Fopt;
%                     Fbb = sqrt(Ns)*Fbb/norm(Frf*Fbb,'fro');
                    Wrf = WrfSelected;
                    Wbb = WbbSelected;
%                     Wmmse = ((Fopt'*(H(sA{iM},:)'*H(sA{iM},:))*Fopt+1/snr*Ns*eye(Ns))\Fopt'*H(sA{iM},:)')';
%                     Ess = 1/Ns*eye(Ns);
%                     Eyy = H(sA{iM},:)*Fopt*Ess*Fopt'*H(sA{iM},:)'+ 1/snr*eye(Nrs);
%                     Wbb = (Wrf'*Eyy*Wrf)\(Wrf'*Eyy*Wmmse);
%                     Wrf = conj(Wrf);
%                     Wbb = conj(Wbb);
% Wbb = (Wrf'*Wrf)\Wrf'*Wopt;
%                                                 Wbb = sqrt(Ns)*Wbb/norm(Wrf*Wbb,'fro');
%% Test
% Wrf2 = conj(Wrf);
% Wbb2 = (Wbb);
% norm(Wopt - Wrf2*Wbb2)
%%
                    sRF{iM}.Frf = Frf;
                    sRF{iM}.Fbb = Fbb;
                    sRF{iM}.Wrf = Wrf;
                    sRF{iM}.Wbb = Wbb;
                case iDASHB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% DAS + DHB.
                    sA{iM} = sA{iDASOPT};
                    %% DHB, Deep Hybrid Beamformer.
                    [Fopt,Wopt] = helperOptimalHybridWeights(H(sA{iM},:),Ns,1/snr);
                    XRFTest(:,:,1) = abs(H(sA{iM},:));
                    XRFTest(:,:,2) = real(H(sA{iM},:));
                    XRFTest(:,:,3) = imag(H(sA{iM},:)); % input data
                    timeCNNRFTemp = tic;
                    [YFRFPred] = double(predict(convnetFRFSelection,XRFTest)); % estimate precoder
                    [YWRFPred] = double(predict(convnetWRFSelection,XRFTest)); % estimate combiner
                    timeCNNRF0(iTrial) = toc(timeCNNRFTemp);
                    if opts.selectOutputAsPhases == 1
                        % find baseband beamformers
                        [Frf,Fbb] = findFrfFbb(H(sA{iM},:),Ns,NtRF,exp(1i*reshape(YFRFPred,[Nt,NtRF])));
                        [Wrf,Wbb] = findWrfWbb(H(sA{iM},:),Ns,NrRF,exp(-1i*reshape(YWRFPred,[Nrs,NrRF])),1/snr);
%                         Frf = exp(1i*reshape(YFRFPred,[Nt,NtRF]));
%                         Fbb = (Frf'*Frf)\Frf'*Fopt;
%                         Fbb = sqrt(Ns)*Fbb/norm(Frf*Fbb,'fro');
%                         Wrf = exp(1i*reshape(YWRFPred,[Nrs,NrRF]));
%                         Wrf = conj(Wrf);
%                         Wmmse = ((Fopt'*(H(sA{iM},:)'*H(sA{iM},:))*Fopt+1/snr*Ns*eye(Ns))\Fopt'*H(sA{iM},:)')';
%                         Ess = 1/Ns*eye(Ns);
%                         Eyy = H(sA{iM},:)*Fopt*Ess*Fopt'*H(sA{iM},:)'+ 1/snr*eye(Nrs);
%                         Wbb = (Wrf'*Eyy*Wrf)\(Wrf'*Eyy*Wmmse);
%                                                 Wbb = (Wrf'*Wrf)\Wrf'*Wopt;
%                                                 Wbb = sqrt(Ns)*Wbb/norm(Wrf*Wbb,'fro');
%                         Wrf = conj(Wrf);
%                         Wbb = conj(Wbb);
                    else
                        % DOA based approach/ not used.
                        txangEst = reshape(YFRFPred, [2, NtRF])*180/pi;
                        rxangEst = reshape(YWRFPred, [2, NrRF])*180/pi;
                        errorDOA(i,iTrial,iM) = rms(rms(DOASelected - txangEst)); % compute DOA error
                        errorDOD(i,iTrial,iM) = rms(rms(DODSelected - rxangEst)); % compute DOA error
                        AtEst = steervec(txpos,txangEst);
                        ArEst = steervec(rxpos,rxangEst);
                        [Frf,Fbb] = findFrfFbb(H(sA{iM},:),Ns,NtRF,AtEst);
                        [Wrf,Wbb] = findWrfWbb(H(sA{iM},:),Ns,NrRF,ArEst(sA{iM},:),1/snr);
                    end
                    % compute LS error.
                    error_F_Prediction0(i,iTrial,iM) = rms(rms ( sRF{iDBest}.Frf*sRF{iDBest}.Fbb - Frf*Fbb))/numel(Fopt);
                    error_W_Prediction0(i,iTrial,iM) = rms(rms ( sRF{iDBest}.Wrf*sRF{iDBest}.Wbb - Wrf*Wbb))/numel(Fopt);
                    sRF{iM}.Frf = Frf;
                    sRF{iM}.Fbb = Fbb;
                    sRF{iM}.Wrf = Wrf;
                    sRF{iM}.Wbb = Wbb;
                case iRASDHB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% RAS + DHB.
                    sA{iM} = randomSelection;
                    %% DHB, Deep Hybrid Beamformer.
                    if opts.selectOutputAsPhases == 1
                        Frf = exp(1i*reshape(YFRFPred,[Nt,NtRF]));
                        Fbb = (Frf'*Frf)\Frf'*Fopt;
                        Fbb = sqrt(Ns)*Fbb/norm(Frf*Fbb,'fro');
                        Wrf = exp(1i*reshape(YWRFPred,[Nrs,NrRF]));
                        Wmmse = ((Fopt'*(H(sA{iM},:)'*H(sA{iM},:))*Fopt+1/snr*Ns*eye(Ns))\Fopt'*H(sA{iM},:)')';
                        Ess = 1/Ns*eye(Ns);
                        Eyy = H(sA{iM},:)*Fopt*Ess*Fopt'*H(sA{iM},:)'+ 1/snr*eye(Nrs);
                        Wbb = (Wrf'*Eyy*Wrf)\(Wrf'*Eyy*Wmmse);
                        Wrf = conj(Wrf);
                        Wbb = conj(Wbb);
                        %                         Wbb = (Wrf'*Wrf)\Wrf'*Wopt;
                        %                         Wbb = sqrt(Ns)*Wbb/norm(Wrf*Wbb,'fro');
                    else
                        txangEst = reshape(YFRFPred, [2, NtRF])*180/pi;
                        rxangEst = reshape(YWRFPred, [2, NrRF])*180/pi;
                        errorDOA(i,iTrial,iM) = rms(rms(DOASelected - txangEst)); % compute DOA error
                        errorDOD(i,iTrial,iM) = rms(rms(DODSelected - rxangEst)); % compute DOA error
                        AtEst = steervec(txpos,txangEst);
                        ArEst = steervec(rxpos,rxangEst);
                        [Frf,Fbb] = findFrfFbb(H(sA{iM},:),Ns,NtRF,AtEst);
                        [Wrf,Wbb] = findWrfWbb(H(sA{iM},:),Ns,NrRF,ArEst(sA{iM},:),1/snr);
                    end
                    sRF{iM}.Frf = Frf;
                    sRF{iM}.Fbb = Fbb;
                    sRF{iM}.Wrf = Wrf;
                    sRF{iM}.Wbb = Wbb;
                case iDASAltMin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% DAS + Alt-Min.
                    sA{iM} = sA{iDASHB};
                    %%
                    [Fopt,Wopt] = helperOptimalHybridWeights(H(sA{iM},:),Ns,1/snr);
                    [Fmo,sRF{iM}.Frf,sRF{iM}.Fbb] = MO_AltMin_F( Fopt, NtRF, 1e-1);
                    [Wmo,sRF{iM}.Wrf,sRF{iM}.Wbb] = MO_AltMin_W( Wopt,Fopt,H(sA{iM},:), NrRF, 1/snr, 1e-1 );
                case iRASAltMin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% RAS + Alt-Min.
                    sA{iM} = randomSelection;
                    %%
                    [Fopt,Wopt] = helperOptimalHybridWeights(H(sA{iM},:),Ns,1/snr);
                    [Fmo,sRF{iM}.Frf,sRF{iM}.Fbb] = MO_AltMin_F( Fopt, NtRF, 1e-1);
                    [Wmo,sRF{iM}.Wrf,sRF{iM}.Wbb] = MO_AltMin_W( Wopt,Fopt,H(sA{iM},:), NrRF, 1/snr, 1e-1 );
                case iDASSOMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% DAS + SOMP
                    sA{iM} = sA{iDASHB};
                    %% SOMP
                    [sRF{iM}.Fbb,sRF{iM}.Frf,sRF{iM}.Wbb,sRF{iM}.Wrf]...
                        = helperOMPHybridWeights(H(sA{iM},:),NtRF,NrRF,Ns,At,Ar(sA{iM},:),1/snr);
                case iRASSOMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% RAS + SOMP
                    sA{iM} = randomSelection;
                    %% SOMP
                    [sRF{iM}.Fbb,sRF{iM}.Frf,sRF{iM}.Wbb,sRF{iM}.Wrf]...
                        = helperOMPHybridWeights(H(sA{iM},:),NtRF,NrRF,Ns,At,Ar(sA{iM},:),1/snr);
            end
            %% Calculate sum-rate
            if iM ~= iFullOPT && iM ~= iRASOPT  && iM ~= iDASOPT
                Ropt0(i,iTrial,iM) = helperComputeSpectralEfficiency(H(sA{iM},:),Fopt,Wopt,Ns,snr);
                Rhyb0(i,iTrial,iM) = helperComputeSpectralEfficiency(H(sA{iM},:),...
                    sRF{iM}.Frf*sRF{iM}.Fbb,sRF{iM}.Wrf*sRF{iM}.Wbb,Ns,snr);
                errorF0(i,iTrial,iM) = norm(Fopt - sRF{iM}.Frf*sRF{iM}.Fbb)/numel(Fopt);
                errorW0(i,iTrial,iM) = norm(Wopt - (sRF{iM}.Wrf*sRF{iM}.Wbb))/numel(Wopt);
            end
            
        end
    end
    accuracy(i,1) = rms(accuracy0(i,trialIndex,1));
    %     accuracyRF(i,1) = rms(accuracyRF0(i,trialIndex,1));
    for iM = selectedMethods
        Ropt(i,iM) = rms(Ropt0(i,trialIndex,iM));
        Rhyb(i,iM) = rms(Rhyb0(i,trialIndex,iM));
        errorF(i,iM) = rms(errorF0(i,trialIndex,iM));
        errorW(i,iM) = rms(errorW0(i,trialIndex,iM));
        if iM == iDASHB
            error_F_Prediction(i,iM) = rms(error_F_Prediction0(i,trialIndex,iM));
            error_W_Prediction(i,iM) = rms(error_W_Prediction0(i,trialIndex,iM));
        end
    end
    i
end
errorF(:,selectedMethods)
errorW(:,selectedMethods)
% Legends{1} = ['CNN_{AS}, N_{RS}=' num2str(Nrs)];
% Legends{2} = ['CNN_{RF}, N_{RS}=' num2str(Nrs)];

% | 'o' | '*' | '.' | 'x' | 'square' | 'diamond' | 'v' | '^' |
% '>' | '<' | 'pentagram' | 'hexagram' | 'none'.

% plotIndex = selectedMethods;


figure(101)
% subplot(211)
for k = selectedMethods
    plotR(k) = plot(xIndex,Rhyb(:,k));
    %         plotR(k).Color = ColorSet(k,:);
    plotR(k).Marker = MarkerSet{k};
    plotR(k).LineStyle = LineStyleSet{k};
    plotR(k).LineWidth = 1;
    hold on
end
hold off
% plot(xIndex,[Ropt])
legend(LegendsR{selectedMethods},'Location','Best')
xlabel(xLabel)
ylabel('Spectral Efficiency [bits/s/Hz]')

% magnifyOnFigure;


% figure(102)
% plot(xIndex,error_F_Prediction(:,iDASHB))
% hold on
% for k = selectedMethods
%     plotR(k) = plot(xIndex,errorF(:,k));
%     %         plotR(k).Color = ColorSet(k,:);
%     plotR(k).Marker = MarkerSet{k};
%     plotR(k).LineStyle = LineStyleSet{k};
%     plotR(k).LineWidth = 1;
% end
% plot(xIndex,error_W_Prediction1(:,iDASHB))
% for k = selectedMethods
%     plotR(k) = plot(xIndex,errorW1(:,k));
%     %         plotR(k).Color = ColorSet(k,:);
%     plotR(k).Marker = MarkerSet{k};
%     plotR(k).LineStyle = LineStyleSet{k};
%     plotR(k).LineWidth = 1;
%     hold on
% end
% hold off
% leg1 = legend('Prediction Error for $\textbf{F}_{RF}\textbf{F}_{BB}$',...
%     LegendsR{selectedMethods},'Prediction Error for $\textbf{W}_{RF}\textbf{W}_{BB}$',LegendsR{selectedMethods});
% set(leg1,'Interpreter','latex');
% xlabel(xLabel)

display(['Figure Plotted: Hybrid Beamforming, N_T=' num2str(Nt) ', N_R=' num2str(Nr) '/N_{RS}=' num2str(Nrs) ', N_S=' num2str(Ns)])
% suggestedFileNameSNRTest = ['SNR_Test' , '_Nt' num2str(Nt) '_Nr' num2str(Nr) ...
%     '_Nrs' num2str(Nrs) '_Ns' num2str(Ns) '_' fileNameLabelsTrainNet ]
beep