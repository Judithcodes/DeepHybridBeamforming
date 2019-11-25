function [XAS,XRF,Y,YFRFr,YWRFr,bestAntennas,subSetA,Z,opts] = generateMIMO(Nr,Nrs,NrRF,NtRF,opts)
[Q, subSetA] =  calculateNumberOfSubsets(Nr,Nrs);
Nray = opts.Nray_param;
Ncl = opts.Ncl_param;
Nscatter = Nray*Ncl;

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
%% Generate User DOA/DOD
angspread = opts.angspread;
snr = db2pow(opts.snr_param);
noisevar = 1/snr;
% Nsnr = numel(opts.snr_param);
NNs = numel(opts.Ns_param);
NNcl = numel(opts.Ncl_param);
Ncl = opts.Ncl_param;
Nray = opts.Nray_param;
Nch = opts.Nchannels;
Nreal = opts.Nreal;
% bestClass = zeros(Nreal,1);
% bestAntennas = zeros(Nreal,Nrs);
Nb = NNcl*NNs*Nch;
N = Nreal*Nb;
XAS = zeros(Nr,Nt,3,N);
XRF = zeros(Nrs,Nt,3,N);
Y = zeros(1,N);
Yb = zeros(1,N);
% YRF = zeros(1,N);
if opts.selectOutputAsPhases == 1
    YFRFr = zeros(N,Nt*NtRF);
    YWRFr = zeros(N,Nrs*NrRF);
else
    YFRFr = zeros(N,2*NtRF);
    YWRFr = zeros(N,2*NrRF);
end
Z = repmat(struct('H',zeros(Nt,Nr)),1,N );
R = zeros(Q,1);
% Ropt  = zeros(Q,1);
j = 1;
for nch = 1:Nch
    for ncl = 1:NNcl
        Ncl = opts.Ncl_param(ncl);
        Nscatter = Nray*Ncl;
        for ns = 1:NNs
            Ns = opts.Ns_param(ns);
            for nr = 1:Nreal
                if length(opts.Ns_param) > 1
                    nn = ns;
                elseif length(opts.Ncl_param) > 1
                    nn = ncl;
                else
                    nn = j;
                end
                %%
                if nn < round(Nreal/3)
                    snrH = opts.noiseLevelHdB(1);
                elseif nn > round(Nreal/3) && nn < round(2*Nreal/3)
                    snrH = opts.noiseLevelHdB(2);
                elseif nn > round(2*Nreal/3)
                    snrH = opts.noiseLevelHdB(3);
                end
                %%
                if nr == 1 % first realization.
                    % compute randomly placed scatterer clusters
                    [H,Ar,At,rxang,txang] = generate_H_Ar_At(Ncl,Nscatter,Nray,angspread,lambda,txarray,rxarray,opts);
                    
                    jFirst = j; % first.
                    H = awgn(H,snrH,'Measured');
                    Z(1,nn).H = H;
                    Z(1,nn).At = At;
                    Z(1,nn).Ar = Ar;
                    Z(1,nn).rxang = rxang;
                    Z(1,nn).txang = txang;
                    %% Antenna Selection.
%                     timeAS = tic;
                    for (qA = 1:Q)
                        R(qA,1) = helperComputeSpectralEfficiencyAS(H(subSetA(qA,:),:),Ns,snr);
                    end
%                     timeAS = toc(timeAS);
                    [~,qAb] = max(R(:,1));
                    %% Best antennas.
                    %                     qAb = 1; % no antenna selection.
                    bestAntennas = subSetA(qAb,:);
                    %% RF Chain Selection.
                    % consider MO Alt-Min results in the feasible set for simplicity. for
                    % a large set they will converge eventually. 
                    [FoptMO,WoptMO] = helperOptimalHybridWeights(H,Ns,1/snr);
                    [~,FrfMO, ~] = MO_AltMin_F( FoptMO, NtRF, 1e-3);
                    [~,WrfMO, ~] = MO_AltMin_W( WoptMO, FoptMO, H, NrRF, 1/snr, 1e-3 );
                    
                    Atb = [At FrfMO];
                    Arb = [Ar WrfMO];
                    % or not.
%                                         Atb = [At ];
%                                         Arb = [Ar ];
                    %%
                    [QF, subSetF] =  calculateNumberOfSubsets(size(Atb,2),NtRF);
                    RhybF = zeros(QF,1);
                    RhybW = zeros(QF,1);
                    Frf = zeros(Nt,NtRF,QF);
                    Fbb = zeros(NtRF,Ns,QF);
                    Wrf = zeros(Nrs,NrRF,QF);
                    Wbb = zeros(NrRF,Ns,QF);
                    
                    Z(1,nn).Atb = Atb;
                    Z(1,nn).Arb = Arb;
                    %% RF precoder design.
                    [Fopt,Wopt] = helperOptimalHybridWeights(H(bestAntennas,:),Ns,1/snr);
                    for qF = 1:QF
                        [Frf(:,:,qF),Fbb(:,:,qF)] = findFrfFbb(H(bestAntennas,:),Ns,NtRF,Atb(:,subSetF(qF,:)));
                        RhybF(qF,1) = helperComputeSpectralEfficiency(H(bestAntennas,:),Frf(:,:,qF)*Fbb(:,:,qF),Wopt,Ns,snr);
                    end
                    [~,qFb] = max(RhybF(:,1)); % best RF precoder.
                    %% RF combiner design.
                    %                     Fbest = Frf(:,:,qFb)*Fbb(:,:,qFb); % worse results.
                    Fbest = Fopt;
                    for qW = 1:QF
                        [Wrf(:,:,qW),Wbb(:,:,qW)] = findWrfWbb(H(bestAntennas,:),Ns,NrRF,Arb(bestAntennas,subSetF(qW,:)),1/snr);
                        RhybW(qW,1) = helperComputeSpectralEfficiency(H(bestAntennas,:),Fbest,Wrf(:,:,qW)*Wbb(:,:,qW),Ns,snr);
                    end
                    [~,qWb] = max(RhybW(:,1)); % best RF combiner.
                    %%
%                     selectedSteeringVectors = [subSetF(qFb,:);subSetF(qWb,:)]
                    %% Select phases for output
                    if opts.selectOutputAsPhases == 1
                        Z(1,nn).FrfSelected = Frf(:,:,qFb);
                        Z(1,nn).WrfSelected = Wrf(:,:,qWb);
                        Z(1,nn).FbbSelected = Fbb(:,:,qFb);
                        Z(1,nn).WbbSelected = Wbb(:,:,qWb);
                    else %% Select DOA/DOD for output
                        DOASelected = txang(:,subSetF(qFb,:));
                        DODSelected = rxang(:,subSetF(qWb,:));
                        Z(1,nn).rxang = rxang;
                        Z(1,nn).txang = txang;
                        Z(1,nn).DOASelected = DOASelected;
                        Z(1,nn).DODSelected = DODSelected;
                    end
                else % other realizations.
                    H = awgn(Z(1,jFirst).H,snrH,'Measured');
                    Z(1,nn).H = H;
                    Z(1,nn).At = Z(1,jFirst).At;
                    Z(1,nn).Ar = Z(1,jFirst).Ar;
                    Z(1,nn).Atb = Z(1,jFirst).Atb;
                    Z(1,nn).Arb = Z(1,jFirst).Arb;
                    if opts.selectOutputAsPhases == 1 %% Select phases for output
                        Z(1,nn).FrfSelected = Z(1,jFirst).FrfSelected;%Frf(:,:,qFb);
                        Z(1,nn).WrfSelected = Z(1,jFirst).WrfSelected;%Wrf(:,:,qWb);
                        Z(1,nn).FbbSelected = Z(1,jFirst).FbbSelected;%Frf(:,:,qFb);
                        Z(1,nn).WbbSelected = Z(1,jFirst).WbbSelected;%Wrf(:,:,qWb);
%                         Z(1,nn).FrfSelected = Frf(:,:,qFb);
%                         Z(1,nn).WrfSelected = Wrf(:,:,qWb);
                    else %% Select DOA/DOD for output
                        DOASelected = txang(:,subSetF(qFb,:));
                        DODSelected = rxang(:,subSetF(qWb,:));
                        Z(1,nn).rxang = rxang;
                        Z(1,nn).txang = txang;
                        Z(1,nn).DOASelected = DOASelected;
                        Z(1,nn).DODSelected = DODSelected;
                    end
                end
                %% output of the network. classification
                Y(1,j) = (qAb);
%                 YRF(1,j) = qFb;
                Yb(1,j) = Y(1,jFirst);
                Z(1,nn).Y = Yb(1,j);
%                 Z(1,nn).YRF = [qFb qWb];
                if opts.selectOutputAsPhases == 1
                    YFRFr(j,:) = angle(vec(Z(1,nn).FrfSelected));
                    YWRFr(j,:) = angle(vec(Z(1,nn).WrfSelected));
%                     if isnan(YFRFr(j,:)) == 1
%                         pause
%                     end
                else % DOA and DOD of users.
                    YFRFr(j,:) = pi/180*reshape(DOASelected,[2*NtRF,1]).';
                    YWRFr(j,:) = pi/180*reshape(DODSelected,[2*NrRF,1]).';
                end
                XAS(:,:,1,j) = abs(H);
                XAS(:,:,2,j) = real(H);
                XAS(:,:,3,j) = imag(H);
                XRF(:,:,1,j) = abs(H(bestAntennas,:));
                XRF(:,:,2,j) = real(H(bestAntennas,:));
                XRF(:,:,3,j) = imag(H(bestAntennas,:));
                
                
                j = j + 1;
                %                 nch
            end
        end
    end
    nch
end
Y = Yb;
% Yu = kron(Yb,ones(1,Nb));
% Y = Yu;