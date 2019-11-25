function [H,Ar,At,rxang,txang] = generate_H_Ar_At(Ncl,Nscatter,Nray,angspread,lambda,txarray,rxarray,opts)
if opts.fixedUsers == 0
    rng('shuffle');
    txclang = [rand(1,Ncl)*120-60;rand(1,Ncl)*60-30]; %+-60, +-30
    rxclang = [rand(1,Ncl)*120-60;rand(1,Ncl)*60-30];
    txang = zeros(2,Nscatter);
    rxang = zeros(2,Nscatter);
    % compute the rays within each cluster
    for m = 1:Ncl
        txang(:,(m-1)*Nray+(1:Nray)) = randn(2,Nray)*sqrt(angspread)+txclang(:,m);
        rxang(:,(m-1)*Nray+(1:Nray)) = randn(2,Nray)*sqrt(angspread)+rxclang(:,m);
    end
    rng(4096);
else
    rng(4096);
    txclang = [rand(1,Ncl)*120-60;rand(1,Ncl)*60-30];
    rxclang = [rand(1,Ncl)*120-60;rand(1,Ncl)*60-30];
    txang = zeros(2,Nscatter);
    rxang = zeros(2,Nscatter);
    % compute the rays within each cluster
    for m = 1:Ncl
        txang(:,(m-1)*Nray+(1:Nray)) = randn(2,Nray)*sqrt(angspread)+txclang(:,m);
        rxang(:,(m-1)*Nray+(1:Nray)) = randn(2,Nray)*sqrt(angspread)+rxclang(:,m);
    end
end
if opts.fixedChannelGain == 0
    rng('shuffle');
    g = (randn(1,Nscatter)+1i*randn(1,Nscatter))/sqrt(Nscatter); % Channel Gain, 1 x Nscatter.
    rng(4096);
else
    rng(4096);
    g = (randn(1,Nscatter)+1i*randn(1,Nscatter))/sqrt(Nscatter); % Channel Gain, 1 x Nscatter.
end
txpos = getElementPosition(txarray)/lambda;
rxpos = getElementPosition(rxarray)/lambda;
%% quantization.
% for q = 1:length(opts.qBits)
%     qb = opts.qBits(q) + 1;
% %     txangQ = double(quantize(fi(txang),numerictype('WordLength',16,'FractionLength',qb),'Nearest', 'Saturate'));
% %     rxangQ = double(quantize(fi(rxang),numerictype('WordLength',16,'FractionLength',qb),'Nearest', 'Saturate'));
%     txangQ = [linspace(-60,60,2^qb);linspace(-30,30,2^qb)];
%     rxangQ = [linspace(-60,60,2^qb);linspace(-30,30,2^qb)];
%     AQ(q).AtQ = steervec(txpos,txangQ);
%     AQ(q).ArQ = steervec(rxpos,rxangQ);
%     AQ(q).txangQ = txangQ;
%     AQ(q).rxangQ = rxangQ;
% end
%%
At = steervec(txpos,txang);
Ar = steervec(rxpos,rxang);

H = scatteringchanmtx(txpos,rxpos,txang,rxang,g);
H = H.';