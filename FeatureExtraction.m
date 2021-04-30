%% Feature Extraction
totalTrials = 75;

%% SUB PRE
signal = sub1PRE_DATA;

% preallocate feature outputs for speed
[sub1_alphaPower_PRE] = zeros(1,32,totalTrials);
[sub1_alphaMean_PRE] = zeros(1,32,totalTrials);
[sub1_alphaSig_PRE] = zeros(length(signal), 32, totalTrials);

[sub1_betaPower_PRE] = zeros(1,32,totalTrials);
[sub1_betaMean_PRE] = zeros(1,32,totalTrials);
[sub1_betaSig_PRE] = zeros(length(signal), 32, totalTrials);



 % DWT, 6 level wavelet db4
for trial = 1:totalTrials
    for channel = 1:32
        [alphaPower,alphaMean,betaPower,betaMean,alphaSig,betaSig] = DWT(signal(:,channel,trial)); 
        sub1_alphaPower_PRE(1,channel,trial) = alphaPower;
        sub1_alphaMean_PRE(1,channel,trial) = alphaMean;

        sub1_betaPower_PRE(1,channel,trial) = betaPower;
        sub1_betaMean_PRE(1,channel,trial) = betaMean;

        
        sub1_alphaSig_PRE(:,channel,trial) = alphaSig;
        sub1_betaSig_PRE(:,channel,trial) = betaSig;
    end
end

%% SUB POST
signal = sub1POST_DATA;

% preallocate feature outputs for speed
[sub1_alphaPower_POST] = zeros(1,32,totalTrials);
[sub1_alphaMean_POST] = zeros(1,32,totalTrials);
[sub1_alphaSig_POST] = zeros(length(signal), 32, totalTrials);

[sub1_betaPower_POST] = zeros(1,32,totalTrials);
[sub1_betaMean_POST] = zeros(1,32,totalTrials);
[sub1_betaSig_POST] = zeros(length(signal), 32, totalTrials);


 % DWT, 6 level wavelet db4
for trial = 1:totalTrials
    for channel = 1:32
        [alphaPower,alphaMean,betaPower,betaMean, alphaSig,betaSig] = DWT(signal(:,channel,trial)); 
        sub1_alphaPower_POST(1,channel,trial) = alphaPower;
        sub1_alphaMean_POST(1,channel,trial) = alphaMean;
        sub1_betaPower_POST(1,channel,trial) = betaPower;
        sub1_betaMean_POST(1,channel,trial) = betaMean;
        
        sub1_alphaSig_POST(:,channel,trial) = alphaSig;
        sub1_betaSig_POST(:,channel,trial) = betaSig;
    end
end


%% Power analysis
MIchannels = [6,9,10,11,12,15,16,17,21,22,25,27];

%% PRE
% flx
% calculates the mean of MI channels for each trial
extAlphaGrandAVGPower = mean(sub1_alphaPower_PRE(1,MIchannels,1:25));

extBetaGrandAVGPower = mean(sub1_betaPower_PRE(1,MIchannels,1:25));


% ext
flxAlphaGrandAVGPower = mean(sub1_alphaPower_PRE(1,MIchannels,26:50));

flxBetaGrandAVGPower = mean(sub1_betaPower_PRE(1,MIchannels,26:50));


% rest
restAlphaGrandAVGPower = mean(sub1_alphaPower_PRE(1,MIchannels,51:75));

restBetaGrandAVGPower = mean(sub1_betaPower_PRE(1,MIchannels,51:75));


%% ERD - percent change of power between rest and action
%% EXT
% Alpha
rest = mean(restAlphaGrandAVGPower);
action = mean(extAlphaGrandAVGPower);
erdAlpha_PRE_ext = 100 * (action - rest)/abs(rest);
Alpha_PRE_ext = action;
Alpha_PRE_rest = rest;
% Beta
rest = mean(restBetaGrandAVGPower);
action = mean(extBetaGrandAVGPower);
erdBeta_PRE_ext = 100 * (action - rest)/abs(rest);
Beta_PRE_ext = action;
Beta_PRE_rest = rest;

%% FLX 
% Alpha
rest = mean(restAlphaGrandAVGPower);
action = mean(flxAlphaGrandAVGPower);
erdAlpha_PRE_flx = 100 * (action - rest)/abs(rest);
Alpha_PRE_flx = action;
% Beta
rest = mean(restBetaGrandAVGPower);
action = mean(flxBetaGrandAVGPower);
erdBeta_PRE_flx = 100 * (action - rest)/abs(rest);
Beta_PRE_flx = action;




%% POST

% ext
% calculates the mean of 32 channels for each trial
mean(sub1_alphaSigEXT_PRE,[2 3])
extAlphaGrandAVGPower = mean(sub1_alphaPower_POST(1,MIchannels,1:25));

extBetaGrandAVGPower = mean(sub1_betaPower_POST(1,MIchannels,1:25));


% flx
flxAlphaGrandAVGPower = mean(sub1_alphaPower_POST(1,MIchannels,26:50));

flxBetaGrandAVGPower = mean(sub1_betaPower_POST(1,MIchannels,26:50));


% rest
restAlphaGrandAVGPower = mean(sub1_alphaPower_POST(1,MIchannels,51:75));

restBetaGrandAVGPower = mean(sub1_betaPower_POST(1,MIchannels,51:75));


%% ERD - percent change of power between rest and action
%% EXT
% Alpha
rest = mean(restAlphaGrandAVGPower);
action = mean(extAlphaGrandAVGPower);
erdAlpha_POST_ext = 100 * (action - rest)/abs(rest);
Alpha_POST_ext = action;
Alpha_POST_rest = rest;
% Beta
rest = mean(restBetaGrandAVGPower);
action = mean(extBetaGrandAVGPower);
erdBeta_POST_ext = 100 * (action - rest)/abs(rest);
Beta_POST_ext = action;
Beta_POST_rest = rest;
%% FLX 
% Alpha
rest = mean(restAlphaGrandAVGPower);
action = mean(flxAlphaGrandAVGPower);
erdAlpha_POST_flx = 100 * (action - rest)/abs(rest);
Alpha_POST_flx = action;
% Beta
rest = mean(restBetaGrandAVGPower);
action = mean(flxBetaGrandAVGPower);
erdBeta_POST_flx = 100 * (action - rest)/abs(rest);
Beta_POST_flx = action;

%% Visualize Power changes

figure;
X = categorical({'PRE TESS','POST TESS'});
X = reordercats(X,{'PRE TESS', 'POST TESS'});
y = [Alpha_PRE_ext Beta_PRE_ext; Alpha_POST_ext Beta_POST_ext];
h = bar(X,y)
set(h, {'DisplayName'}, {'Alpha','Beta'}')
legend()
title("Average power of alpha and beta bands - EXT");

figure;
X = categorical({'PRE TESS','POST TESS'});
X = reordercats(X,{'PRE TESS', 'POST TESS'});
y = [Alpha_PRE_flx Beta_PRE_flx; Alpha_POST_flx Beta_POST_flx];
h = bar(X,y)
set(h, {'DisplayName'}, {'Alpha','Beta'}')
legend()
title("Average power of alpha and beta bands - FLX");

figure;
X = categorical({'PRE TESS','POST TESS'});
X = reordercats(X,{'PRE TESS', 'POST TESS'});
y = [Alpha_PRE_rest Beta_PRE_rest; Alpha_POST_rest Beta_POST_rest];
h = bar(X,y)
set(h, {'DisplayName'}, {'Alpha','Beta'}')
legend()
title("Average power of alpha and beta bands - REST");



figure;
X = categorical({'PRE TESS','POST TESS'});
X = reordercats(X,{'PRE TESS', 'POST TESS'});
y = [erdAlpha_PRE_ext erdBeta_PRE_ext; erdAlpha_POST_ext erdBeta_POST_ext];
h = bar(X,y)
set(h, {'DisplayName'}, {'Alpha','Beta'}')
legend()
title("Percent change in average power of alpha and beta bands between ext and rest");

figure;
X = categorical({'PRE TESS','POST TESS'});
X = reordercats(X,{'PRE TESS', 'POST TESS'});
y = [erdAlpha_PRE_flx erdBeta_PRE_flx; erdAlpha_POST_flx erdBeta_POST_flx];
h = bar(X,y)
set(h, {'DisplayName'}, {'Alpha','Beta'}')
legend()
title("Percent change in average power of alpha and beta bands between flx and rest");

%% Visualize Means

% Channel selection
% FZ, FC5, FC1, FC2, FC6, C3, Cz, C4, CP1, CP2, P3, P4
% MIchannels = find(channels == 'FZ' || channels == 'FC5' || channels == 'FC1' ...
%     || channels == 'FC2' || channels == 'FC6' || channels == 'C3' ...
%     || channels == 'Cz' || channels == 'C4' || channels == 'CP1' ...
%     || channels == 'CP2' || channels == 'P3' || channels == 'P4');


%% EXT
figure();
for trial = 1:25
    for channel = MIchannels
        scatter(trial,sub1_betaMean_PRE(1,channel,trial),'filled')
        hold on;
    end
end
title("BETA means PRE - EXT")

figure();
for trial = 1:25
    for channel = MIchannels
        scatter(trial,sub1_betaMean_POST(1,channel,trial),'filled')
        hold on;
    end
end
title("BETA means POST - EXT")

figure();
for trial = 1:25
    for channel = MIchannels
        scatter(trial,sub1_alphaMean_PRE(1,channel,trial),'filled')
        hold on;
    end
end
title("ALPHA means PRE - EXT")

figure();
for trial = 1:25
    for channel = MIchannels
        scatter(trial,sub1_alphaMean_POST(1,channel,trial),'filled')
        hold on;
    end
end
title("ALPHA means POST - EXT")

%% Visualize ERD 
% Only looking at channels specifically focused for IM


% PRE

% ext 
sub1_alphaSigEXT_PRE = sub1_alphaSig_PRE(:,MIchannels,1:25);
sub1_GrandAlphaExt_PRE = mean(sub1_alphaSigEXT_PRE,[2 3]);

sub1_betaSigEXT_PRE = sub1_betaSig_PRE(:,MIchannels,1:25);
sub1_GrandBetaExt_PRE = mean(sub1_betaSigEXT_PRE,[2 3]);

% flx 
sub1_alphaSigFLX_PRE = sub1_alphaSig_PRE(:,MIchannels,26:50);
sub1_GrandAlphaFlx_PRE = mean(sub1_alphaSigFLX_PRE,[2 3]);

sub1_betaSigFLX_PRE = sub1_betaSig_PRE(:,MIchannels,26:50);
sub1_GrandBetaFlx_PRE = mean(sub1_betaSigFLX_PRE,[2 3]);

% rest
sub1_alphaSigREST_PRE = sub1_alphaSig_PRE(:,MIchannels,51:75);
sub1_GrandAlphaRest_PRE = mean(sub1_alphaSigREST_PRE,[2 3]);

sub1_betaSigREST_PRE = sub1_betaSig_PRE(:,MIchannels,51:75);
sub1_GrandBetaRest_PRE = mean(sub1_betaSigREST_PRE,[2 3]);



% change in powers

figure();

ERD_alpha_PRE = zeros(length(length(sub1_GrandAlphaRest_PRE)));

for sample = 1:length(sub1_GrandAlphaRest_PRE)

    sampleAction = sub1_GrandAlphaExt_PRE(sample);
    instPowerAction = abs(sampleAction).^2;
    
    sampleRest = sub1_GrandAlphaRest_PRE(sample);
    instPowerRest = abs(sampleRest).^2;
    
    deltaPOW = 100 * (instPowerAction - instPowerRest)/abs(instPowerRest);
    ERD_alpha_PRE(sample) = deltaPOW;  
    
end

ERD_alpha_PRE = filloutliers((ERD_alpha_PRE),'nearest','mean');

plot( (1:length(ERD_alpha_PRE))./fs, ERD_alpha_PRE)

hold on;

ERD_beta_PRE = zeros(length(length(sub1_GrandBetaRest_PRE)));

for sample = 1:length(sub1_GrandBetaRest_PRE)

    sampleAction = sub1_GrandBetaExt_PRE(sample);
    instPowerAction = abs(sampleAction) .^2;
    
    sampleRest = sub1_GrandBetaRest_PRE(sample);
    instPowerRest = abs(sampleRest).^2;
    
    deltaPOW = 100 * (instPowerAction - instPowerRest)/abs(instPowerRest);
    ERD_beta_PRE(sample) = deltaPOW;  
end

ERD_beta_PRE = filloutliers((ERD_beta_PRE),'nearest','mean');

plot( (1:length(ERD_beta_PRE))./fs, ERD_beta_PRE)
% ylim([-10 100])

%%
figure;
plot( (1:length(sub1_GrandBetaExt_PRE))./fs, sub1_GrandBetaExt_PRE, "blue")
hold on;
plot( (1:length(sub1_GrandBetaRest_PRE))./fs, sub1_GrandBetaRest_PRE, "red")
hold on;



%% POST

% ext 
sub1_alphaSigEXT_POST = sub1_alphaSig_POST(:,MIchannels,1:25);
sub1_GrandAlphaExt_POST = mean(sub1_alphaSigEXT_POST,[2 3]);

sub1_betaSigEXT_POST = sub1_betaSig_POST(:,MIchannels,1:25);
sub1_GrandBetaExt_POST = mean(sub1_betaSigEXT_POST,[2 3]);

% flx 
sub1_alphaSigFLX_POST = sub1_alphaSig_POST(:,MIchannels,26:50);
sub1_GrandAlphaFlx_POST = mean(sub1_alphaSigFLX_POST,[2 3]);

sub1_betaSigFLX_POST = sub1_betaSig_POST(:,MIchannels,26:50);
sub1_GrandBetaFlx_POST = mean(sub1_betaSigFLX_POST,[2 3]);

% rest
sub1_alphaSigREST_POST = sub1_alphaSig_POST(:,MIchannels,51:75);
sub1_GrandAlphaRest_POST = mean(sub1_alphaSigREST_POST,[2 3]);

sub1_betaSigREST_POST = sub1_betaSig_POST(:,MIchannels,51:75);
sub1_GrandBetaRest_POST = mean(sub1_betaSigREST_POST,[2 3]);


figure;
plot( (1:length(sub1_GrandBetaExt_POST))./fs, sub1_GrandBetaExt_POST, "blue")
hold on;
plot( (1:length(sub1_GrandBetaRest_POST))./fs, sub1_GrandBetaRest_POST, "red")
hold on;


%% SNR


% EXT 
% PRE
Alpha_SNR_PRE_EXT = snr(sub1_GrandAlphaExt_PRE,sub1_GrandAlphaRest_PRE)
Beta_SNR_PRE_EXT = snr(sub1_GrandBetaExt_PRE,sub1_GrandBetaRest_PRE)
% POST
Alpha_SNR_POST_EXT = snr(sub1_GrandAlphaExt_POST,sub1_GrandAlphaRest_POST)
Beta_SNR_POST_EXT = snr(sub1_GrandBetaExt_POST,sub1_GrandBetaRest_POST)

% FLX 
% PRE
Alpha_SNR_PRE_FLX = snr(sub1_GrandAlphaFlx_PRE,sub1_GrandAlphaRest_PRE)
Beta_SNR_PRE_FLX = snr(sub1_GrandBetaFlx_PRE,sub1_GrandBetaRest_PRE)
% POST
Alpha_SNR_POST_FLX = snr(sub1_GrandAlphaFlx_POST,sub1_GrandAlphaRest_POST)
Beta_SNR_POST_FLX = snr(sub1_GrandBetaFlx_POST,sub1_GrandBetaRest_POST)


figure;
X = categorical({'PRE TESS','POST TESS'});
X = reordercats(X,{'PRE TESS', 'POST TESS'});
y = [Alpha_SNR_PRE_EXT Beta_SNR_PRE_EXT; Alpha_SNR_POST_EXT Beta_SNR_POST_EXT];
h = bar(X,y)
set(h, {'DisplayName'}, {'Alpha','Beta'}')
legend()
title("SNR EXT");

figure;
X = categorical({'PRE TESS','POST TESS'});
X = reordercats(X,{'PRE TESS', 'POST TESS'});
y = [Alpha_SNR_PRE_FLX Beta_SNR_PRE_FLX; Alpha_SNR_POST_FLX Beta_SNR_POST_FLX];
h = bar(X,y)
set(h, {'DisplayName'}, {'Alpha','Beta'}')
legend()
title("SNR FLX");



