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


%% ERD/ERP analysis


%% PRE
% flx
% calculates the mean of 32 channels for each trial
extAlphaGrandAVGPower = mean(sub1_alphaPower_PRE(1,:,1:25));

extBetaGrandAVGPower = mean(sub1_betaPower_PRE(1,:,1:25));


% ext
flxAlphaGrandAVGPower = mean(sub1_alphaPower_PRE(1,:,26:50));

flxBetaGrandAVGPower = mean(sub1_betaPower_PRE(1,:,26:50));


% rest
restAlphaGrandAVGPower = mean(sub1_alphaPower_PRE(1,:,51:75));

restBetaGrandAVGPower = mean(sub1_betaPower_PRE(1,:,51:75));


%% ERD - percent change of power between rest and action
%% EXT
% Alpha
rest = mean(restAlphaGrandAVGPower);
action = mean(extAlphaGrandAVGPower);
erdAlpha_PRE_ext = 100 * (action - rest)/abs(rest);
% Beta
rest = mean(restBetaGrandAVGPower);
action = mean(extBetaGrandAVGPower);
erdBeta_PRE_ext = 100 * (action - rest)/abs(rest);

%% FLX 
% Alpha
rest = mean(restAlphaGrandAVGPower);
action = mean(flxAlphaGrandAVGPower);
erdAlpha_PRE_flx = 100 * (action - rest)/abs(rest);
% Beta
rest = mean(restBetaGrandAVGPower);
action = mean(flxBetaGrandAVGPower);
erdBeta_PRE_flx = 100 * (action - rest)/abs(rest);




%% POST

% ext
% calculates the mean of 32 channels for each trial
extAlphaGrandAVGPower = mean(sub1_alphaPower_POST(1,:,1:25));

extBetaGrandAVGPower = mean(sub1_betaPower_POST(1,:,1:25));


% flx
flxAlphaGrandAVGPower = mean(sub1_alphaPower_POST(1,:,26:50));

flxBetaGrandAVGPower = mean(sub1_betaPower_POST(1,:,26:50));


% rest
restAlphaGrandAVGPower = mean(sub1_alphaPower_POST(1,:,51:75));

restBetaGrandAVGPower = mean(sub1_betaPower_POST(1,:,51:75));


%% ERD - percent change of power between rest and action
%% EXT
% Alpha
rest = mean(restAlphaGrandAVGPower);
action = mean(extAlphaGrandAVGPower);
erdAlpha_POST_ext = 100 * (action - rest)/abs(rest);
% Beta
rest = mean(restBetaGrandAVGPower);
action = mean(extBetaGrandAVGPower);
erdBeta_POST_ext = 100 * (action - rest)/abs(rest);

%% FLX 
% Alpha
rest = mean(restAlphaGrandAVGPower);
action = mean(flxAlphaGrandAVGPower);
erdAlpha_POST_flx = 100 * (action - rest)/abs(rest);
% Beta
rest = mean(restBetaGrandAVGPower);
action = mean(flxBetaGrandAVGPower);
erdBeta_POST_flx = 100 * (action - rest)/abs(rest);

%% Visualize Power changes
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
MIchannels = [6,9,10,11,12,15,16,17,21,22,25,27];

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
ylim([-10 200])


