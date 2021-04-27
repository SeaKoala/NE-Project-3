%% Feature Extraction
totalTrials = 75;

%% SUB1 PRE
signal = sub2PRE_DATA;

% preallocate feature outputs for speed
[sub1_alphaPower_PRE] = zeros(1,32,totalTrials);
[sub1_alphaMax_PRE] = zeros(1,32,totalTrials);
[sub1_alphaMin_PRE] = zeros(1,32,totalTrials);
[sub1_alphaSig_PRE] = zeros(length(signal), 32, totalTrials);


[sub1_betaPower_PRE] = zeros(1,32,totalTrials);
[sub1_betaMax_PRE] = zeros(1,32,totalTrials);
[sub1_betaMin_PRE] = zeros(1,32,totalTrials);
[sub1_betaSig_PRE] = zeros(length(signal), 32, totalTrials);



 % DWT, 6 level wavelet db4
for trial = 1:totalTrials
    for channel = 1:32
        [alphaPower,alphaMax,alphaMin,betaPower,betaMax,betaMin] = DWT(signal(:,channel,trial)); 
        sub1_alphaPower_PRE(1,channel,trial) = alphaPower;
        sub1_alphaMax_PRE(1,channel,trial) = alphaMax;
        sub1_alphaMin_PRE(1,channel,trial) = alphaMin;
        sub1_betaPower_PRE(1,channel,trial) = betaPower;
        sub1_betaMax_PRE(1,channel,trial) = betaMax;
        sub1_betaMin_PRE(1,channel,trial) = betaMin;
    end
end

%% SUB1 POST
signal = sub2POST_DATA;

% preallocate feature outputs for speed
[sub1_alphaPower_POST] = zeros(1,32,totalTrials);
[sub1_alphaMax_POST] = zeros(1,32,totalTrials);
[sub1_alphaMin_POST] = zeros(1,32,totalTrials);
[sub1_alphaSig_POST] = zeros(length(signal), 32, totalTrials);

[sub1_betaPower_POST] = zeros(1,32,totalTrials);
[sub1_betaMax_POST] = zeros(1,32,totalTrials);
[sub1_betaMin_POST] = zeros(1,32,totalTrials);
[sub1_betaSig_POST] = zeros(length(signal), 32, totalTrials);


 % DWT, 6 level wavelet db4
for trial = 1:totalTrials
    for channel = 1:32
        [alphaPower,alphaMax,alphaMin,betaPower,betaMax,betaMin] = DWT(signal(:,channel,trial)); 
        sub1_alphaPower_POST(1,channel,trial) = alphaPower;
        sub1_alphaMax_POST(1,channel,trial) = alphaMax;
        sub1_alphaMin_POST(1,channel,trial) = alphaMin;
        sub1_betaPower_POST(1,channel,trial) = betaPower;
        sub1_betaMax_POST(1,channel,trial) = betaMax;
        sub1_betaMin_POST(1,channel,trial) = betaMin;
    end
end


%% ERD/ERP analysis

% Grand Average Power for each movement

%% PRE
% flx
extAlphaGrandAVGPower = mean(sub1_alphaPower_PRE(1,:,1:25));
extAlphaGrandMIN = mean(sub1_alphaMin_PRE(1,:,1:25));
extAlphaGrandMAX = mean(sub1_alphaMax_PRE(1,:,1:25));
extBetaGrandAVGPower = mean(sub1_betaPower_PRE(1,:,1:25));
extBetaGrandMIN = mean(sub1_betaMin_PRE(1,:,1:25));
extBetaGrandMAX = mean(sub1_betaMax_PRE(1,:,1:25));

% ext
flxAlphaGrandAVGPower = mean(sub1_alphaPower_PRE(1,:,26:50));
flxAlphaGrandMIN = mean(sub1_alphaMin_PRE(1,:,26:50));
flxAlphaGrandMAX = mean(sub1_alphaMax_PRE(1,:,26:50));
flxBetaGrandAVGPower = mean(sub1_betaPower_PRE(1,:,26:50));
flxBetaGrandMIN = mean(sub1_betaMin_PRE(1,:,26:50));
flxBetaGrandMAX = mean(sub1_betaMax_PRE(1,:,26:50));

% rest
restAlphaGrandAVGPower = mean(sub1_alphaPower_PRE(1,:,51:75));
restAlphaGrandMIN = mean(sub1_alphaMin_PRE(1,:,51:75));
restAlphaGrandMAX = mean(sub1_alphaMax_PRE(1,:,51:75));
restBetaGrandAVGPower = mean(sub1_betaPower_PRE(1,:,51:75));
restBetaGrandMIN = mean(sub1_betaMin_PRE(1,:,51:75));
restBetaGrandMAX = mean(sub1_betaMax_PRE(1,:,51:75));

%% ERD - percent change of power between rest and action
%% EXT
% Alpha
rest = mean(restAlphaGrandAVGPower);
action = mean(extAlphaGrandAVGPower);
erdAlpha_PRE_ext = 100 * abs(rest - action)/rest;
% Beta
rest = mean(restBetaGrandAVGPower);
action = mean(extBetaGrandAVGPower);
erdBeta_PRE_ext = 100 * abs(rest - action)/rest;

%% FLX 
% Alpha
rest = mean(restAlphaGrandAVGPower);
action = mean(flxAlphaGrandAVGPower);
erdAlpha_PRE_flx = 100 * abs(rest - action)/rest;
% Beta
rest = mean(restBetaGrandAVGPower);
action = mean(flxBetaGrandAVGPower);
erdBeta_PRE_flx = 100 * abs(rest - action)/rest;


%% POST

% ext
% calculates the mean of 32 channels for each trial
extAlphaGrandAVGPower = mean(sub1_alphaPower_POST(1,:,1:25));
extAlphaGrandMIN = mean(sub1_alphaMin_POST(1,:,1:25));
extAlphaGrandMAX = mean(sub1_alphaMax_POST(1,:,1:25));
extBetaGrandAVGPower = mean(sub1_betaPower_POST(1,:,1:25));
extBetaGrandMIN = mean(sub1_betaMin_POST(1,:,1:25));
extBetaGrandMAX = mean(sub1_betaMax_POST(1,:,1:25));

% flx
flxAlphaGrandAVGPower = mean(sub1_alphaPower_POST(1,:,26:50));
flxAlphaGrandMIN = mean(sub1_alphaMin_POST(1,:,26:50));
flxAlphaGrandMAX = mean(sub1_alphaMax_POST(1,:,26:50));
flxBetaGrandAVGPower = mean(sub1_betaPower_POST(1,:,26:50));
flxBetaGrandMIN = mean(sub1_betaMin_POST(1,:,26:50));
flxBetaGrandMAX = mean(sub1_betaMax_POST(1,:,26:50));

% rest
restAlphaGrandAVGPower = mean(sub1_alphaPower_POST(1,:,51:75));
restAlphaGrandMIN = mean(sub1_alphaMin_POST(1,:,51:75));
restAlphaGrandMAX = mean(sub1_alphaMax_POST(1,:,51:75));
restBetaGrandAVGPower = mean(sub1_betaPower_POST(1,:,51:75));
restBetaGrandMIN = mean(sub1_betaMin_POST(1,:,51:75));
restBetaGrandMAX = mean(sub1_betaMax_POST(1,:,51:75));

%% ERD - percent change of power between rest and action
%% EXT
% Alpha
rest = mean(restAlphaGrandAVGPower);
action = mean(extAlphaGrandAVGPower);
erdAlpha_POST_ext = 100 * abs(rest - action)/rest;
% Beta
rest = mean(restBetaGrandAVGPower);
action = mean(extBetaGrandAVGPower);
erdBeta_POST_ext = 100 * abs(rest - action)/rest;

%% FLX 
% Alpha
rest = mean(restAlphaGrandAVGPower);
action = mean(flxAlphaGrandAVGPower);
erdAlpha_POST_flx = 100 * abs(rest - action)/rest;
% Beta
rest = mean(restBetaGrandAVGPower);
action = mean(flxBetaGrandAVGPower);
erdBeta_POST_flx = 100 * abs(rest - action)/rest;

%% Visualize Features
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



figure;

