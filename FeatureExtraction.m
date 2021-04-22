%% Feature Extraction
totalTrials = 75;
totalLength = length(signal);

%% SUB1 PRE
% preallocate feature outputs for speed
[sub1_alphaPower_PRE] = zeros(1,32,totalTrials);
[sub1_alphaMax_PRE] = zeros(1,32,totalTrials);
[sub1_alphaMin_PRE] = zeros(1,32,totalTrials);

[sub1_betaPower_PRE] = zeros(1,32,totalTrials);
[sub1_betaMax_PRE] = zeros(1,32,totalTrials);
[sub1_betaMin_PRE] = zeros(1,32,totalTrials);

signal = sub1PRE_DATA;
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
% preallocate feature outputs for speed
[sub1_alphaPower_POST] = zeros(1,32,totalTrials);
[sub1_alphaMax_POST] = zeros(1,32,totalTrials);
[sub1_alphaMin_POST] = zeros(1,32,totalTrials);

[sub1_betaPower_POST] = zeros(1,32,totalTrials);
[sub1_betaMax_POST] = zeros(1,32,totalTrials);
[sub1_betaMin_POST] = zeros(1,32,totalTrials);

signal = sub1POST_DATA;
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



%% Visualize Features
figure;
for trial = 1:25
    scatter(sub1_alphaMax_PRE(:,:,trial),sub1_betaMax_PRE(:,:,trial), "blue")
    hold on;
    scatter(sub1_alphaMax_POST(:,:,trial),sub1_betaMax_POST(:,:,trial), "red")
end
 

xlabel('mean Alpha power') 
ylabel('mean Beta power') 


