%% Feature Extraction

totalTrials = 75;

%% SUB1 PRE
% preallocate outputs for speed
[sub1_meanAlpha_PRE] = zeros(1,32,totalTrials);
[sub1_stdAlpha_PRE] = zeros(1,32,totalTrials);
[sub1_meanBeta_PRE] = zeros(1,32,totalTrials);
[sub1_stdBeta_PRE] = zeros(1,32,totalTrials);

 % DWT, 5 level wavelet db2
 signal = sub1PRE_DATA;
for trial = 1:totalTrials
    for channel = 1:32
        [meanAlpha,stdAlpha,meanBeta,stdBeta] = DWT(signal(:,channel,trial));
        sub1_meanAlpha_PRE(1,channel,trial) = meanAlpha;
        sub1_stdAlpha_PRE(1,channel,trial) = stdAlpha;
        sub1_meanBeta_PRE(1,channel,trial) = meanBeta;
        sub1_stdBeta_PRE(1,channel,trial) = stdBeta;
    end
end

%% SUB1 POST
% preallocate outputs for speed
[sub1_meanAlpha_POST] = zeros(1,32,totalTrials);
[sub1_stdAlpha_POST] = zeros(1,32,totalTrials);
[sub1_meanBeta_POST] = zeros(1,32,totalTrials);
[sub1_stdBeta_POST] = zeros(1,32,totalTrials);

 % DWT, 5 level wavelet db2
signal = sub1POST_DATA;
for trial = 1:totalTrials
    for channel = 1:32
        [meanAlpha,stdAlpha,meanBeta,stdBeta] = DWT(signal(:,channel,trial));
        sub1_meanAlpha_POST(1,channel,trial) = meanAlpha;
        sub1_stdAlpha_POST(1,channel,trial) = stdAlpha;
        sub1_meanBeta_POST(1,channel,trial) = meanBeta;
        sub1_stdBeta_POST(1,channel,trial) = stdBeta;
    end
end

%% SUB2 PRE
% preallocate outputs for speed
[sub2_meanAlpha_PRE] = zeros(1,32,totalTrials);
[sub2_stdAlpha_PRE] = zeros(1,32,totalTrials);
[sub2_meanBeta_PRE] = zeros(1,32,totalTrials);
[sub2_stdBeta_PRE] = zeros(1,32,totalTrials);

 % DWT, 5 level wavelet db2
 signal = sub2PRE_DATA;
for trial = 1:totalTrials
    for channel = 1:32
        [meanAlpha,stdAlpha,meanBeta,stdBeta] = DWT(signal(:,channel,trial));
        sub2_meanAlpha_PRE(1,channel,trial) = meanAlpha;
        sub2_stdAlpha_PRE(1,channel,trial) = stdAlpha;
        sub2_meanBeta_PRE(1,channel,trial) = meanBeta;
        sub2_stdBeta_PRE(1,channel,trial) = stdBeta;
    end
end

%% SUB2 POST
% preallocate outputs for speed
[sub2_meanAlpha_POST] = zeros(1,32,totalTrials);
[sub2_stdAlpha_POST] = zeros(1,32,totalTrials);
[sub2_meanBeta_POST] = zeros(1,32,totalTrials);
[sub2_stdBeta_POST] = zeros(1,32,totalTrials);

 % DWT, 5 level wavelet db2
signal = sub2POST_DATA;
for trial = 1:totalTrials
    for channel = 1:32
        [meanAlpha,stdAlpha,meanBeta,stdBeta] = DWT(signal(:,channel,trial));
        sub2_meanAlpha_POST(1,channel,trial) = meanAlpha;
        sub2_stdAlpha_POST(1,channel,trial) = stdAlpha;
        sub2_meanBeta_POST(1,channel,trial) = meanBeta;
        sub2_stdBeta_POST(1,channel,trial) = stdBeta;
    end
end



%% Visualize Features
figure;
for trial = 1:25
    scatter(sub1_meanBeta_PRE(:,:,trial),sub1_stdBeta_PRE(:,:,trial), "blue")
    hold on;
    scatter(sub1_meanBeta_POST(:,:,trial),sub1_stdBeta_POST(:,:,trial), "red")
end




