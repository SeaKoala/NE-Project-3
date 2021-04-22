clear
close all
clc

Ext = 100;
Flx = 300;
Rest = 400;

% build inital Datasets
load p3_subjectData.mat;

% constants
fs = subjectData(1).pre(1).hdr.fs;
channels = subjectData(1).pre(1).hdr.Label;

% SUBJECT 1 PRE TESS
% remove mean to normalize
pre1 = subjectData(1).pre(1).eeg  - mean(subjectData(1).pre(1).eeg);
pre2 = subjectData(1).pre(2).eeg - mean(subjectData(1).pre(2).eeg);
% get position offset to store position of all samples
posOFFSET = size(pre1, 1);
% get event labels
typ1 = subjectData(1).pre(1).hdr.EVENT.TYP;
typ2 = subjectData(1).pre(2).hdr.EVENT.TYP;
% get position of labels
pos1 = subjectData(1).pre(1).hdr.EVENT.POS;
pos2 = subjectData(1).pre(2).hdr.EVENT.POS;
pos2 = pos2+posOFFSET;
% build dataset for sub1 pre TESS 
sub1PRE_EEG = [pre1; pre2];
sub1PRE_TYP = [typ1; typ2];
sub1PRE_POS = [pos1; pos2];

% SUBJECT 2 PRE TESS
% remove mean to normalize
pre1 = subjectData(2).pre(1).eeg- mean(subjectData(2).pre(1).eeg);
pre2 = subjectData(2).pre(2).eeg- mean(subjectData(2).pre(2).eeg);
% get position offset to store position of all samples
posOFFSET = size(pre1, 1);
% get event labels
typ1 = subjectData(2).pre(1).hdr.EVENT.TYP;
typ2 = subjectData(2).pre(2).hdr.EVENT.TYP;
% get position of labels
pos1 = subjectData(2).pre(1).hdr.EVENT.POS;
pos2 = subjectData(2).pre(2).hdr.EVENT.POS;
pos2 = pos2+posOFFSET;
% build dataset for sub2 pre TESS 
sub2PRE_EEG = [pre1; pre2];
sub2PRE_TYP = [typ1; typ2];
sub2PRE_POS = [pos1; pos2];

% SUBJECT 1 POST TESS
% remove mean to normalize
post1 = subjectData(1).post(1).eeg- mean(subjectData(1).post(1).eeg);
post2 = subjectData(1).post(2).eeg- mean(subjectData(1).post(2).eeg);
% get position offset to store position of all samples
posOFFSET = size(post1, 1);
% get event labels
typ1 = subjectData(1).post(1).hdr.EVENT.TYP;
typ2 = subjectData(1).post(2).hdr.EVENT.TYP;
% get position of labels
pos1 = subjectData(1).post(1).hdr.EVENT.POS;
pos2 = subjectData(1).post(2).hdr.EVENT.POS;
pos2 = pos2+posOFFSET;
% build dataset for sub1 post TESS 
sub1POST_EEG = [post1; post2];
sub1POST_TYP = [typ1; typ2];
sub1POST_POS = [pos1; pos2];

% SUBJECT 2 POST TESS
% remove mean to normalize
post1 = subjectData(2).post(1).eeg- mean(subjectData(2).post(1).eeg);
post2 = subjectData(2).post(2).eeg- mean(subjectData(2).post(2).eeg);
% get position offset to store position of all samples
posOFFSET = size(post1, 1);
% get event labels
typ1 = subjectData(2).post(1).hdr.EVENT.TYP;
typ2 = subjectData(2).post(2).hdr.EVENT.TYP;
% get position of labels
pos1 = subjectData(2).post(1).hdr.EVENT.POS;
pos2 = subjectData(2).post(2).hdr.EVENT.POS;
pos2 = pos2+posOFFSET;
% build dataset for sub1 post TESS 
sub2POST_EEG = [post1; post2];
sub2POST_TYP = [typ1; typ2];
sub2POST_POS = [pos1; pos2];

%% Remove Outliers and Filter

% Remove outlier
sub1PRE_EEG = filloutliers((sub1PRE_EEG),'nearest','mean');

% Filter: ~[8-30] Hz
BPF = getBPFilter;
sub1PRE_EEG = BPF(sub1PRE_EEG);

% SUBJECT 1 POST:

% Signal already has offset removed

% Remove outlier
sub1POST_EEG = filloutliers((sub1POST_EEG),'nearest','mean');

% Filter: ~[8-30] Hz
BPF = getBPFilter;
sub1POST_EEG = BPF(sub1POST_EEG);
 
% Rectify
%sub1POST_EEG = abs(sub1POST_EEG);

% SUBJECT 2 PRE:

% Signal already has offset removed

% Remove outlier
sub2PRE_EEG = filloutliers((sub2PRE_EEG),'nearest','mean');

% Filter: ~[8-30] Hz
BPF = getBPFilter;
sub2PRE_EEG = BPF(sub2PRE_EEG);
 
% Rectify
%sub2PRE_EEG = abs(sub2PRE_EEG);

% SUBJECT 2 POST:

% Signal already has offset removed

% Remove outlier
sub2POST_EEG = filloutliers((sub2POST_EEG),'nearest','mean');

% Filter: ~[8-30] Hz
BPF = getBPFilter;
sub2POST_EEG = BPF(sub2POST_EEG);
 
% Rectify
%sub2POST_EEG = abs(sub2POST_EEG);

%% Build dataset by Task
% 
% totalTrials = 75;
% 
% extStartIndex_PRE = find(sub1PRE_TYP == 100);
% extStopIndex_PRE = find(sub1PRE_TYP == 102);
% flxStartIndex_PRE = find(sub1PRE_TYP == 300);
% flxStopIndex_PRE = find(sub1PRE_TYP == 302);
% restStartIndex_PRE = find(sub1PRE_TYP == 400);
% restStopIndex_PRE = find(sub1PRE_TYP == 402);
% 
% startIndexs = [extStartIndex_PRE; flxStartIndex_PRE; restStartIndex_PRE];
% stopIndexs = [extStopIndex_PRE; flxStopIndex_PRE; restStopIndex_PRE];
% 
% sizes = zeros(25,1);
% 
% for i = 1:75
%     sizes(i) = size(sub1PRE_EEG(sub1PRE_POS(startIndexs(i))-1024:sub1PRE_POS(stopIndexs(i)), :), 1);
% end
% max_samples = (max(sizes))
% 
% sub1PRE_DATA = zeros(max_samples, 32,totalTrials);
% 
% for i = 1:75
%     sub1PRE_DATA(1:sizes(i),:,i) = sub1PRE_EEG(sub1PRE_POS(startIndexs(i))-1024:sub1PRE_POS(stopIndexs(i)), :); 
% end
% 
% 
% extStartIndex_POST = find(sub1POST_TYP == 101);
% extStopIndex_POST = find(sub1POST_TYP == 102);
% flxStartIndex_POST = find(sub1POST_TYP == 301);
% flxStopIndex_POST = find(sub1POST_TYP == 302);
% restStartIndex_POST = find(sub1POST_TYP == 401);
% restStopIndex_POST = find(sub1POST_TYP == 402);
% 
% startIndexs = [extStartIndex_POST; flxStartIndex_POST; restStartIndex_POST];
% stopIndexs = [extStopIndex_POST; flxStopIndex_POST; restStopIndex_POST];
% 
% sizes = zeros(25,1);
% 
% for i = 1:75
%     sizes(i) = size(sub1POST_EEG(sub1POST_POS(startIndexs(i)):sub1POST_POS(stopIndexs(i)), :), 1);
% end
% max_samples = (max(sizes));
% 
% sub1POST_DATA = zeros(max_samples, 32,totalTrials);
% 
% for i = 1:75
%     sub1POST_DATA(1:sizes(i),:,i) = sub1POST_EEG(sub1POST_POS(startIndexs(i)):sub1POST_POS(stopIndexs(i)), :); 
% end
% 
% extStartIndex_PRE = find(sub2PRE_TYP == 101);
% extStopIndex_PRE = find(sub2PRE_TYP == 102);
% flxStartIndex_PRE = find(sub2PRE_TYP == 301);
% flxStopIndex_PRE = find(sub2PRE_TYP == 302);
% restStartIndex_PRE = find(sub2PRE_TYP == 401);
% restStopIndex_PRE = find(sub2PRE_TYP == 402);
% 
% startIndexs = [extStartIndex_PRE; flxStartIndex_PRE; restStartIndex_PRE];
% stopIndexs = [extStopIndex_PRE; flxStopIndex_PRE; restStopIndex_PRE];
% 
% sizes = zeros(25,1);
% 
% for i = 1:75
%     sizes(i) = size(sub2PRE_EEG(sub2PRE_POS(startIndexs(i)):sub2PRE_POS(stopIndexs(i)), :), 1);
% end
% max_samples = (max(sizes));
% 
% sub2PRE_DATA = zeros(max_samples, 32,totalTrials);
% 
% for i = 1:75
%     sub2PRE_DATA(1:sizes(i),:,i) = sub2PRE_EEG(sub2PRE_POS(startIndexs(i)):sub2PRE_POS(stopIndexs(i)), :); 
% end
% 
% 
% extStartIndex_POST = find(sub2POST_TYP == 101);
% extStopIndex_POST = find(sub2POST_TYP == 102);
% flxStartIndex_POST = find(sub2POST_TYP == 301);
% flxStopIndex_POST = find(sub2POST_TYP == 302);
% restStartIndex_POST = find(sub2POST_TYP == 401);
% restStopIndex_POST = find(sub2POST_TYP == 402);
% 
% startIndexs = [extStartIndex_POST; flxStartIndex_POST; restStartIndex_POST];
% stopIndexs = [extStopIndex_POST; flxStopIndex_POST; restStopIndex_POST];
% 
% sizes = zeros(25,1);
% 
% for i = 1:75
%     sizes(i) = size(sub2POST_EEG(sub2POST_POS(startIndexs(i)):sub2POST_POS(stopIndexs(i)), :), 1);
% end
% max_samples = (max(sizes));
% 
% sub2POST_DATA = zeros(max_samples, 32,totalTrials);
% 
% for i = 1:75
%     sub2POST_DATA(1:sizes(i),:,i) = sub2POST_EEG(sub2POST_POS(startIndexs(i)):sub2POST_POS(stopIndexs(i)), :); 
% end



%%
sub1PRE_DATA = slice_n_dice(sub1PRE_EEG, sub1PRE_TYP, sub1PRE_POS); 
sub1POST_DATA = slice_n_dice(sub1POST_EEG, sub1POST_TYP, sub1POST_POS); 
sub2PRE_DATA = slice_n_dice(sub2PRE_EEG, sub2PRE_TYP, sub2PRE_POS); 
sub2POST_DATA = slice_n_dice(sub2POST_EEG, sub2POST_TYP, sub2POST_POS); 

taskType  = [{'EXT', 'FLX', 'REST'}];
results(1:25) = taskType(1);
results(26:50) = taskType(2);
results(51:75) = taskType(3);

%% Grand Average
% GAVG_1_PRE = zeros(size(sub1PRE_DATA, 1), 32, 3);
% GAVG_1_POST = zeros(size(sub1POST_DATA, 1), 32, 3);
% GAVG_2_PRE = zeros(size(sub2PRE_DATA, 1), 32, 3);
% GAVG_2_POST = zeros(size(sub2POST_DATA, 1), 32, 3);
% 
% for i = 1:size(sub1PRE_DATA, 1)
%     for t =1:25
%         GAVG_1_PRE(i,:, 1) = sub1PRE_DATA(i, :, t);
%         GAVG_1_PRE(i,:, 2) = sub1PRE_DATA(i, :, t+25);
%         GAVG_1_PRE(i,:, 3) = sub1PRE_DATA(i, :, t+50);
%     end
% end
% for i = 1:size(sub1POST_DATA, 1)
%     for t =1:25
%         GAVG_1_POST(i,:, 1) = sub1POST_DATA(i, :, t);
%         GAVG_1_POST(i,:, 2) = sub1POST_DATA(i, :, t+25);
%         GAVG_1_POST(i,:, 3) = sub1POST_DATA(i, :, t+50);
%     end
% end
% for i = 1:size(sub2PRE_DATA, 1)
%     for t =1:25
%         GAVG_2_PRE(i,:, 1) = sub2PRE_DATA(i, :, t);
%         GAVG_2_PRE(i,:, 2) = sub2PRE_DATA(i, :, t+25);
%         GAVG_2_PRE(i,:, 3) = sub2PRE_DATA(i, :, t+50);
%     end
% end
% for i = 1:size(sub2POST_DATA, 1)
%     for t =1:25
%         GAVG_2_POST(i,:, 1) = sub2POST_DATA(i, :, t);
%         GAVG_2_POST(i,:, 2) = sub2POST_DATA(i, :, t+25);
%         GAVG_2_POST(i,:, 3) = sub2POST_DATA(i, :, t+50);
%     end
% end
GAVG_1_PRE = Gavg(sub1PRE_DATA);
GAVG_1_POST = Gavg(sub1POST_DATA);
GAVG_2_PRE = Gavg(sub2PRE_DATA);
GAVG_2_POST = Gavg(sub2POST_DATA);
%%
wsize = 1;

[Pavg_1PRE] = segmentation(GAVG_1_PRE, wsize, 0, 1, 'alpha', channels);
[Pavg_1POST] = segmentation(GAVG_1_POST, wsize, 0);
[Pavg_2PRE] = segmentation(GAVG_2_PRE, wsize, 0);
[Pavg_2POST] = segmentation(GAVG_2_POST, wsize, 0);


%% Feature Extraction
[var_1PRE, kurt_1PRE, skew_1PRE, mse_1PRE, mae_1PRE, SNR_1PRE, PSNR_1PRE,  Pavg_1PRE]  = feat_extract(sub1PRE_DATA);
[var_1POST, kurt_1POST, skew_1POST, mse_1POST, mae_1POST, SNR_1POST, PSNR_1POST,  Pavg_1POST]  = feat_extract(sub1POST_DATA);
[var_2PRE, kurt_2PRE, skew_2PRE, mse_2PRE, mae_2PRE, SNR_2PRE, PSNR_2PRE,   Pavg_2PRE]  = feat_extract(sub2PRE_DATA);
[var_2POST, kurt_2POST, skew_2POST, mse_2POST, mae_2POST, SNR_2POST, PSNR_2POST,   Pavg_2POST]  = feat_extract(sub2POST_DATA);

%% 
function [var, kurt, skew, mse, mae, SNR, PSNR, Pavg]  = feat_extract(data) 
    
    samp = 512; % number of samples for power averaging
    x = data;
    size(x);
    var = zeros(32, 75);
    kurt = zeros(32, 75);
    skew = zeros(32, 75);
    mse = zeros(32,50);
    mae = zeros(32,50);
    SNR = zeros(32,50);
    PSNR = zeros(32,50);
    Pavg  = zeros(32, 75);  % avg power for first 2 seconds of task
    
    for c = 1:32
        for j = 1:75
            var(c,j) = jfeeg('var', x(:,c,j));
            kurt(c,j) = jfeeg('kurt', x(:,c,j));
            skew(c,j) = jfeeg('skew', x(:,c,j));
            Pavg(c,j) = 1/samp * sum((x(1:samp,c,j)).^2);


            if j<= 25
            temp = x(:, c, j); % signal
            y = x(:, c, j+50); % rest

            mse(c,j)=0;
            for i=1:length(temp)
            mse(c,j)=mse(c,j)+(y(i)-temp(i))^2;
            end
            mse(c,j)=mse(c,j)/length(temp);
            %MAE %Mean absolute error
            mae(c,j)=0;
            for i=1:length(temp)
            mae(c,j)=mae(c,j)+abs(y(i)-temp(i));
            end
            mae(c,j)=mae(c,j)/length(temp);

            num=0;
            den=0;
            for i=1:length(temp)
            den=den+(y(i)-temp(i))^2;
            end
            for i=1:length (temp)
            num=num+temp(i)^2;
            end
            SNR(c,j) = 20*log10(sqrt(num)/sqrt(den));
            PSNR(c, j)= 20*log10(max(temp)/sqrt(mse(c,j)));
            
            temp = x(:, c, j+25); % signal
            y = x(:, c, j+50); % rest

            mse(c,j+25)=0;
            for i=1:length(temp)
            mse(c,j+25)=mse(c,j+25)+(y(i)-temp(i))^2;
            end
            mse(c,j+25)=mse(c,j+25)/length(temp);
            %MAE %Mean absolute error
            mae(c,j+25)=0;
            for i=1:length(temp)
            mae(c,j+25)=mae(c,j+25)+abs(y(i)-temp(i));
            end
            mae(c,j+25)=mae(c,j+25)/length(temp);

            num=0;
            den=0;
            for i=1:length(temp)
            den=den+(y(i)-temp(i))^2;
            end
            for i=1:length (temp)
            num=num+temp(i)^2;
            end
            SNR(c,j+25) = 20*log10(sqrt(num)/sqrt(den));
            PSNR(c, j+25)= 20*log10(max(temp)/sqrt(mse(c,j+25)));
            end
        end
    end
end


function [x_DATA] = slice_n_dice(eeg, typ, pos) 

    totalTrials = 75;

    extStartIndex_PRE = find(typ == 100);
    extStopIndex_PRE = find(typ == 102);
    flxStartIndex_PRE = find(typ == 300);
    flxStopIndex_PRE = find(typ == 302);
    restStartIndex_PRE = find(typ == 400);
    restStopIndex_PRE = find(typ == 402);

    startIndexs = [extStartIndex_PRE; flxStartIndex_PRE; restStartIndex_PRE];
    stopIndexs = [extStopIndex_PRE; flxStopIndex_PRE; restStopIndex_PRE];

    sizes = zeros(25,1);

    for i = 1:75
        sizes(i) = size(eeg(pos(startIndexs(i))-1024:pos(stopIndexs(i)), :), 1);
    end
    max_samples = (max(sizes))

    x_DATA = zeros(max_samples, 32,totalTrials);

    for i = 1:75
        x_DATA(1:sizes(i),:,i) = eeg(pos(startIndexs(i))-1024:pos(stopIndexs(i)), :); 
    end
end

function [Pavg] = segmentation(data, Wsize, Olap, task, name, channels);
    
    fs = 512;
    WSize = floor(Wsize*fs);	    % length of each data frame, 30ms
    nOlap = floor(Olap*WSize);  % overlap of successive frames, half of WSize
    hop = WSize-nOlap;	    % amount to advance for next data frame
    nx = length(data(:,1,1));            % length of input vector
    len = fix((nx - (WSize-hop))/hop);	%length of output vector = total frames
    
%     [MAV, WT] = deal(zeros(tasks, len, 4));
%     [MAV_max, WT_max] = deal(zeros(tasks, 4));
    Pavg = zeros(len, 32, 3);
        
    for T = 1:3
        for c = 1:32
            for i = 1:len
                segment = data(         (((i-1)*hop+1):((i-1)*hop+WSize)),  c,          T);
                Pavg(i, c, T) = 1/length(segment) * sum(      segment(1:length(segment)).^2          );
    %             MAV(t, i, p) = mean(abs(segment));
    %             WT(t, i, p) = jfemg('ewl', segment);
            end
        end
    end
end

function gavg = Gavg(x)
    gavg = zeros(size(x, 1), 32, 3);

    for i = 1:size(x, 1)
        for t =1:25
            gavg(i,:, 1) = gavg(i,:, 1) + x(i, :, t);
            gavg(i,:, 2) = gavg(i,:, 2) +  x(i, :, t+25);
            gavg(i,:, 3) = gavg(i,:, 3)+ x(i, :, t+50);
        end
    end
        gavg = gavg(:,:, :)./25;
end