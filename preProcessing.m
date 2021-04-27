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

%% remove offset against the different runs
% sub1PRE_EEG = sub1PRE_EEG - mean(sub1PRE_EEG);
% sub2PRE_EEG = sub2PRE_EEG - mean(sub2PRE_EEG);
% sub1POST_EEG = sub1POST_EEG - mean(sub1POST_EEG); 
% sub2POST_EEG = sub2POST_EEG - mean(sub2POST_EEG);

%% Remove Outliers and Filter

% Signal already has offset removed
signal = sub1PRE_EEG;
%figure;
% for i = 1:32
%     subplot(8,4,i);
%     plot((1:length(signal))./fs,(signal(:,i)))
%     title('Ch',i)
%     hold on;  
% end
%sgtitle('Remove Offset') 

% Remove outlier
sub1PRE_EEG = filloutliers((sub1PRE_EEG),'nearest','mean');
signal = sub1PRE_EEG;
%figure;
% for i = 1:32
%     subplot(8,4,i);
%     plot((1:length(signal))./fs,(signal(:,i)))
%     title('Ch',i)
%     hold on;  
% end
%sgtitle('Remove Outliers') 

% Filter: ~[8-30] Hz
BPF = getBPFilter;
sub1PRE_EEG = BPF(sub1PRE_EEG);
signal = sub1PRE_EEG;
%figure;
% for i = 1:32
%     subplot(8,4,i);
%     plot((1:length(signal))./fs,(signal(:,i)))
%     title('Ch',i)
%     hold on;  
% end
%sgtitle('Bandpass Filter') 
% Rectify
%sub1PRE_EEG = abs(sub1PRE_EEG);
%signal = sub1PRE_EEG;
%figure;
%for i = 1:32
    %subplot(8,4,i);
    %plot((1:length(signal))./fs,(signal(:,i)))
    %title('Ch',i)
    %hold on;  
%end
%sgtitle('Rectified') 

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


%% Extract task periods

% Subject 1 PRE 
%EXT
extStartIndex_PRE = find(sub1PRE_TYP == 101);
extStopIndex_PRE = find(sub1PRE_TYP == 102);
% Get the first trial of ext
extStartSample_PRE = sub1PRE_POS(extStartIndex_PRE(1));
extStopSample_PRE = sub1PRE_POS(extStopIndex_PRE(1));

% FLX
flxStartIndex_PRE = find(sub1PRE_TYP == 301);
flxStopIndex_PRE = find(sub1PRE_TYP == 302);
% Get the first trial of flx
flxStartSample_PRE = sub1PRE_POS(flxStartIndex_PRE(1));
flxStopSample_PRE = sub1PRE_POS(flxStopIndex_PRE(1));
% REST
restStartIndex_PRE = find(sub1PRE_TYP == 401);
restStopIndex_PRE = find(sub1PRE_TYP == 402);
% Get the first trial of rest
restStartSample_PRE = sub1PRE_POS(restStartIndex_PRE(1));
restStopSample_PRE = sub1PRE_POS(restStopIndex_PRE(1));

% Subject 1 POST
%EXT
extStartIndex_POST = find(sub1POST_TYP == 101);
extStopIndex_POST = find(sub1POST_TYP == 102);
% Get the first trial of ext
extStartSample_POST = sub1POST_POS(extStartIndex_POST(1));
extStopSample_POST = sub1POST_POS(extStopIndex_POST(1));
% FLX
flxStartIndex_POST = find(sub1POST_TYP == 301);
flxStopIndex_POST = find(sub1POST_TYP == 302);
% Get the first trial of flx
flxStartSample_POST = sub1POST_POS(flxStartIndex_POST(1));
flxStopSample_POST = sub1POST_POS(flxStopIndex_POST(1));
% REST
restStartIndex_POST = find(sub1POST_TYP == 401);
restStopIndex_POST = find(sub1POST_TYP == 402);
% Get the first trial of rest
restStartSample_POST = sub1POST_POS(restStartIndex_POST(1));
restStopSample_POST = sub1POST_POS(restStopIndex_POST(1));


%% Build dataset by Task

totalTrials = 75;

extStartIndex_PRE = find(sub1PRE_TYP == 101);
extStopIndex_PRE = find(sub1PRE_TYP == 102);
flxStartIndex_PRE = find(sub1PRE_TYP == 301);
flxStopIndex_PRE = find(sub1PRE_TYP == 302);
restStartIndex_PRE = find(sub1PRE_TYP == 401);
restStopIndex_PRE = find(sub1PRE_TYP == 402);

startIndexs = [extStartIndex_PRE; flxStartIndex_PRE; restStartIndex_PRE];
stopIndexs = [extStopIndex_PRE; flxStopIndex_PRE; restStopIndex_PRE];

sizes = zeros(25,1);

for i = 1:75
    sizes(i) = size(sub1PRE_EEG(sub1PRE_POS(startIndexs(i)):sub1PRE_POS(stopIndexs(i)), :), 1);
end
max_samples = (max(sizes));

sub1PRE_DATA = zeros(max_samples, 32,totalTrials);

for i = 1:75
    sub1PRE_DATA(1:sizes(i),:,i) = sub1PRE_EEG(sub1PRE_POS(startIndexs(i)):sub1PRE_POS(stopIndexs(i)), :); 
end


extStartIndex_POST = find(sub1POST_TYP == 101);
extStopIndex_POST = find(sub1POST_TYP == 102);
flxStartIndex_POST = find(sub1POST_TYP == 301);
flxStopIndex_POST = find(sub1POST_TYP == 302);
restStartIndex_POST = find(sub1POST_TYP == 401);
restStopIndex_POST = find(sub1POST_TYP == 402);

startIndexs = [extStartIndex_POST; flxStartIndex_POST; restStartIndex_POST];
stopIndexs = [extStopIndex_POST; flxStopIndex_POST; restStopIndex_POST];

sizes = zeros(25,1);

for i = 1:75
    sizes(i) = size(sub1POST_EEG(sub1POST_POS(startIndexs(i)):sub1POST_POS(stopIndexs(i)), :), 1);
end
max_samples = (max(sizes));

sub1POST_DATA = zeros(max_samples, 32, totalTrials);

for i = 1:75
    sub1POST_DATA(1:sizes(i),:,i) = sub1POST_EEG(sub1POST_POS(startIndexs(i)):sub1POST_POS(stopIndexs(i)), :); 
end

extStartIndex_PRE = find(sub2PRE_TYP == 101);
extStopIndex_PRE = find(sub2PRE_TYP == 102);
flxStartIndex_PRE = find(sub2PRE_TYP == 301);
flxStopIndex_PRE = find(sub2PRE_TYP == 302);
restStartIndex_PRE = find(sub2PRE_TYP == 401);
restStopIndex_PRE = find(sub2PRE_TYP == 402);

startIndexs = [extStartIndex_PRE; flxStartIndex_PRE; restStartIndex_PRE];
stopIndexs = [extStopIndex_PRE; flxStopIndex_PRE; restStopIndex_PRE];

sizes = zeros(25,1);

for i = 1:75
    sizes(i) = size(sub2PRE_EEG(sub2PRE_POS(startIndexs(i)):sub2PRE_POS(stopIndexs(i)), :), 1);
end
max_samples = (max(sizes));

sub2PRE_DATA = zeros(max_samples, 32,totalTrials);

for i = 1:75
    sub2PRE_DATA(1:sizes(i),:,i) = sub2PRE_EEG(sub2PRE_POS(startIndexs(i)):sub2PRE_POS(stopIndexs(i)), :); 
end


extStartIndex_POST = find(sub2POST_TYP == 101);
extStopIndex_POST = find(sub2POST_TYP == 102);
flxStartIndex_POST = find(sub2POST_TYP == 301);
flxStopIndex_POST = find(sub2POST_TYP == 302);
restStartIndex_POST = find(sub2POST_TYP == 401);
restStopIndex_POST = find(sub2POST_TYP == 402);

startIndexs = [extStartIndex_POST; flxStartIndex_POST; restStartIndex_POST];
stopIndexs = [extStopIndex_POST; flxStopIndex_POST; restStopIndex_POST];

sizes = zeros(25,1);

for i = 1:75
    sizes(i) = size(sub2POST_EEG(sub2POST_POS(startIndexs(i)):sub2POST_POS(stopIndexs(i)), :), 1);
end
max_samples = (max(sizes));

sub2POST_DATA = zeros(max_samples, 32, totalTrials);

for i = 1:75
    sub2POST_DATA(1:sizes(i),:,i) = sub2POST_EEG(sub2POST_POS(startIndexs(i)):sub2POST_POS(stopIndexs(i)), :); 
end

taskType  = [{'EXT', 'FLX', 'REST'}];
results(1:25) = taskType(1);
results(26:50) = taskType(2);
results(51:75) = taskType(3);


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


    
    
    
    
    
    
    
    
    
    
    
    
    