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

sub1POST_DATA = zeros(max_samples, 32,totalTrials);

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

sub2POST_DATA = zeros(max_samples, 32,totalTrials);

for i = 1:75
    sub2POST_DATA(1:sizes(i),:,i) = sub2POST_EEG(sub2POST_POS(startIndexs(i)):sub2POST_POS(stopIndexs(i)), :); 
end

taskType  = [{'EXT', 'FLX', 'REST'}];
results(1:25) = taskType(1);
results(26:50) = taskType(2);
results(51:75) = taskType(3);


%% Grand Average
GAVG_1_PRE = zeros(size(sub1PRE_DATA, 1), 32, 3);
GAVG_1_POST = zeros(size(sub1POST_DATA, 1), 32, 3);
GAVG_2_PRE = zeros(size(sub2PRE_DATA, 1), 32, 3);
GAVG_2_POST = zeros(size(sub2POST_DATA, 1), 32, 3);

for i = 1:size(sub1PRE_DATA, 1)
    for t =1:25
        GAVG_1_PRE(i,:, 1) = sub1PRE_DATA(i, :, t);
        GAVG_1_PRE(i,:, 2) = sub1PRE_DATA(i, :, t+25);
        GAVG_1_PRE(i,:, 3) = sub1PRE_DATA(i, :, t+50);
    end
end
for i = 1:size(sub1POST_DATA, 1)
    for t =1:25
        GAVG_1_POST(i,:, 1) = sub1POST_DATA(i, :, t);
        GAVG_1_POST(i,:, 2) = sub1POST_DATA(i, :, t+25);
        GAVG_1_POST(i,:, 3) = sub1POST_DATA(i, :, t+50);
    end
end
for i = 1:size(sub2PRE_DATA, 1)
    for t =1:25
        GAVG_2_PRE(i,:, 1) = sub2PRE_DATA(i, :, t);
        GAVG_2_PRE(i,:, 2) = sub2PRE_DATA(i, :, t+25);
        GAVG_2_PRE(i,:, 3) = sub2PRE_DATA(i, :, t+50);
    end
end
for i = 1:size(sub2POST_DATA, 1)
    for t =1:25
        GAVG_2_POST(i,:, 1) = sub2POST_DATA(i, :, t);
        GAVG_2_POST(i,:, 2) = sub2POST_DATA(i, :, t+25);
        GAVG_2_POST(i,:, 3) = sub2POST_DATA(i, :, t+50);
    end
end