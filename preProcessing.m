clear
close all
clc

Ext = 100;
Flx = 300;
Rest = 400;

%% build inital Datasets
load p3_subjectData.mat;

% constants
fs = subjectData(1).pre(1).hdr.fs;
channels = subjectData(1).pre(1).hdr.Label;

%% SUBJECT 1 PRE TESS
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

%% SUBJECT 2 PRE TESS
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

%% SUBJECT 1 POST TESS
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

%% SUBJECT 2 POST TESS
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
sub1PRE_EEG = sub1PRE_EEG - mean(sub1PRE_EEG);
sub2PRE_EEG = sub2PRE_EEG - mean(sub2PRE_EEG);
sub1POST_EEG = sub1POST_EEG - mean(sub1POST_EEG); 
sub2POST_EEG = sub2POST_EEG - mean(sub2POST_EEG);

%% SUBJECT 1 PRE:

% Signal already has offset removed
signal = sub1PRE_EEG;
figure;
for i = 1:32
    subplot(8,4,i);
    plot((1:length(signal))./fs,(signal(:,i)))
    title('Ch',i)
    hold on;  
end
sgtitle('Remove Offset') 
% Remove outlier
sub1PRE_EEG = filloutliers((sub1PRE_EEG),'nearest','mean');
signal = sub1PRE_EEG;
figure;
for i = 1:32
    subplot(8,4,i);
    plot((1:length(signal))./fs,(signal(:,i)))
    title('Ch',i)
    hold on;  
end
sgtitle('Remove Outliers') 
% Filter: ~[8-30] Hz
BPF = getBPFilter;
sub1PRE_EEG = BPF(sub1PRE_EEG);
signal = sub1PRE_EEG;
figure;
for i = 1:32
    subplot(8,4,i);
    plot((1:length(signal))./fs,(signal(:,i)))
    title('Ch',i)
    hold on;  
end
sgtitle('Bandpass Filter') 
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

%% SUBJECT 1 POST:

% Signal already has offset removed

% Remove outlier
sub1POST_EEG = filloutliers((sub1POST_EEG),'nearest','mean');

% Filter: ~[8-30] Hz
BPF = getBPFilter;
sub1POST_EEG = BPF(sub1POST_EEG);
 
% Rectify
%sub1POST_EEG = abs(sub1POST_EEG);

%% SUBJECT 2 PRE:

% Signal already has offset removed

% Remove outlier
sub2PRE_EEG = filloutliers((sub2PRE_EEG),'nearest','mean');

% Filter: ~[8-30] Hz
BPF = getBPFilter;
sub2PRE_EEG = BPF(sub2PRE_EEG);
 
% Rectify
%sub2PRE_EEG = abs(sub2PRE_EEG);

%% SUBJECT 2 POST:

% Signal already has offset removed

% Remove outlier
sub2POST_EEG = filloutliers((sub2POST_EEG),'nearest','mean');

% Filter: ~[8-30] Hz
BPF = getBPFilter;
sub2POST_EEG = BPF(sub2POST_EEG);
 
% Rectify
%sub2POST_EEG = abs(sub2POST_EEG);


