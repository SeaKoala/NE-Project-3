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

%%
sub1PRE_DATA = slice_n_dice(sub1PRE_EEG, sub1PRE_TYP, sub1PRE_POS); 
sub1POST_DATA = slice_n_dice(sub1POST_EEG, sub1POST_TYP, sub1POST_POS); 
sub2PRE_DATA = slice_n_dice(sub2PRE_EEG, sub2PRE_TYP, sub2PRE_POS); 
sub2POST_DATA = slice_n_dice(sub2POST_EEG, sub2POST_TYP, sub2POST_POS); 

taskType  = [{'EXT', 'FLX', 'REST'}];
results(1:25) = taskType(1);
results(26:50) = taskType(2);
results(51:75) = taskType(3);
% save sub1PRE_DATA sub1PRE_DATA;
% save sub1POST_DATA sub1POST_DATA;
% save sub2PRE_DATA sub2PRE_DATA;
save sub2POST_DATA sub2POST_DATA;

%% Grand Average
GAVG_1_PRE = Gavg(sub1PRE_DATA);
GAVG_1_POST = Gavg(sub1POST_DATA);
GAVG_2_PRE = Gavg(sub2PRE_DATA);
GAVG_2_POST = Gavg(sub2POST_DATA);
%%
wsize = 1;

[Pavg_1PRE] = segmentation(GAVG_1_PRE, wsize, 0);
[Pavg_1POST] = segmentation(GAVG_1_POST, wsize, 0);
[Pavg_2PRE] = segmentation(GAVG_2_PRE, wsize, 0);
[Pavg_2POST] = segmentation(GAVG_2_POST, wsize, 0);

%% Grand Variance
Gvar_1PRE = Gvar(sub1PRE_DATA);
Gvar_1POST = Gvar(sub1POST_DATA);
Gvar_2PRE = Gvar(sub2PRE_DATA);
Gvar_2POST = Gvar(sub2POST_DATA);


% %% ERS
% wsize = 1;
% GAVGAlphaPRE = Gavg(sub1_alphaSig_PRE);
% GAVGAlphaPOST = Gavg(sub1_alphaSig_POST);
% GAVGBetaPRE = Gavg(sub1_betaSig_PRE);
% GAVGBetaPOST = Gavg(sub1_betaSig_POST);
% 
% [alpha_Pavg_1PRE] = segmentation(GAVGAlphaPRE, wsize, 0);
% [alpha_Pavg_1POST] = segmentation(GAVGAlphaPOST, wsize, 0);
% [beta_Pavg_1PRE] = segmentation(GAVGBetaPRE, wsize, 0);
% [beta_Pavg_1POST] = segmentation(GAVGBetaPOST, wsize, 0);
% 
% ERS_a_1PRE = zeros(size(alpha_Pavg_1PRE, 1), 32, 2);
% ERS_a_1POST = zeros(size(alpha_Pavg_1PRE, 1), 32, 2);
% ERS_b_1PRE = zeros(size(alpha_Pavg_1PRE, 1), 32, 2);
% ERS_b_1POST = zeros(size(alpha_Pavg_1PRE, 1), 32, 2);
% for c = 1:32
%     for i = 1:size(alpha_Pavg_1PRE, 1)
%         ERS_a_1PRE(i, c, 1) = 100 * abs(alpha_Pavg_1PRE(i, c, 3) - alpha_Pavg_1PRE(i, c, 1))/abs(alpha_Pavg_1PRE(i, c, 3));
%         ERS_a_1POST(i, c, 1) = 100 * abs(alpha_Pavg_1POST(i, c, 3) - alpha_Pavg_1POST(i, c, 1))/abs(alpha_Pavg_1POST(i, c, 3));
%         ERS_b_1PRE(i, c, 1) = 100 * abs(beta_Pavg_1PRE(i, c, 3) - beta_Pavg_1PRE(i, c, 1))/abs(beta_Pavg_1PRE(i, c, 3));
%         ERS_b_1POST(i, c, 1) = 100 * abs(beta_Pavg_1POST(i, c, 3) - beta_Pavg_1POST(i, c, 1))/abs(beta_Pavg_1POST(i, c, 3));
%         ERS_a_1PRE(i, c, 2) = 100 * abs(alpha_Pavg_1PRE(i, c, 3) - alpha_Pavg_1PRE(i, c, 2))/abs(alpha_Pavg_1PRE(i, c, 3));
%         ERS_a_1POST(i, c, 2) = 100 * abs(alpha_Pavg_1POST(i, c, 3) - alpha_Pavg_1POST(i, c, 2))/abs(alpha_Pavg_1POST(i, c, 3));
%         ERS_b_1PRE(i, c, 2) = 100 * abs(beta_Pavg_1PRE(i, c, 3) - beta_Pavg_1PRE(i, c, 2))/abs(beta_Pavg_1PRE(i, c, 3));
%         ERS_b_1POST(i, c, 2) = 100 * abs(beta_Pavg_1POST(i, c, 3) - beta_Pavg_1POST(i, c, 2))/abs(beta_Pavg_1POST(i, c, 3));
%     end
% end
% for i = 1:size(alpha_Pavg_1PRE, 1)
%     for t= 1:2
%     Avg_ERS_a_1PRE(i,t) = mean(ERS_a_1PRE(i, :, t));
%     Avg_ERS_a_1POST(i,t) = mean(ERS_a_1POST(i, :, t));
%     Avg_ERS_b_1PRE(i,t) = mean(ERS_b_1PRE(i, :, t));
%     Avg_ERS_b_1POST(i,t) = mean(ERS_b_1POST(i, :, t));
%     end
% end
% %% 
% GAVGAlphaPRE = Gavg(sub1_alphaPower_PRE);
% GAVGAlphaPOST = Gavg(sub1_alphaPower_POST);
% GAVGBetaPRE = Gavg(sub1_betaPower_PRE);
% GAVGBetaPOST = Gavg(sub1_betaPower_POST);
% 
% [alpha_Pavg_1PRE] = segmentation(GAVGAlphaPRE, 1, 0);
% [alpha_Pavg_1POST] = segmentation(GAVGAlphaPOST, 1, 0);
% [beta_Pavg_1PRE] = segmentation(GAVGBetaPRE, 1, 0);
%[beta_Pavg_1POST] = segmentation(GAVGBetaPOST, 1, 0);

% erdAlpha_POST_flx = 100 * abs(rest - action)/rest;
%% 
function [x_DATA] = slice_n_dice(eeg, typ, pos) 
    extraSamp = 0
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
        sizes(i) = size(eeg(pos(startIndexs(i))-extraSamp:pos(stopIndexs(i)), :), 1);
    end
    max_samples = (max(sizes))

    x_DATA = zeros(max_samples, 32,totalTrials);

    for i = 1:75
        x_DATA(1:sizes(i),:,i) = eeg(pos(startIndexs(i))-extraSamp:pos(stopIndexs(i)), :); 
    end
end

function [Pavg] = segmentation(data, Wsize, Olap);
    
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


function gvar = Gvar(x)
    gvar = zeros(size(x, 1), 32, 3);

    for i = 1:size(x, 1)
        for c =1:32
            gvar(i,c, 1) = var(x(i, c, 1:25));
            gvar(i,c, 2) = var(x(i, c, 26:50));
            gvar(i,c, 3) = var(x(i, c, 51:75));
        end
    end
end