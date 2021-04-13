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
channels = subjectData(1).pre(1).hdr.Label

%% build data for sub 1 pre TESS

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

%% build data for sub 2 pre TESS

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

%% build data for sub 1 post TESS
post1 = subjectData(1).post(1).eeg- mean(subjectData(1).post(1).eeg);
post2 = subjectData(1).post(2).eeg- mean(subjectData(1).post(2).eeg);
posOFFSET = size(post1, 1);

typ1 = subjectData(1).post(1).hdr.EVENT.TYP;
typ2 = subjectData(1).post(2).hdr.EVENT.TYP;

pos1 = subjectData(1).post(1).hdr.EVENT.POS;
pos2 = subjectData(1).post(2).hdr.EVENT.POS;
pos2 = pos2+posOFFSET;

sub1POST_EEG = [post1; post2];
sub1POST_TYP = [typ1; typ2];
sub1POST_POS = [pos1; pos2];

%% build data for sub 2 post TESS
post1 = subjectData(2).post(1).eeg- mean(subjectData(2).post(1).eeg);
post2 = subjectData(2).post(2).eeg- mean(subjectData(2).post(2).eeg);
posOFFSET = size(post1, 1);

typ1 = subjectData(2).post(1).hdr.EVENT.TYP;
typ2 = subjectData(2).post(2).hdr.EVENT.TYP;

pos1 = subjectData(2).post(1).hdr.EVENT.POS;
pos2 = subjectData(2).post(2).hdr.EVENT.POS;
pos2 = pos2+posOFFSET;

sub2POST_EEG = [post1; post2];
sub2POST_TYP = [typ1; typ2];
sub2POST_POS = [pos1; pos2];

%% remove mean
%sub1PRE_EEG = sub1PRE_EEG - mean(sub1PRE_EEG);
%sub2PRE_EEG = sub2PRE_EEG - mean(sub2PRE_EEG);
%sub1POST_EEG = sub1POST_EEG - mean(sub1POST_EEG); 
%sub2POST_EEG = sub2POST_EEG - mean(sub2POST_EEG);

%% 3 task all channels

% ext
sub1PRE_TYP(20:21) % purpose???

% Plot 1 trial of each task
figure('units','normalized','Position',[0.1,0.5,0.3,0.3])

% Subject 1 flex 
s = 12;
e = 16;
subplot(3,1,1)
hold on
plot((sub1PRE_POS(s):sub1PRE_POS(e))./fs, sub1PRE_EEG(sub1PRE_POS(s):sub1PRE_POS(e), :));
stem(sub1PRE_POS(s:e)/fs, sub1PRE_TYP(s:e), 'filled');
xlabel('Seconds');
legend(channels);
title('Subject-1 FLX');
hold off

% Subject 1 rest
s = 7;
e = 11;
subplot(3,1,2)
hold on
plot((sub1PRE_POS(s):sub1PRE_POS(e))./fs, sub1PRE_EEG(sub1PRE_POS(s):sub1PRE_POS(e), :));
stem(sub1PRE_POS(s:e)/fs, sub1PRE_TYP(s:e), 'filled');
xlabel('Seconds');
title('Subject-1 REST');
hold off

% Subject 1 ext
s = 17;
e = 21;
subplot(3,1,3)
hold on
plot((sub1PRE_POS(s):sub1PRE_POS(e))./fs, sub1PRE_EEG(sub1PRE_POS(s):sub1PRE_POS(e), :));
stem(sub1PRE_POS(s:e)/fs, sub1PRE_TYP(s:e), 'filled');
xlabel('Seconds');
title('Subject-1 EXT');
hold off

%% PSD of all channels for each type of movement
h = spectrum.welch; % creates the Welch spectrum estimator
flx = psd(h,sub1PRE_EEG(sub1PRE_POS(5):sub1PRE_POS(6), :),'Fs',fs); 
rest = psd(h,sub1PRE_EEG(sub1PRE_POS(10):sub1PRE_POS(11), :),'Fs',fs);
ext = psd(h,sub1PRE_EEG(sub1PRE_POS(20):sub1PRE_POS(21), :),'Fs',fs);
figure('units','normalized','Position',[0.1,0.5,0.3,0.3])

subplot(3,1,1)
title('FLX');
plot(flx);
subplot(3,1,2)
title('EXT');
plot(ext);
subplot(3,1,3)
title('REST');
plot(rest);
% plot(rest);
% plot(ext);
% temp =get(gca);
% temp.Children(1).Color = 'r'; %rest
% temp.Children(2).Color = 'b'; %flx_dist
hold off
% legend( 'Distal Flexor', 'Rest');


%% PSD of Flex vs Rest for all channels
figure('units','normalized','Position',[0.1,0.5,0.3,0.3])
for i = 1:32
    flx = psd(h,sub1PRE_EEG(sub1PRE_POS(5):sub1PRE_POS(6), i),'Fs',fs); 
    rest = psd(h,sub1PRE_EEG(sub1PRE_POS(10):sub1PRE_POS(11), i),'Fs',fs);
    subplot(8,4,i)
    hold on
    plot(flx)
    plot(rest)
    temp =get(gca);
    temp.Children(1).Color = 'r'; 
    temp.Children(2).Color = 'b'; 
    hold off
end
title("PSD of Flex vs Rest for all channels")

%% PSD of Flex vs Rest for single channel
figure('units','normalized','Position',[0.1,0.5,0.3,0.3])
flx = psd(h,sub1PRE_EEG(sub1PRE_POS(5):sub1PRE_POS(6), 6),'Fs',fs); 
rest = psd(h,sub1PRE_EEG(sub1PRE_POS(10):sub1PRE_POS(11), 6),'Fs',fs);
hold on
plot(flx)
plot(rest)
temp =get(gca);
temp.Children(1).Color = 'r'; 
temp.Children(2).Color = 'b'; 
hold off
title("PSD of Flex vs Rest for channel 6")

