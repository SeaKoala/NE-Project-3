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

% build data for sub 1 pre TESS
pre1 = subjectData(1).pre(1).eeg;
pre2 = subjectData(1).pre(2).eeg;
posOFFSET = size(pre1, 1);

typ1 = subjectData(1).pre(1).hdr.EVENT.TYP;
typ2 = subjectData(1).pre(2).hdr.EVENT.TYP;

pos1 = subjectData(1).pre(1).hdr.EVENT.POS;
pos2 = subjectData(1).pre(2).hdr.EVENT.POS;
pos2 = pos2+posOFFSET;

sub1PRE_EEG = [pre1; pre2];
sub1PRE_TYP = [typ1; typ2];
sub1PRE_POS = [pos1; pos2];

% build data for sub 2 pre TESS
pre1 = subjectData(2).pre(1).eeg;
pre2 = subjectData(2).pre(2).eeg;
posOFFSET = size(pre1, 1);

typ1 = subjectData(2).pre(1).hdr.EVENT.TYP;
typ2 = subjectData(2).pre(2).hdr.EVENT.TYP;

pos1 = subjectData(2).pre(1).hdr.EVENT.POS;
pos2 = subjectData(2).pre(2).hdr.EVENT.POS;
pos2 = pos2+posOFFSET;

sub2PRE_EEG = [pre1; pre2];
sub2PRE_TYP = [typ1; typ2];
sub2PRE_POS = [pos1; pos2];

% build data for sub 1 post TESS
post1 = subjectData(1).post(1).eeg;
post2 = subjectData(1).post(2).eeg;
posOFFSET = size(post1, 1);

typ1 = subjectData(1).post(1).hdr.EVENT.TYP;
typ2 = subjectData(1).post(2).hdr.EVENT.TYP;

pos1 = subjectData(1).post(1).hdr.EVENT.POS;
pos2 = subjectData(1).post(2).hdr.EVENT.POS;
pos2 = pos2+posOFFSET;

sub1POST_EEG = [post1; post2];
sub1POST_TYP = [typ1; typ2];
sub1POST_POS = [pos1; pos2];

% build data for sub 2 post TESS
post1 = subjectData(2).post(1).eeg;
post2 = subjectData(2).post(2).eeg;
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
sub1PRE_EEG = sub1PRE_EEG - mean(sub1PRE_EEG);
sub2PRE_EEG = sub2PRE_EEG - mean(sub2PRE_EEG);
sub1POST_EEG = sub1POST_EEG - mean(sub1POST_EEG); 
sub2POST_EEG = sub2POST_EEG - mean(sub2POST_EEG);

%% 3 tasks all channels
sub1PRE_TYP(20:21)


figure('units','normalized','Position',[0.2,0.65,0.3,0.3])

s = 12;
e = 16;
subplot(3,1,1)
hold on
plot((sub1PRE_POS(s):sub1PRE_POS(e))./fs, sub1PRE_EEG(sub1PRE_POS(s):sub1PRE_POS(e), :));
stem(sub1PRE_POS(s:e)/fs, sub1PRE_TYP(s:e), 'filled');
xlabel('Seconds');
% legend(channels);
title('Subject-1 FLX');
hold off

s = 7;
e = 11;
subplot(3,1,2)
hold on
plot((sub1PRE_POS(s):sub1PRE_POS(e))./fs, sub1PRE_EEG(sub1PRE_POS(s):sub1PRE_POS(e), :));
stem(sub1PRE_POS(s:e)/fs, sub1PRE_TYP(s:e), 'filled');
xlabel('Seconds');
% legend(channels);
title('Subject-1 REST');
hold off

s = 17;
e = 21;
subplot(3,1,3)
hold on
plot((sub1PRE_POS(s):sub1PRE_POS(e))./fs, sub1PRE_EEG(sub1PRE_POS(s):sub1PRE_POS(e), :));
stem(sub1PRE_POS(s:e)/fs, sub1PRE_TYP(s:e), 'filled');
xlabel('Seconds');
% legend(channels);
title('Subject-1 EXT');
hold off
%% PSD
h = spectrum.welch; % creates the Welch spectrum estimator
flx = psd(h,sub1PRE_EEG(sub1PRE_POS(5):sub1PRE_POS(6), :),'Fs',fs); 
rest = psd(h,sub1PRE_EEG(sub1PRE_POS(10):sub1PRE_POS(11), :),'Fs',fs);
ext = psd(h,sub1PRE_EEG(sub1PRE_POS(20):sub1PRE_POS(21), :),'Fs',fs);
figure('units','normalized','Position',[0.2,0.65,0.3,0.3])

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

