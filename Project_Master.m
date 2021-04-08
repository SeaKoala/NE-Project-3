clear
close all
clc

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

%%


