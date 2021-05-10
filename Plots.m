%% PLOTS 
clc
close all
%%
channels(13) = {'T9'};
channels(19) = {'T10'};
figure
sgtitle('AVG Power Sub 1');
for i = 1:5
    subplot(2,5,i)
    plot_topography(channels, Pavg_1PRE(:,i+25));
end

for i = 1:5
    subplot(2,5,i+5)
    plot_topography(channels, Pavg_1POST(:,i+25));
end

%%
figure
hold on
plot(sub1PRE_DATA(:,1,25))
plot(GAVG_1_PRE(:,1,1))
plot(gavg(:,1,1))
hold off
legend

%%
task = 1;
% 
% channels(13) = {'T9'};
% channels(19) = {'T10'};
figure
sgtitle('G AVG Power Sub 1 PRE TESS EXT');

% plot_topography.make_contour
% for i = 1:2
max = max(max(Pavg_1PRE(:,:, task)));
min =  min(min(Pavg_1PRE(:,:, task)));
% norm = max(max(Pavg_1PRE(:,:, 2))) - min(min(Pavg_1PRE(:,:, 2)));

for i = 1:size(Pavg_1PRE,1)
%     Pavg_1PRE(i,13, task) = max;
%     Pavg_1PRE(i,19, task) = min;
    subplot(4,2/wsize,i)
    hold on
    str = ['TIME: ' , num2str((i-1)*wsize) , ' :: ' , num2str((i)*wsize)];
    title(str)
    plot_topography(channels, Pavg_1PRE(i,:, task), true);
    hold off
end

%% Power grand Avg topo EXT
% plotz(Pavg_1PRE, wsize, 1, 'Avg Power sub 1 PRE' , channels)
% plotz(Pavg_1POST, wsize, 1, 'Avg Power sub 1 POST' , channels)
plotz(Pavg_2PRE, wsize, 1, 'Avg Power sub 2 PRE' , channels)
plotz(Pavg_2POST, wsize, 1, 'Avg Power sub 2 POST' , channels)

%% Power grand Avg topo FLX
% plotz(Pavg_1PRE, wsize, 2, 'Avg Power sub 1 PRE' , channels)
% plotz(Pavg_1POST, wsize, 2, 'Avg Power sub 1 POST' , channels)
% plotz(Pavg_2PRE(:,MIchannels,:), wsize, 2, 'Avg Power sub 2 PRE' , channels(MIchannels))
plotz(Pavg_2POST(:,MIchannels,:), wsize, 2, 'Avg Power sub 2 POST' , channels(MIchannels))
%% Power grand Avg topo Rest
% plotz(Pavg_1PRE, wsize, 3, 'Avg Power sub 1 PRE' , channels)
% plotz(Pavg_1POST, wsize, 3, 'Avg Power sub 1 POST' , channels)
plotz(Pavg_2PRE, wsize, 3, 'Avg Power sub 2 PRE' , channels)
plotz(Pavg_2POST, wsize, 3, 'Avg Power sub 2 POST' , channels)
%% Grand Varaiance Topo

plotz(Gvar_1PRE, wsize, 1, 'Gvar Sub 1 PRE' , channels)
plotz(Gvar_1POST, wsize, 1, 'Gvar Sub 1 POST' , channels)

%% Alpha & Beta Topo Plots
plotz(beta_Pavg_1PRE, wsize, 1, 'Avg Beta Power Sub 11 Pre' , channels)
plotz(beta_Pavg_1POST, wsize, 1, 'Avg Beta Power Sub 11 Post' , channels)

plotz(alpha_Pavg_1PRE, wsize, 1, 'Avg Alpha Power Sub 11 Pre' , channels)
plotz(alpha_Pavg_1POST, wsize, 1, 'Avg Alpha Power Sub 11 Post' , channels)

%% ERS Topo
plotz(ERS_a_1PRE, wsize, 1, 'ERS Alpha Sub 1 Pre' , channels)
plotz(ERS_a_1PRE, wsize, 1, 'ERS Alpha Sub 1 Post' , channels)

plotz(ERS_b_1PRE, wsize, 1, 'ERS Beta Sub 1 Pre' , channels)
plotz(ERS_b_1POST, wsize, 1, 'ERS Beta Sub 1 Post' , channels)

% plotz(ERS_a_1PRE, wsize, 2, 'ERS Alpha Sub 1 Pre' , channels)
% plotz(ERS_a_1PRE, wsize, 2, 'ERS Alpha Sub 1 Post' , channels)
% 
% plotz(ERS_b_1PRE, wsize, 2, 'ERS Beta Sub 1 Pre' , channels)
% plotz(ERS_b_1POST, wsize, 2, 'ERS Beta Sub 1 Post' , channels)

%% PSD Plots
figure('units','normalized','Position',[0.2,0.65,0.3,0.3])
h = spectrum.welch; % creates the Welch spectrum estimator
% 2 sec segment length -> 1024 samples
% more frequency resolution at higher sample rate 
h.SegmentLength = 1024;
SOIf3=psd(h,sub1_GrandAlphaExt_PRE(:),'Fs',fs); % calculates and plot the one sided PSD
plot(SOIf3); % Plot the one-sided PSD. 
temp =get(gca);
temp.Children(1).Color = 'r';
hold on;

SOIf3=psd(h,sub1_GrandAlphaRest_PRE(:),'Fs',fs); % calculates and plot the one sided PSD
plot(SOIf3); % Plot the one-sided PSD. 
temp =get(gca);
temp.Children(1).Color = 'r';
temp.Children(1).LineStyle = '--';
hold on;

SOIf3=psd(h,sub1_GrandAlphaExt_POST(:),'Fs',fs); % calculates and plot the one sided PSD
plot(SOIf3); % Plot the one-sided PSD. 
temp =get(gca);
temp.Children(1).Color = 'g';
hold on;

SOIf3=psd(h,sub1_GrandAlphaRest_POST(:),'Fs',fs); % calculates and plot the one sided PSD
plot(SOIf3); % Plot the one-sided PSD. 
temp =get(gca);
temp.Children(1).Color = 'g';
temp.Children(1).LineStyle = '--';
hold on;
legend('ext pre', 'rest pre', 'ext post', 'rest post');
% 


figure('units','normalized','Position',[0.2,0.65,0.3,0.3])
% 2 sec segment length -> 1024 samples
% more frequency resolution at higher sample rate 
h.SegmentLength = 1024;
SOIf3=psd(h,sub1_GrandBetaExt_PRE(:),'Fs',fs); % calculates and plot the one sided PSD
plot(SOIf3); % Plot the one-sided PSD. 
temp =get(gca);
temp.Children(1).Color = 'r';
hold on;

SOIf3=psd(h,sub1_GrandBetaRest_PRE(:),'Fs',fs); % calculates and plot the one sided PSD
plot(SOIf3); % Plot the one-sided PSD. 
temp =get(gca);
temp.Children(1).Color = 'r';
temp.Children(1).LineStyle = '--';
hold on;

SOIf3=psd(h,sub1_GrandBetaExt_POST(:),'Fs',fs); % calculates and plot the one sided PSD
plot(SOIf3); % Plot the one-sided PSD. 
temp =get(gca);
temp.Children(1).Color = 'g';
hold on;

SOIf3=psd(h,sub1_GrandBetaRest_POST(:),'Fs',fs); % calculates and plot the one sided PSD
plot(SOIf3); % Plot the one-sided PSD. 
temp =get(gca);
temp.Children(1).Color = 'g';
temp.Children(1).LineStyle = '--';
hold on;
legend('ext pre', 'rest pre', 'ext post', 'rest post');
% SOIf3=psd(h,sub1PRE_DATA(:, c, 75),'Fs',fs); % calculates and plot the one sided PSD
% plot(SOIf3); % Plot the one-sided PSD. 
% temp =get(gca);
% temp.Children(1).Color = 'b';
% legend('flex','ext','rest');
% title("Trial 1 PSD Flex vs Rest Single Channel Sub1 PRE")'






%%
labelA = [{'Ext PRE', 'Ext POST', 'Flex PRE', 'Flex POST'}];
labelB = [{'beta pre Ext', 'beta post Ext', 'beta pre Flex', 'post Flex'}];
time = [2, 4] ;
lab = [30, 50];

figure
subplot(1,2,1)
hold on
title('Alpha ERS')
plot((1:length(Avg_ERS_a_1PRE(:,1)))*wsize, Avg_ERS_a_1PRE(:,1))
plot((1:length(Avg_ERS_a_1PRE(:,1)))*wsize,Avg_ERS_a_1POST(:,1))
plot((1:length(Avg_ERS_a_1PRE(:,1)))*wsize,Avg_ERS_a_1PRE(:,2))
plot((1:length(Avg_ERS_a_1PRE(:,1)))*wsize,Avg_ERS_a_1POST(:,2))
stem(time, lab);
legend(labelA); 
hold off


subplot(1,2,2)
hold on
title('Beta ERS')
plot((1:length(Avg_ERS_a_1PRE(:,1)))*wsize, Avg_ERS_b_1PRE(:,1))
plot((1:length(Avg_ERS_a_1PRE(:,1)))*wsize, Avg_ERS_b_1POST(:,1))
plot((1:length(Avg_ERS_a_1PRE(:,1)))*wsize, Avg_ERS_b_1PRE(:,2))
plot((1:length(Avg_ERS_a_1PRE(:,1)))*wsize, Avg_ERS_b_1POST(:,2))
legend(labelA);
hold off
%%

function plotz(Pavg, Wsize, task, name, channels)
    taskType  = [{'EXT', 'FLX', 'REST'}];
%     timePeriod = [{'Fixation' , 'Task Cue', 'Tast Exct'}];
% 
%         channels(13) = {'T9'};
%         channels(19) = {'T10'};
        figure
        sgtitle([name,  taskType(task)]);
%         labCount =1;
        % plot_topography.make_contour
        % for i = 1:2
%         max = max(max(Pavg_1PRE(:,:, task)));
%         min =  min(min(Pavg_1PRE(:,:, task)));
        % norm = max(max(Pavg_1PRE(:,:, 2))) - min(min(Pavg_1PRE(:,:, 2)));

        for i = 1:size(Pavg,1)
        %     Pavg_1PRE(i,13, task) = max;
        %     Pavg_1PRE(i,19, task) = min;
            subplot(3,2/Wsize,i)
            hold on
%             if (mod(i-1, 2/Wsize) == 0)
%                 ylabel(timePeriod(labCount))
%                 labCount = labCount +1;
%             end
            str = ['TIME: ' , num2str((i-1)*Wsize) , ' :: ' , num2str((i)*Wsize)];
            title(str)
%             ylabel('Fix');
            plot_topography(channels, Pavg(i,:, task), true);
            hold off
        end
end
