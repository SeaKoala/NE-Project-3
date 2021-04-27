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

channels(13) = {'T9'};
channels(19) = {'T10'};
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

%% Power grand Avg topo
% plotz(Pavg_1PRE, wsize, 1, '1 PRE' , channels)
% plotz(Pavg_1POST, wsize, 1, '1 POST' , channels)
plotz(Pavg_2PRE, wsize, 2, '2 PRE' , channels)
plotz(Pavg_2POST, wsize, 2, '2 POST' , channels)
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

        channels(13) = {'T9'};
        channels(19) = {'T10'};
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
            subplot(4,2/Wsize,i)
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
