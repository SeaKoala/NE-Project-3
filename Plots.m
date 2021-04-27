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

%%
% plotz(Pavg_1PRE, 1, 1, '1 PRE' , channels)
% plotz(Pavg_1POST, 1, 1, '1 POST' , channels)
plotz(Pavg_2PRE, 1, 2, '2 PRE' , channels)
plotz(Pavg_2POST, 1, 2, '2 POST' , channels)
%%

% plotz(Gvar_1PRE, 1, 1, 'Gvar Sub 1 PRE' , channels)
plotz(Gvar_1POST, 1, 1, 'Gvar Sub 1 POST' , channels)
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
