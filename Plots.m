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