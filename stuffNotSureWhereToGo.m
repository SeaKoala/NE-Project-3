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


%% GAVG Plots
c = 1
figure
subplot(2,2,1)
plot(GAVG_1_PRE(:, c, 1));
subplot(2,2,2)
plot(GAVG_1_POST(:, c, 1));
subplot(2,2,3)
plot(GAVG_2_PRE(:, c, 1));
subplot(2,2,4)
plot(GAVG_2_POST(:, c, 1));

%%

figure('units','normalized','Position',[0.2,0.65,0.3,0.3])
h = spectrum.welch; % creates the Welch spectrum estimator
% 2 sec segment length -> 1024 samples
% more frequency resolution at higher sample rate 
h.SegmentLength = 512;
SOIf3=psd(h,sub1PRE_DATA(:, c, 25),'Fs',fs); % calculates and plot the one sided PSD
plot(SOIf3); % Plot the one-sided PSD. 
temp =get(gca);
temp.Children(1).Color = 'r';
hold on;
SOIf3=psd(h,sub1PRE_DATA(:, c, 50),'Fs',fs); % calculates and plot the one sided PSD
plot(SOIf3); % Plot the one-sided PSD. 
temp =get(gca);
temp.Children(1).Color = 'g';
hold on;
SOIf3=psd(h,sub1PRE_DATA(:, c, 75),'Fs',fs); % calculates and plot the one sided PSD
plot(SOIf3); % Plot the one-sided PSD. 
temp =get(gca);
temp.Children(1).Color = 'b';
legend('flex','ext','rest');
title("Trial 1 PSD Flex vs Rest Single Channel Sub1 PRE")'

%% Feature Extraction
x = sub1PRE_DATA;
size(x);
var = zeros(32, 75);
kurt = zeros(32, 75);
skew = zeros(32, 75);
for i = 1:32
    for j = 1:75
        var(i,j) = jfeeg('var', x(:,i,j));
        kurt(i,j) = jfeeg('kurt', x(:,i,j));
        skew(i,j) = jfeeg('skew', x(:,i,j));

    end
end

%%
c = 9;
figure('units','normalized','Position',[0.2,0.65,0.3,0.3])
hold on
plot(var(c,:));
plot(kurt(c,:));
plot(skew(c,:));
hold off
legend( 'var', 'kurt', 'skew');


%% SNR
mse = zeros(32,25);
mae = zeros(32,25);
SNR = zeros(32,25);
PSNR = zeros(32,25);


for c = 1:32
    for j = 1:25
        temp = sub1PRE_DATA(:, c, j); % signal
        y = sub1PRE_DATA(:, c, j+50); % rest

        mse(c,j)=0;
        for i=1:length(temp)
        mse(c,j)=mse(c,j)+(y(i)-temp(i))^2;
        end
        mse(c,j)=mse(c,j)/length(temp);
        %MAE %Mean absolute error
        mae(c,j)=0;
        for i=1:length(temp)
        mae(c,j)=mae(c,j)+abs(y(i)-temp(i));
        end
        mae(c,j)=mae(c,j)/length(temp);

        num=0;
        den=0;
        for i=1:length(temp)
        den=den+(y(i)-temp(i))^2;
        end
        for i=1:length (temp)
        num=num+temp(i)^2;
        end
        SNR(c,j) = 20*log10(sqrt(num)/sqrt(den));
        PSNR(c, j)= 20*log10(max(temp)/sqrt(mse(c,j)));

    end
end

%%
figure
sgtitle('Filt EXT vs REST');
for i =1:32
    subplot(8,4,i)
    hold on
%     plot(mse(i,:));
%     plot(mae(i,:));
    plot(SNR(i,:));
    plot(PSNR(i,:));
    hold off
end
% legend('mse', 'mae', 'SNR', 'PSNR');
legend('SNR', 'PSNR');
    
%%
clc
channels(13) = {'T9'};
channels(19) = {'T10'};
for i = 1:25
    subplot(5,5,i)
    plot_topography(channels, mse(:,i));
end
    
    
    
    

   
    

 
    
    
    
    
    
    
    
    
    
    