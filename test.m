s = sub2PRE_DATA(:,19,32);

figure;
plot(s)

waveletFunction = 'db4';
                [C,L] = wavedec(s,6,waveletFunction);
                
                cD1 = detcoef(C,L,1); % NOISY
                cD2 = detcoef(C,L,2); % NOISY
                cD3 = detcoef(C,L,3); %GAMA
                cD4 = detcoef(C,L,4); %BETA
                cD5 = detcoef(C,L,5); %ALPHA
                cD6 = detcoef(C,L,6); %THETA
                cA6 = appcoef(C,L,waveletFunction,6); %DELTA
                
                
                D1 = wrcoef('d',C,L,waveletFunction,1); % NOISY
                D2 = wrcoef('d',C,L,waveletFunction,2); % NOISY
                D3 = wrcoef('d',C,L,waveletFunction,3); %GAMMA
                D4 = wrcoef('d',C,L,waveletFunction,4); %BETA
                D5 = wrcoef('d',C,L,waveletFunction,5); %ALPHA
                D6 = wrcoef('d',C,L,waveletFunction,6); %THETA
                A6 = wrcoef('a',C,L,waveletFunction,6); %DELTA
                
                Gamma = D3;
                figure; subplot(5,1,1); plot(1:1:length(Gamma),Gamma);title('GAMMA');
               
                Beta = D4;
                subplot(5,1,2); plot(1:1:length(Beta), Beta); title('BETA');
                betaEnergy = sum(abs(Beta).^2);
                
                Alpha = D5;
                subplot(5,1,3); plot(1:1:length(Alpha),Alpha); title('ALPHA');
                alphaEnergy = sum(abs(Alpha).^2);
                
                Theta = D6;
                subplot(5,1,4); plot(1:1:length(Theta),Theta);title('THETA');
                
                Delta = A6;
                subplot(5,1,5);plot(1:1:length(Delta),Delta);title('DELTA');
                
                % Energy 
                totalE = alphaEnergy + betaEnergy + sum(abs(Gamma).^2) +sum(abs(Theta).^2) + sum(abs(Delta).^2) + sum(abs(D1).^2) + sum(abs(D2).^2);
                
       
                
                
                
                