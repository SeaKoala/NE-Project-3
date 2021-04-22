function [meanAlpha,stdAlpha, meanBeta, stdBeta] = DWT(signal)

    waveletFunction = 'db2';
    [C,L] = wavedec(signal,5,waveletFunction);
    cD11 = detcoef(C,L,1);                   %Gamma
    cD21 = detcoef(C,L,2);                   %Beta
    cD31 = detcoef(C,L,3);                   %Alpha
    cD41 = detcoef(C,L,4);                   %Theta
    cD51 = detcoef(C,L,5);                   
    cA51 = appcoef(C,L,waveletFunction,5);   %Delta
    D11 = wrcoef('d',C,L,waveletFunction,1); %Gamma
    D21 = wrcoef('d',C,L,waveletFunction,2); %Beta
    D31 = wrcoef('d',C,L,waveletFunction,3); %Alpha
    D41 = wrcoef('d',C,L,waveletFunction,4); %Theta
    D51 = wrcoef('d',C,L,waveletFunction,5); 
    A51 = wrcoef('a',C,L,waveletFunction,5); %Delta
        
    % Extract features: mean and std. of Alpha and Beta
    meanAlpha = mean(D31);
    stdAlpha = std(D31);
    
    meanBeta = mean(D21);
    stdBeta = std(D21);
    
    % also implement energy of alpha and beta
    
   
end