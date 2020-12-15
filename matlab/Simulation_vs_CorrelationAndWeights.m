%% Simulation_vs_Correlations & Weights
% clear all;
coeffs = -1 : 0.2 : 1;
coeffs(1) = -0.9999;
coeffs(end) = 0.9999;
rounds = 100;
W = [0.5;0.5];
AttriNum = 2;
n = 100;   


r = -1:0.2:1;
weightRatio = 10.^(r);
w1 = weightRatio ./ (1 + weightRatio);
w2 = 1 ./ (1 + weightRatio);
Weights = [w1;w2];
Weights = Weights(:,1:6);

 M = [1 -3 3 -1; 0, 3, -6, 3; 0, 0, 3, -3; 0, 0, 0, 1];
 
 %%

coeffNum = length(coeffs);
weightNum = size(Weights,2);

cnt_total = weightNum * coeffNum*rounds ;
simResults(cnt_total).Coeff = [];
simResults(cnt_total).Weight = [];
simResults(cnt_total).Round = [];
simResults(cnt_total).X = [];
simResults(cnt_total).controlPoints = [];
simResults(cnt_total).score = [];
simResults(cnt_total).PC = [];

cnt = 0;
% c = 6;
for c = 1 : coeffNum
    coeff = coeffs(c);
    for w = 1 : weightNum
        W = Weights(:,w);
        disp(['coeff = ', num2str(coeff), ' weight1 = ', num2str(W(1))])
        for r = 1 : rounds
            % generate data
            cnt = cnt+1;
            [u, ~] = random2dAttributeValues(coeff, n);
            X = u';
            
            alpha = [1, 1]';
            if coeff < 0
                alpha(2) = -1;
            end       
            
            tic
            [~, controlPoints, score] = PCA_Ranking(X,alpha, W, AttriNum);
            toc
            
            [~, pc, s] = nearestPointsOnCurve(X, controlPoints, M, W);
            
            simResults(cnt).Coeff = coeff;
            simResults(cnt).Round = r;
            simResults(cnt).X = X';
            simResults(cnt).controlPoints = controlPoints;
            simResults(cnt).score = score;
            simResults(cnt).PC = pc;
        end
    end
end

% simResults_vs_CorrelationsAndWeights = simResults;
% save('simResults_vs_CorrelationsAndWeights.mat', 'simResults_vs_CorrelationsAndWeights')


%% calculate Prediction Intervels for Imma Sort results

simResults_vs_attributeNum_PICells = cell(5,6);
for i = 1 : size(simResults,2)
    
    iX = simResults(i).X;
    iW = simResults(i).Weight;
    iAttriNum = simResults(i).AttriNum;
    
    iScore = simResults(i).score;
    [~, immaIdx] = sort(iScore);
    a = 1;
    for a = 1 : max(AttriNum)
        if a <= iAttriNum
            iattri = iX(:,a);
            iattri_sort = iattri(immaIdx);
            [iattri_PI_imma, f1_imma] = rankPredictionInterval(iattri_sort, ' ', 0.95);
            iattri_PI_avglength_imma = mean(iattri_PI_imma(:,2) - iattri_PI_imma(:,1));
            simResults_vs_attributeNum_PICells{iAttriNum-1,a} = [simResults_vs_attributeNum_PICells{iAttriNum-1,a}, iattri_PI_avglength_imma];
        end
    end
end
