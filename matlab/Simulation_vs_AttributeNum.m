%% Simulation Imma vs AttributeNum
rounds = 100; % number of round for each set of parameters
AttriNum = 2 : 1 : 6; % attribute number
paramNum = length(AttriNum);

% generate weights
Weights = cell( length(AttriNum), 1);
weightRatios = 1: 6;%[1, 2, 4, 8, 16, 32];
i = 6;
for i = 1 : length(AttriNum)
iw = weightRatios(1 : AttriNum(i));
iw_norm = flip( iw /(sum(iw)));
Weights{i} = iw_norm';
end


cnt_total = paramNum*rounds;
simResults(cnt_total).Coeff = []; % storing result for each run
simResults(cnt_total).Round = [];
simResults(cnt_total).Weight = [];
simResults(cnt_total).AttriNum = [];
simResults(cnt_total).X = [];
simResults(cnt_total).controlPoints = [];
simResults(cnt_total).score = [];
simResults(cnt_total).PC = [];
i = 4;
r = 1;
n = 100;
coeff = 0.4;
cnt = 0;
for i = 1 : paramNum
    attriNum = AttriNum(i);
%     disp(i)
    for r = 1 : rounds
        % generate data
        disp(r)
        cnt = cnt+1;
        [u, curcoeff] = randomMoreThan2AttributeValues(coeff, n, attriNum); % randomly generate samples with > 2 attributes with fixed coeff
        
        alpha = ones(attriNum,1); % all the attributes are positive correlation in this simulation
        W = Weights{i};
        X = u';
%         tic
        [~, controlPoints, score] = PCA_Ranking(X,alpha, W, attriNum); % run PCA algorithm
%         toc
        
        t = 0:0.0001:1;
        Z = [ones(1,length(t)); t; t.^2; t.^3];
        M = [1 -3 3 -1; 0, 3, -6, 3; 0, 0, 3, -3; 0, 0, 0, 1];
        [msedist, pc, s] = nearestPointsOnCurve(X, controlPoints, M, W); % analysis result details
               
        simResults(cnt).Coeff = curcoeff;
        simResults(cnt).Round = r;
        simResults(cnt).Weight = W;
        simResults(cnt).AttriNum = attriNum;
        simResults(cnt).X = X';
        simResults(cnt).controlPoints = controlPoints;
        simResults(cnt).score = score;
        simResults(cnt).PC = pc;

    end
end


%% calculate PI for each results:

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
