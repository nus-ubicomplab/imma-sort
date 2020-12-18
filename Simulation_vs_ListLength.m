%% Simulation_vs_ListLength

% simulation settings
rounds = 100;
ListLen = 10 : 10 : 100;
coeff = -0.4;
i = 3;
paramNum = length(ListLen);
W = [0.7; 0.3];
attriNum = 2;

%%
cnt_total = paramNum*rounds;
simResults(cnt_total).Coeff = [];
simResults(cnt_total).Round = [];
simResults(cnt_total).Weight = [];
simResults(cnt_total).ListLen = [];
simResults(cnt_total).X = [];
simResults(cnt_total).controlPoints = [];
simResults(cnt_total).score = [];
simResults(cnt_total).PC = [];
i = 4;
r = 1;
% attriNum = 2;
cnt = 0;
for i = 1 : paramNum

    n = ListLen(i);
    disp(i)
    for r = 1 : rounds
        % generate data
        disp(r)
        cnt = cnt+1;
        [u, curcoeff] = random2dAttributeValues(coeff, n);
        X = u';
        alpha = ones(attriNum,1);
        if coeff < 0
            alpha(2) = -1;
        end
%         tic
        [xx, controlPoints, score] = PCA_Ranking(X,alpha, W, attriNum);
%         toc
        
        t = 0:0.0001:1;
        Z = [ones(1,length(t)); t; t.^2; t.^3];
        M = [1 -3 3 -1; 0, 3, -6, 3; 0, 0, 3, -3; 0, 0, 0, 1];
        bzc = controlPoints * M * Z;
        
        [mseA2, pc, s] = nearestPointsOnCurve(X, controlPoints, M, W);
        
        
        simResults(cnt).Coeff = coeff;
        simResults(cnt).Round = r;
        simResults(cnt).Weight = W;
        simResults(cnt).ListLen = n;
        simResults(cnt).X = X';
        simResults(cnt).controlPoints = controlPoints;
        simResults(cnt).score = score;
        simResults(cnt).PC = pc;

    end
end

% simResults_vs_listLen_v3 = simResults;
% save('simResults_vs_listLen_v3.mat', 'simResults_vs_listLen_v3')



%% calculate prediction intervals for the results

% load('simResults_vs_listLen_v3.mat')
% simResults = simResults_vs_listLen_v3;

simResults_vs_listLen_PICells_imma = cell(attriNum,length(ListLen));
simResults_vs_listLen_PICells_saw = cell(attriNum,length(ListLen));
simResults_vs_listLen_PICells_A1 = cell(attriNum,length(ListLen));
simResults_vs_listLen_PICells_A2 = cell(attriNum,length(ListLen));
i = 1;

for i = 1 : size(simResults,2)
    
    iX = simResults(i).X;
    iW = simResults(i).Weight;
    iListLen =  simResults(i).ListLen;
    iListLenIdx = find(ListLen == iListLen);
    
    iScore = simResults(i).score;
    [~, immaIdx] = sort(iScore);
    
    isaw = iX *iW;
    [~, isawIdx] = sort(isaw);
    
    [~, iA1Idx] = sort(iX(:,1));
    
    [~, iA2Idx] = sort(iX(:,2));
    
%     a = 1;
    for a = 1 : attriNum
            iattri = iX(:,a);
            
            iattri_sort = iattri(immaIdx);
            [iattri_PI_imma, f1_imma] = rankPredictionInterval(iattri_sort, ' ', 0.95);
            iattri_PI_avglength_imma = mean(iattri_PI_imma(:,2) - iattri_PI_imma(:,1));
            simResults_vs_listLen_PICells_imma{a,iListLenIdx} = [simResults_vs_listLen_PICells_imma{a,iListLenIdx}, iattri_PI_avglength_imma];
            
            % for saw
            isaw_sort = iattri(isawIdx);
            [iattri_PI_saw, ~] = rankPredictionInterval(isaw_sort, ' ', 0.95);
            iattri_PI_avglength_saw = mean(iattri_PI_saw(:,2) - iattri_PI_saw(:,1));
            simResults_vs_listLen_PICells_saw{a,iListLenIdx} = [simResults_vs_listLen_PICells_saw{a,iListLenIdx}, iattri_PI_avglength_saw];
            
            % for single sort
            isingle_sort = iattri(iA1Idx);
            [iattri_PI_single, ~] = rankPredictionInterval(isingle_sort, ' ', 0.95);
            iattri_PI_avglength_single = mean(iattri_PI_single(:,2) - iattri_PI_single(:,1));
            simResults_vs_listLen_PICells_A1{a,iListLenIdx} = [simResults_vs_listLen_PICells_A1{a,iListLenIdx}, iattri_PI_avglength_single];           
            
            % for single sort
            isingle_sort = iattri(iA2Idx);
            [iattri_PI_single, ~] = rankPredictionInterval(isingle_sort, ' ', 0.95);
            iattri_PI_avglength_single = mean(iattri_PI_single(:,2) - iattri_PI_single(:,1));
            simResults_vs_listLen_PICells_A2{a,iListLenIdx} = [simResults_vs_listLen_PICells_A2{a,iListLenIdx}, iattri_PI_avglength_single];           

    end
end