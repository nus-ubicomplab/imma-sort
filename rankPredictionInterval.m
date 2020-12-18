function [actual_PI, f1, ci, rmse] = rankPredictionInterval(attributeValueSorted, isNormalize, confidencelevel)

% attributeValueSorted = X_dist(1,:)';
% isNormalize = '';

% confidencelevel = 0.95;
listLen = length(attributeValueSorted);
rank = 1 : listLen;
rank = rank';

if strcmp(isNormalize,'normalize')
    attributeValueSorted = (attributeValueSorted - min(attributeValueSorted)) / (max(attributeValueSorted) - min(attributeValueSorted));
end

f1=fit(rank,attributeValueSorted,'poly2');
y_preds = f1(rank);

% prediction interval
se_sqrt = 1 / listLen + (rank - mean(rank)).^2 ./ (sum((rank - mean(rank)).^2));
mse =   sum((attributeValueSorted - y_preds).^2) / listLen;
percentile = tinv(confidencelevel,(listLen-2));
PI = percentile * sqrt(mse * ( 1 + se_sqrt) );
actual_PI = [y_preds - PI, y_preds + PI];

% confidence interval
ci = predint(f1,rank, confidencelevel, 'functional','off');

% rmse
sse = sum((attributeValueSorted - y_preds).^2);
rmse = sqrt(sse/listLen);

% ci_1 = predint(f1,rank, 0.95,'observation','off');

% figure;
% plot(f1,rank,price_imma); hold on;
% plot(rank, ci_1, 'm--'); hold on;
% plot(rank, actual_PI, 'k--'); hold on;
