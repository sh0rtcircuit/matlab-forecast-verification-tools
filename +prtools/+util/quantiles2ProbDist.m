function [probDists] = quantiles2ProbDist(quantileForecasts,quantiles,varargin)
%QUANTILES2PROBDIST Summary of this function goes here
%   Detailed explanation goes here

    % parse inputs
    p = inputParser;
    p.addRequired('quantileForecasts',@ismatrix);
    p.addRequired('quantiles',@isvector);
    p.addOptional('lowerboundary',0,@isscalar);
    p.addOptional('upperboundary',1,@isscalar);
    
    % parse inputs
    p.parse(quantileForecasts,quantiles,varargin{:});
    
    %get parsed variables
    quantileForecasts = p.Results.quantileForecasts;
    quantiles = p.Results.quantiles;
    lowerboundary = p.Results.lowerboundary;
    upperboundary = p.Results.upperboundary;
    
    %check and conform inputs
    if min(quantiles) < 0
        error(['Minimum quantile value has to be 0.']);
    elseif max(quantiles) > 1
        error(['Maximum quantile value has to be 1.']);
    end
    if size(quantileForecasts,2) ~= length(quantiles)
        error(['Error: Number of quantile forecasts has to equal numbero f quantiles.']);
    end
    quantiles = reshape(quantiles,1,[]);
    
    % sort quantiles ascending
    [quantiles, permutIdx] = sort(quantiles);
    % sort the corresponding forecasts accordingly
    quantileForecasts = quantileForecasts(:,permutIdx);
    
    % add boundaries if necessary
    if min(quantiles) ~= 0
        quantiles = [0 quantiles];
        quantileForecasts = [repmat(lowerboundary,size(quantileForecasts,1),1) quantileForecasts];
    end
    if max(quantiles) ~= 1
        quantiles = [quantiles 1];
        quantileForecasts = [quantileForecasts repmat(upperboundary,size(quantileForecasts,1),1)];
    end
    
    % define conversion function
    convertFunc = @(quantileForecast,quantile) ...
        makedist('PiecewiseLinear','x',quantileForecast','Fx',quantile');
    
    %convert to probdist structure
    quantileForecastsCell = mat2cell(quantileForecasts,...
        ones(size(quantileForecasts,1),1), ...
        size(quantileForecasts,2));
    quantilesCell{1} = quantiles;
    quantilesCell = repmat(quantilesCell,size(quantileForecastsCell,1),1);
    probDists = cellfun(convertFunc,quantileForecastsCell,...
        quantilesCell);
end

