function [ meanQS, singleQS ] = quantileScore( probDists, measurements, varargin)
%quantileScore Computes the quantile score / pinball score / tick function.
%   INPUT:
%       probDists: N x 1 cell array with probability density class function
%           supports the <a
%           href="matlab:web('https://de.mathworks.com/help/stats/makedist.html')">Matlab distributions</a>.
%           If direct quantile forecasts are to be used, convert using
%           prtools.util.quantiles2ProbDist
%       measurements: N x 1 vector with "true" measurements.
%
%   PARAMETERS (Num-Value Pair):
%       'quantileValues': 1 x Q vector containing all desired quantile levels in interval [0..1].
%
%   OUTPUT:
%       meanQS: 1 x Q Mean quantile score for each
%       quantile.
%       singleQS: N x Q QS values for each fc-meas pair.
%
%   SEE prtools.util.quantiles2ProbDist

    % parse inputs
    p = inputParser;
    p.addRequired('probDists',@iscell);
    p.addRequired('measurements',@isvector);
    
    % all the quantile values to estimate
    p.addParameter('quantileValues',0.1:0.1:0.9,@isvector);
    
    % parse inputs
    p.parse(probDists,measurements,varargin{:});
    
    %get parsed variables
    probDists = p.Results.probDists;
    measurements = p.Results.measurements;
    quantileValues = reshape(p.Results.quantileValues,length(p.Results.quantileValues),1);

    % check input arguments and conform inputs
    [probDists, measurements] = prtools.util.checkConformScoreInput(probDists,measurements);
    
    % get quantile forecasts for all desired quantiles (quantileValues)
    quantileForecasts = cellfun(@(pd)pd.icdf(quantileValues),probDists,'UniformOutput',false);
    quantileForecasts = cell2mat(quantileForecasts')';
    
    % prepare data for matrix operation
    measRep = repmat(measurements,1,length(quantileValues));
    
    %compute the individual quantilescores
    singleQS = pPinball(measRep-quantileForecasts,quantileValues);
    
    %compute aggregated measure for each quantile
    meanQS = mean(singleQS,1);
end

function [ error ] = pPinball( distances, tau )
    
    %preallocate matrix
    error = zeros(size(distances));
    
    % repeat quantile
    if isscalar(tau)
       tau = repmat(tau,size(error));
    else
       tau = repmat(tau,size(error,1),1); 
    end
    
    %smaller 0
    error(distances<=0) = (tau(distances<=0) -1 ) .* distances(distances<=0);
    %bigger 0
    error(distances>0) = tau(distances>0) .* distances(distances>0);
end
