function [ meanCrps, singleCrps ] = crps( probDists, measurements ,varargin)
%crps Compute the Continuous ranked probability score.
%   INPUT:
%       probDists: N x 1 cell array with probability density class function
%           supports the distributions of https://de.mathworks.com/help/stats/makedist.html
%       measurements: N x 1 vector with "true" measurements
%   OUTPUT:
%       meanCrps: The mean CRPS over all evaluated prob density and meas.
%           pairs
%       singleCrps: N x 1 CRPS values for each pair.
%

    % parse inputs
    p = inputParser;
    p.addRequired('probDists',@iscell);
    p.addRequired('measurements',@isvector);
    
    % use parallel processing (only active if 'trapz'=false)
    p.addParameter('useParallel',false,@isscalar);
    
    % use faster (but less precise) trapz integration
    p.addParameter('trapz',false,@isscalar);
    
    % no. of sampling points for trapz integration
    p.addParameter('trapzSampling',1000,@isscalar);
    
    % estimation of sensible minimum and maximum for trapz integration
    p.addParameter('trapzMin',min(cellfun(@(pdf)pdf.icdf(eps),probDists)),@isscalar);
    p.addParameter('trapzMax',max(cellfun(@(pdf)pdf.icdf(1-eps),probDists)),@isscalar);
    
    % parse inputs
    p.parse(probDists,measurements,varargin{:});
    
    %get parsed variables
    probDists = p.Results.probDists;
    measurements = p.Results.measurements;
    useParallel = p.Results.useParallel;
    trapzIntegration = p.Results.trapz;
    trapzSampling = p.Results.trapzSampling;
    trapzMin = p.Results.trapzMin;
    trapzMax = p.Results.trapzMax;

    % check input arguments and conform inputs
    [probDists, measurements] = prtools.util.checkConformScoreInput(probDists,measurements);
    
    %compute the CRPS    
    if ~trapzIntegration
        % initialize resultstore
        singleCrps = nan(size(measurements));
        
        % init workers
        if useParallel workers=inf;, else workers=0;, end
        
        %compute individual probfc-meas pairs
        parfor (i=1:numel(probDists),workers)
            singleCrps(i) = computeCRPS(probDists{i},measurements(i));
        end
    else
        %compute invidivual crps values
        singleCrps = computeCRPSAll(probDists,measurements,trapzSampling,trapzMin,trapzMax);
    end
    
    %build mean crps score
    meanCrps = mean(singleCrps);
end

% this is the precise (but slower) variant for computation
function res = computeCRPS(pdf,measurement)

    %integrate from -inf to current measurement
    fixCdfL = @(x) cdf(pdf,x).^2;
    S1 = integral(fixCdfL,-inf,measurement);
    
    %integrate from measurement to inf
    fixCdfU = @(x) (1-cdf(pdf,x)).^2;
    S2 = integral(fixCdfU,measurement,inf);
    
    %return sum
    res = (S1 + S2);
end

% this is a more approximative (but faster version)
function singleCrps = computeCRPSAll(pdfs,measurements,noSamplingPoints,trapzMin,trapzMax)
    
    % create sampling points, and replicate for matrix multiplication
    samplingPts = linspace(trapzMin,trapzMax,noSamplingPoints);
    samplingPtsRepMat = repmat(samplingPts,numel(pdfs),1);
    samplingPtsRepCell = mat2cell(samplingPtsRepMat,ones(numel(pdfs),1),size(samplingPtsRepMat,2));
    
    % abstract integration function
    f = @(singlepdf,x) singlepdf.cdf(x);
    
    % get all cdf values & convert to matrix
    cdfCell = cellfun(f,pdfs,samplingPtsRepCell,'UniformOutput',false);
    cdfMat = cell2mat(cdfCell);
    
    % compute CRPS components (up to measurement and above)
    crpsOverall = cdfMat.^2;
    crpsOverall(samplingPtsRepMat>measurements) = (1-cdfMat(samplingPtsRepMat>measurements)).^2;
    
    % compute integrals and store
    singleCrps = trapz(samplingPts,crpsOverall')';    
end
