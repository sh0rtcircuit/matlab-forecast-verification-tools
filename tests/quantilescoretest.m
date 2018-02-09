classdef quantilescoretest < matlab.unittest.TestCase
    
    properties
        probForecastsNorm
        crpsAnalytical
        measurements
    end
    
    methods (TestMethodSetup)
        function addnamespace(obj)
            addpath(fullfile('..'));
        end
        
        function loadDataset(obj)
            pd = makedist('Normal');
            obj.probForecastsNorm = cell(3,1);
            obj.probForecastsNorm(:) = {pd};
            obj.measurements = [-1 0 1]';
            
        end
    end
    
    methods (Test)
        function testQsComputation(obj)
            
            % check that mean computation is correct
            [meanqs, singleqs] = prtools.eval.quantileScore(obj.probForecastsNorm,obj.measurements);
            obj.assertEqual(meanqs,mean(singleqs,1));
            
            % test conversion 
            quantiles = 0.1:0.1:0.9;
            quantileForecasts = cellfun(@(pd)pd.icdf(quantiles),obj.probForecastsNorm,'UniformOutput',false);
            quantileForecastsMat = cell2mat(quantileForecasts);
            probFc = prtools.util.quantiles2ProbDist(quantileForecastsMat,quantiles,'lowerboundary',-10,'upperboundary',10);
            
            [meanqs2, singleqs2] = prtools.eval.quantileScore(obj.probForecastsNorm,obj.measurements);
            obj.assertEqual(meanqs,meanqs2);
            obj.assertEqual(singleqs,singleqs2);
        end
        
    end
end

