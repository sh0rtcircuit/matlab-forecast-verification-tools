classdef crpstest < matlab.unittest.TestCase
    
    properties
        probForecasts
        measurements
    end
    
    methods (TestMethodSetup)
        function loadDataset(obj)
            pd = makedist('Normal');
            obj.probForecasts = cell(3,1);
            obj.probForecasts(:) = {pd};
            obj.measurements = [-1 0 1]';
        end
    end
    
    methods (Test)
        function testUtilClasses(obj)
            val = prtools.eval.crps(obj.probForecasts,obj.measurements);
        end
    end
end

