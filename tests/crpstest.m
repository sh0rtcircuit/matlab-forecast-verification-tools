classdef crpstest < matlab.unittest.TestCase
    
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
            
            % compare with closed form solution
            % from https://cran.r-project.org/web/packages/scoringRules/vignettes/crpsformulas.html#Normal
            crpsClosedForm = @(pd,meas) meas .* (2.* pd.cdf(meas) -1) + 2.* pd.pdf(meas) - 1./sqrt(pi);
            obj.crpsAnalytical = cellfun(crpsClosedForm,obj.probForecastsNorm,num2cell(obj.measurements));
        end
    end
    
    methods (Test)
        function testCrpsComputation(obj)
            
            [~,singleCrps] = prtools.eval.crps(obj.probForecastsNorm,obj.measurements);
            %error smaller than 1%
            obj.assertLessThan(obj.crpsAnalytical-singleCrps,0.01.*(obj.crpsAnalytical));
        end
        
        function testCrpsParallel(obj)
            if license('test', 'distrib_computing_toolbox') % only if parallel computing toolbox available
                [~,singleCrps] = prtools.eval.crps(obj.probForecastsNorm,obj.measurements,'useParallel',true);
                obj.assertLessThan(obj.crpsAnalytical-singleCrps,0.01.*(obj.crpsAnalytical));
            end
        end
        
        function testCrpsTrapz(obj)
            [~,singleCrps] = prtools.eval.crps(obj.probForecastsNorm,obj.measurements,'trapz',true);
            % error less than 1%
            obj.assertLessThan((obj.crpsAnalytical-singleCrps)./obj.crpsAnalytical,0.01);
            
            [~,singleCrps] = prtools.eval.crps(obj.probForecastsNorm,obj.measurements,'trapz',true,'trapzSampling',10000);
            % error now less than 0.1%
            obj.assertLessThan((obj.crpsAnalytical-singleCrps)./obj.crpsAnalytical,0.01);
            
            [~,singleCrps] = prtools.eval.crps(obj.probForecastsNorm,obj.measurements,'trapz',true,'trapzSampling',1000,'trapzMin',-30,'trapzMax',30);
            % now larger interval, thus higher errow (error now less than
            % 3%)
            obj.assertLessThan((obj.crpsAnalytical-singleCrps)./obj.crpsAnalytical,0.03);
        end
    end
end

