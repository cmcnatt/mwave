classdef RandomValue < IRandomSample
    
    properties (SetAccess = protected, GetAccess = protected)
        mean;
        std;
        distType;
    end
    
    properties (Dependent)
       Mean;
       Std;
       DistType;
    end
    
    methods    
        
        function [rnum] = RandomValue(varargin)
            rnum.mean = 0;
            rnum.std = 1;
            rnum.distType = 'normal';
            if ~isempty(varargin)
                rnum.mean = varargin{1};
                rnum.std = varargin{2};
                if length(varargin) > 2
                    rnum.distType = varargin{3};
                end
            end
        end
        
        function [val] = get.Mean(rnum)
            % the mean value
            val = rnum.mean;
        end
        function [] = set.Mean(rnum, val)
            % the mean value
            if ~isscalar(val)
                error('Mean must be a scalar value');
            end
            rnum.mean = val;
        end
        
        function [val] = get.Std(rnum)
            % the standard deviation
            val = rnum.std;
        end
        function [] = set.Std(rnum, val)
            % the mean value
            if ~isscalar(val)
                error('Std must be a scalar value');
            end
            rnum.std = val;
        end
        
        function [val] = get.DistType(rnum)
            % the distribution type ('uniform' or 'normal')
            val = rnum.distType;
        end
        function [] = set.DistType(rnum, val)
            % the distribution type ('uniform' or 'normal')
            if ~strcmp(val, 'normal')
                if ~strcmp(val, 'uniform')
                    error('DistTyhpe must be a ''normal'' or ''uniform''');
                end
            end
            rnum.distType = val;
        end
        
        function [val] = Sample(rnum, varargin)
            % get a sample random number
            a = 1;
            if ~isempty(varargin)
                a = varargin{1};
            end
            if strcmp(rnum.distType, 'normal')
                val = rnum.std*randn(a) + rnum.mean;
            else
                low = rnum.mean - rnum.std/2;
                hi = rnum.mean + rnum.std/2;
                val = (hi-low).*rand(a) + low;
            end
        end
    end
end