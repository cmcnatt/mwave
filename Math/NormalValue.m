classdef NormalValue < IRandomSample
    
    properties (SetAccess = protected, GetAccess = protected)
        mean;
        std;
    end
    
    properties (Dependent)
       Mean;
       Std;
    end
    
    methods    
        
        function [rnum] = NormalValue(varargin)
            rnum.mean = 0;
            rnum.std = 1;
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
        
        function [val] = Sample(rnum, varargin)
            % get a sample random number
            a = 1;
            if ~isempty(varargin)
                a = varargin{1};
            end
            
            val = rnum.std*randn(a) + rnum.mean;
        end
    end
end