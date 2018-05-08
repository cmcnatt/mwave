classdef NormalVariable < IRandomVariable
    
    properties (SetAccess = protected, GetAccess = protected)
        lims;
        mean;
        std;
    end
    
    properties (Dependent)
       Mean;
       Std;
       Limits;
    end
    
    methods    
        
        function [rnum] = NormalVariable(mean, std, lims)
            if nargin < 1
                rnum.mean = 0;
            else
                rnum.mean = mean;
            end
            if nargin < 2
                rnum.std = 1;
            else
                rnum.std = std;
            end
            if nargin < 3
                rnum.lims = [-Inf, Inf];
            else
                rnum.lims = lims;
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
        
        function [val] = get.Limits(rnum)
            % limits on the distribution
            val = rnum.lims;
        end
        function [] = set.Limits(rnum, val)
            % limits on the distribution
            if length(val) ~= 2
                error('The limits must be a vector of a lower and upper limit');
            end
            rnum.lims = val;
        end
        
        function [val] = Sample(rnum, varargin)
            % get a sample random number
            a = 1;
            if ~isempty(varargin)
                a = varargin{1};
            end
            
            val = rnum.std*randn(a) + rnum.mean;
            
            val(val < rnum.lims(1)) = rnum.lims(1);
            val(val > rnum.lims(2)) = rnum.lims(2);
        end
    end
end