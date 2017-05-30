classdef UniformValue < IRandomSample
    
    properties (SetAccess = protected, GetAccess = protected)
        lims;
    end
    
    properties (Dependent)
        Mean;
        Limits;
    end
    
    methods    
        
        function [rnum] = UniformValue(varargin)
            rnum.lims = [0 1];
            if ~isempty(varargin)
                rnum.lims = varargin{1};
            end
        end
        
        function [val] = get.Mean(rnum)
            % the mean value
            val = mean(rnum.lims);
        end
        
        function [val] = get.Limits(rnum)
            % the lower and upper bounds of the uniform distribution
            val = rnum.lims;
        end
        function [] = set.Limits(rnum, val)
            % the lower and upper bounds of the uniform distribution
            if length(val) > 2
                error('The limits must be a 1x2 or 2x1 vector');
            end
            if val(2) < val(1)
                error('The first value must be less than or equal to the second');
            end
            rnum.lims = val;
        end
        
        function [val] = Sample(rnum, varargin)
            % get a sample random number
            a = 1;
            if ~isempty(varargin)
                a = varargin{1};
            end

            low = rnum.lims(1);
            hi = rnum.lims(2);
            val = (hi-low).*rand(a) + low;
        end
    end
end