classdef IRandomValue < handle
    
    properties (SetAccess = protected, GetAccess = protected)
        mean;
        distType;
    end
    
    properties (Dependent)
       Mean;
       Std;
       DistType;
    end
    
    methods (Abstract)
        Sample(rnum, varargin);
    end
    
    methods (Abstract, Access = protected)   
        computeIfNot(hcomp);
    end
    
    methods        
        function [val] = get.Mean(rnum)
            % get the mean value
            val = rnum.mean;
        end
        
        function [val] = get.Std(rnum)
            % get the standard deviation
            val = rnum.std;
        end
        
        function [val] = get.DistType(rnum)
            % get the distribution type ('uniform' or 'normal')
            val = rnum.distType;
        end
        
        function [val] = Sample(rnum, varargin)
            a = 1;
            if ~isempty(varargin)
                a = varargin{1};
            end
            val = rnum.std*rand(a) + r
        end
    end
end