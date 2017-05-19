classdef IRandomSample < handle
        
    properties (Abstract)
       Mean;
       Std;
    end
    
    methods (Abstract)
        Sample(rnum, varargin);
    end
end