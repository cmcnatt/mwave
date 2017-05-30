classdef IRandomSample < matlab.mixin.Copyable
        
    properties (Abstract)
       Mean;
    end
    
    methods (Abstract)
        Sample(rnum, varargin);
    end
end