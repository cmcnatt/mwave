classdef IRandomVariable < matlab.mixin.Copyable
        
    properties (Abstract)
       Mean;
    end
    
    methods (Abstract)
        Sample(rnum, varargin);
    end
    
    methods
        
        function [val] = double(samp)
            val = samp.Mean;
        end
    end
    
    methods (Static)
        
        function [valOut] = TrySample(valIn, varargin)
            a = [];
            if ~isempty(varargin)
                a = varargin{1};
            end
            
            valOut = double(valIn);
            if ~isempty(a)
                if isa(valIn, 'IRandomVariable');
                    valOut = valIn.Sample(a);
                end
            end
        end
    end
end