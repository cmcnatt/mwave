classdef IRandomVariable < matlab.mixin.Copyable
        
    properties (Abstract)
       Mean;
    end
    
    methods (Abstract)
        Sample(var, varargin);
    end
    
    methods
        
        function [val] = double(var)
            val = var.Mean;
        end
    end
    
    methods (Static)
        
        function [valOut] = TrySample(valIn, varargin)
            a = [];
            if ~isempty(varargin)
                if strcmp(varargin{1}, 'sample')
                    a = varargin{2};
                else
                    a = varargin{1};
                end
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