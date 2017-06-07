classdef LookupTable < handle
    
    properties
        x;
        y;
    end
    
    methods
        
        function [table] = LookupTable(varargin)
            if length(varargin) > 1
                table.x = varargin{1};
                table.y = varargin{2};
            end
        end
        
        function [] = Add(table, xin, yin)
            N = length(xin);
            if length(yin) ~= N
                error('The length of xin and yin must be the same');
            end
            
            table.x(end+1:end+N) = xin;
            table.y(end+1:end+N) = yin;
        end
        
        function [val] = Get(table, xq)
            
            % first check if there are exact matches
            inds = table.x == xq;
            if sum(inds) > 1
                warning('More than one match');
            elseif sum(inds) == 1
                val = table.y(inds);
            else
                N = length(table.x);
                val = NaN;
                for n = 1:N
                    if xq <= table.x(n)
                        val = table.y(n);
                        break;
                    end
                end
            end
            
        end
    end
end