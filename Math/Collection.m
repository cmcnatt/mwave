classdef Collection < CollectionBase
    
    methods
        
        function [list, count] = GetList(coll)
            % beware these are shallow copies!
            [list, count] = coll.getList;
        end
        
        function [] = Add(coll, item, varargin)
            coll.add(item, varargin{:});
        end
    end
end