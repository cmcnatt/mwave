classdef CollectionBase < handle
    
    properties (Access = protected)
        list;
        count;
    end
    
    methods (Access = protected)
        
        function [list, count] = getList(coll)
            % beware these are shallow copies!
            list = coll.list;
            count = coll.count;
        end
        
        function [] = add(coll, item, varargin)
            Nc = 1;
            if ~isempty(varargin)
                Nc = varargin{1};
            end
            
            if ~isscalar(Nc)
                error('The number of items must be scalar');
            end
                                    
            if isempty(coll.list)
                coll.list = cell(0);
                coll.count = [];
            end
            
            coll.list{end+1} = item;
            coll.count(end+1) = Nc;
        end
    end
end