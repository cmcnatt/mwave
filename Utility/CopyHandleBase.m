%{ 
mwave - A water wave and wave energy converter computation package 
Copyright (C) 2014  Cameron McNatt

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contributors:
    C. McNatt
%}
classdef CopyHandleBase < handle
    
    methods (Abstract, Access = protected)
        setAnyProp(obj, prop, val);
        % Implement this function in the class:
        %   function [] = setAnyProp(obj, prop, val)
        %       obj.(prop) = val;
        %   end
    end
    
    methods
        function [newObj] = copy(thisObj)
            % Instantiate new object of the same class.
            newObj = feval(class(thisObj));

            % Create a metaclass object
            mc = metaclass(thisObj);
            p = mc.PropertyList;
            
            warning('off');
            thisObjStruct = struct(thisObj);
            warning('on');
            
            for n = 1:length(p)
                if ~p(n).NonCopyable
                    newObj.setAnyProp(p(n).Name, thisObjStruct.(p(n).Name));
                end
            end
        end
    end
end