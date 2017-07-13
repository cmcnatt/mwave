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
classdef CopyLoadHandleBase < handle
    
    properties (Access = protected)
        unloadedProps;
    end
    
    methods (Abstract, Access = protected)
        setAnyProp(obj, prop, val);
        % Implement this function in the class:
        %   function [] = setAnyProp(obj, prop, val)
        %       obj.(prop) = val;
        %   end
    end
    
    methods (Abstract, Static)
        loadobj(a)
        % Implement something like this in the class:
        % function [newObj] = loadobj(a)
        %     if isstruct(a)
        %         newObj = feval(class(mfilename('class')));
        %         
        %         [currentUnset, loadedUnset] = CopyLoadHandleBase.trySetProps(newObj, a);
        %     else
        %         newObj = a;
        %     end
        % end
    end
    
    methods
        function [newObj] = copy(thisObj)
            % Instantiate new object of the same class.
            newObj = feval(class(thisObj));
            
            warning('off');
            thisObjStruct = struct(thisObj);
            warning('on');
            
            CopyLoadHandleBase.trySetProps(newObj, thisObjStruct);
        end
    end
    
    methods (Static, Access = protected)
        function [currentUnset, loadedUnset] = trySetProps(newObj, astruct)
            % Create a metaclass object
            mc = metaclass(newObj);
            p = mc.PropertyList;
            fields = fieldnames(astruct);
            
            currentUnset = struct;
            currn = 0;
            
            for n = 1:length(p)
                if ~p(n).NonCopyable
                    if isfield(astruct, p(n).Name)
                        newObj.setAnyProp(p(n).Name, astruct.(p(n).Name));
                        astruct.(p(n).Name) = 'valueloaded';
                    else
                        currn = currn + 1;
                        currentUnset.(p(n).Name) = [];
                    end
                else
                    astruct.(p(n).Name) = 'valueloaded';
                end
            end
            
            loadedUnset = struct;
            
            for n = 1:length(fields)
                if ~strcmp(astruct.(fields{n}), 'valueloaded')
                    loadedUnset.(fields{n}) = astruct.(fields{n});
                end
            end
        end
    end
end