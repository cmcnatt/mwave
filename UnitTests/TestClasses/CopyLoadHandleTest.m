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
classdef CopyLoadHandleTest < CopyLoadHandleBase
    
    properties (Access = public)
        PubProp1;
        PubProp2;
%         PubProp3;   % delete
    end
    
    properties (Access = private)
        priProp1;
        priProp2;
%         priProp3;   % delete
        unSetProps;
    end
    
    properties (Access = protected)
        proProp1;
        proProp2;
%         proProp3;   % delete
    end
    
    properties (Dependent)
        PriProp1;
        ProProp1;
        % to delete
%         PriProp3;
%         ProProp3;
        % end to delete
    end
    
    methods
        
        function [this] = CopyLoadHandleTest(varargin)
            if ~isempty(varargin)
                this.priProp2 = varargin{1};
                this.proProp2 = varargin{2};
            end
        end
        
        function [val] = get.PriProp1(this)
            val = this.priProp1;
        end
        function [] = set.PriProp1(this, val)
            this.priProp1 = val;
        end
        
        function [val] = get.ProProp1(this)
            val = this.proProp1;
        end
        function [] = set.ProProp1(this, val)
            this.proProp1 = val;
        end
        
        function [val] = GetPriProp2(this)
            val = this.priProp2;
        end
        
        function [] = SetPriProp2(this, val)
            this.priProp2 = val;
        end
        
        function [val] = GetProProp2(this)
            val = this.proProp2;
        end
        
        function [] = SetProProp2(this, val)
            this.proProp2 = val;
        end
        
        function [val] = GetUnsetProps(this)
            val = this.unSetProps;
        end
        
        % To delete (comment)
%         function [val] = get.PriProp3(this)
%             val = this.priProp3;
%         end
%         function [] = set.PriProp3(this, val)
%             this.priProp3 = val;
%         end
%         
%         function [val] = get.ProProp3(this)
%             val = this.proProp3;
%         end
%         function [] = set.ProProp3(this, val)
%             this.proProp3 = val;
%         end
        % end to delete
    end
        
    methods (Access = protected)
        function [] = setAnyProp(obj, prop, val)
            obj.(prop) = val;
        end
    end
    
    methods (Static)
        function [newObj] = loadobj(a)
            if isstruct(a)
                newObj = feval(mfilename('class'));
                
                [~, loadedUnset] = CopyLoadHandleBase.trySetProps(newObj, a);
                
                newObj.unSetProps = loadedUnset;
            else
                newObj = a;
            end
        end
    end
end