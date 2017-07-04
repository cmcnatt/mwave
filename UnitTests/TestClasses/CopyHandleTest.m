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
classdef CopyHandleTest < CopyHandleBase
    
    properties (Access = public)
        PubProp1;
        PubProp2;
    end
    
    properties (Access = private)
        priProp1;
        priProp2;
    end
    
    properties (Access = protected)
        proProp1;
        proProp2;
    end
    
    properties (Dependent)
        PriProp1;
        ProProp1;
    end
    
    methods
        
        function [this] = CopyHandleTest(varargin)
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
    end
        
    methods (Access = protected)
        function [] = setAnyProp(obj, prop, val)
            obj.(prop) = val;
        end
    end
end