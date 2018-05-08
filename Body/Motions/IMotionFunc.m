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
classdef IMotionFunc < matlab.mixin.Heterogeneous & handle
    % Defines motions for modes of motion
    
    properties(Access = protected)
        cg;
        isym;
        isGen;
        evalCheckFunc;
    end
    
    properties (Dependent)
        Cg;
        ISym;   % Indicates the type of Symmetry of the motion
        IsGen;
        EvalCheckFunc;
    end
    
    properties (Abstract)
        MotionIn;
    end
    
    methods (Abstract)
        Evaluate(pos);
        Divergence(pos);
        GravityForce(pos);
    end
    
    methods
        function [c] = get.Cg(motF)
            c = motF.cg;
        end
        function [] = set.Cg(motF, c)
            [N, M] = size(c);
            
            if ((N ~= 1) || (M ~= 3))
                error('Cg must be a 1x3 vector');
            end
            
            motF.cg = c;
        end
        
        function [isy] = get.ISym(motF)
            % Indicates the type of Symmetry of the motion. 
            % 1, 5 - Surge, pitch like symmetry
            % 2, 4 - Sway, roll like symmetry
            % 3 - Heave like symmetry
            % 6 - Yaw like symmetry
            isy = motF.isym;
        end
        
        function [isg] = get.IsGen(motF)
            isg = motF.isGen;
        end
        
        function [ef] = get.EvalCheckFunc(motF)
            ef = motF.evalCheckFunc;
        end
        function [motF] = set.EvalCheckFunc(motF, ef)
            motF.evalCheckFunc = ef;
        end
    end
    
    methods (Access = protected)
        function [] = initCg(motF, varargin)
            if (isempty(varargin))
                motF.cg = [0 0 0];
            elseif (length(varargin) == 1)
                c = varargin{1};
                
                [N, M] = size(c);
            
                if ((N ~= 1) || (M ~= 3))
                    error('Cg must be a 1x3 vector');
                end
                
                motF.cg = c;
            else
                error('MotionFunc input not recognize');               
            end
        end
        
        function [eval] = checkEvalPoint(motF, pos)
            if (isempty(motF.evalCheckFunc))
                eval = true;
            else
                eval = motF.evalCheckFunc(pos);
            end
        end
    end
end