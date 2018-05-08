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
classdef ModesOfMotion < handle
    
    properties (Access = private)
        cg;
        vector;
        genMotFuncs;
    end

    properties (Dependent)
        Surge;
        Sway;
        Heave;
        Roll;
        Pitch;
        Yaw;
        Generalized;
        HasGen;
        NGen;
        Motions;        % Cell array of strings of active modes
        Cg;
        DoF;            % Number of degrees of freedom
        Vector;         % Vector of active modes
        MotionFuncs;    
        GenMotFuncs;
    end
    
    methods
         % Constructor 
        function [modes] = ModesOfMotion(mods)
            % Constructor
            if (nargin == 0)
                modes.cg = [0 0 0];
                modes.vector = [1 1 1 1 1 1];
                modes.genMotFuncs = [];
            else
                N = length(mods);
                if (N < 6)
                    error('Constructor input must be a vector of length 6 or greater');
                end
                
                modes.vector = zeros(1, N);
                
                for n = 1:6
                    if(~isBool(mods(n)))
                        error('The first 6 entries in the modes vector must be boolean');
                    end
                    modes.vector(n) = mods(n);
                end
                
                for n = 7:N
                    if (~isInt(mods(n)) || (mods(n) < 0))
                        error('The generalized modes entry must be a non-negative integer');
                    end
                    modes.vector(n) = mods(n);
                end
            end
        end
        
        function [val] = get.Cg(modes)
            val = modes.cg;
        end
        function [] = set.Cg(modes, val)
            modes.cg = val;
        end
        
        function [su] = get.Surge(modes)
            su = boolean(modes.vector(1));
        end
        function [modes] = set.Surge(modes, su)
            isBool(su);           
            modes.vector(1) = su;
        end
        
        function [sw] = get.Sway(modes)
            sw = boolean(modes.vector(2));
        end
        function [modes] = set.Sway(modes, sw)
            isBool(sw);           
            modes.vector(2) = sw;
        end
        
        function [he] = get.Heave(modes)
            he = boolean(modes.vector(3));
        end
        function [modes] = set.Heave(modes, he)
            isBool(he);           
            modes.vector(3) = he;
        end
        
        function [ro] = get.Roll(modes)
            ro = boolean(modes.vector(4));
        end
        function [modes] = set.Roll(modes, ro)
            isBool(ro);           
            modes.vector(4) = ro;
        end
        
        function [pit] = get.Pitch(modes)
            pit = boolean(modes.vector(5));
        end
        function [modes] = set.Pitch(modes, pit)
            isBool(pit);           
            modes.vector(5) = pit;
        end
        
        function [ya] = get.Yaw(modes)
            ya = boolean(modes.vector(6));
        end
        function [modes] = set.Yaw(modes, ya)
            isBool(ya);           
            modes.vector(6) = ya;
        end
        
        function [gen] = get.Generalized(modes)
            % The number of generalized modes
            if (modes.HasGen)
                gen = modes.vector(7:end);
            else
                gen = 0;
            end
        end
        function [modes] = set.Generalized(modes, gen)
            % The number of generalized modes
            if (~all(isInt(gen)) || all(gen < 0))
                error('The generalized modes entry must be a non-negative integer');
            end
            
            [M, N] = size(gen);
            if (M ~= 1)
                error('Generalized modes vector must be a 1xN vector');
            end
            
            modes.vector = [modes.vector(1:6), gen];
        end
        
        function [yn] = get.HasGen(modes)
            % Indicates whether it has generalized
            if (length(modes.vector) < 7)
                yn = false;
            else
                yn = true;
            end
        end
        
        function [ngen] = get.NGen(modes)
            % The number of generalized modes
            ngen = 0;
            N = length(modes.vector);
            for n = 7:N
                ngen = ngen + modes.vector(n);
            end
        end
        
        function [mo] = get.Motions(modes)
            % List of the mode names
            dof = modes.DoF;
            mo = cell(dof, 1);
            m = 1;
            for n = 1:6
                if (modes.vector(n))
                    switch n
                        case 1
                            mo{m} = 'Surge';
                        case 2
                            mo{m} = 'Sway';
                        case 3
                            mo{m} = 'Heave';
                        case 4
                            mo{m} = 'Roll';
                        case 5
                            mo{m} = 'Pitch';
                        case 6
                            mo{m} = 'Yaw';
                    end
                    m = m + 1;
                end
            end
            for n = m:dof
                mo{n} = ['Gen' num2str(n-m+1)];
            end
        end
        
        function [dof] = get.DoF(modes)
            % Number of degrees of freedom
            dof = sum(modes.vector);
        end
        
        function [v] = get.Vector(modes)
            % Mode Vector
            v = modes.vector;
        end
        function [modes] = set.Vector(modes, v)         
            % Mode Vector
            N = length(v);
            if(N < 6)
                error('Modes vector must be at least of length 6');
            end
            
            for n = 1:6
                if(~isBool(v(n)))
                    error('The first 6 entries in the modes vector must be boolean');
                end
                modes.vector(n) = v(n);
            end
            
            for n = 7:N
                if (~isInt(v(n)) || (v(n) < 0))
                    error('The generalized modes entry must be a non-negative integer');
                end
                modes.vector(n) = v(n);
            end
        end
        
        function [mfuncs] = get.MotionFuncs(modes)
            mind = 1;
            for n = 1:6
                if (modes.vector(n))
                    switch(n)
                        case 1
                            mfuncs(mind) = SurgeFunc;
                        case 2
                            mfuncs(mind) = SwayFunc;
                        case 3
                            mfuncs(mind) = HeaveFunc;
                        case 4
                            mfuncs(mind) = RollFunc;
                        case 5
                            mfuncs(mind) = PitchFunc;
                        case 6
                            mfuncs(mind) = YawFunc;
                    end
                    mfuncs(mind).Cg = modes.cg;
                    mind = mind + 1;
                end
            end
            
            for n = 7:length(modes.vector)
                nGenn = modes.vector(n);
                for m = 1:nGenn
                    mfuncs(mind) = modes.genMotFuncs(n-6);
                    mind = mind + 1;
                end
            end
        end
        
        function [mfuncs] = get.GenMotFuncs(modes)
            mfuncs = modes.genMotFuncs;
        end
        function [modes] = set.GenMotFuncs(modes, funcs)
            Nf = length(funcs);
            if (Nf ~= (length(modes.vector) - 6))
                error('The number of functions must equal the number of different types of generalized modes');
            end
            
            for n = 1:Nf
                if(~isa(funcs(n), 'IMotionFunc'))
                    error('All functions must be IMotionFunc');
                end
            end
            
            modes.genMotFuncs = funcs;
        end
    end
        
end