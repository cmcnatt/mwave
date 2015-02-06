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
classdef FloatingBody < matlab.mixin.Heterogeneous & handle
    % Defines a wave energy converter geometry
    
    properties (Access = protected)
        handle;
        geoFile;
        panelGeo;
        cg;
        m;
        d;
        k;
        position;
        wpSec;
        angle;
        modes;
        nGen;
        iGenMds;
        iLowHi;
        iSurfPan;
        writeFileMeth;
        writeParams;
    end

    properties (Dependent)
        Handle;         % Descriptive name of floating body
        GeoFile;        % Name of .cfg geometry file (do not include ".cfg")
        PanelGeo;       % The actual panel geometry  
        Cg;             % Centor of Gravity in body coordinates
        M;              % Mass matrix
        Dpto;           % PTO Damping matrix
        K;              % Stiffness matrix
        XYpos;          % XY Position of body origin in global coordinates 
        Zpos;           % Z position of body origin in global coordinates - seperate from the XYPosition because it controls the point about which the body pitches
        WaterPlaneSec;  % (x, y) coordinates of the cross-section of the body at the water plane in the bodies coordinate system
        Rcir;           % Radius of circumscribed cirle
        Angle;          % Rotation angle (deg) of body x-axis relative to global x-axis
        Modes;          % Modes to be evaluated on this floating body
        Ngen;           % The number of generalized modes;
        WamIGenMds;     % Integer that specifies which NEWMODES subroutine to call.  All bodies in a given run must have the same value.
        WamILowHi;      % Indicates whether the geometry file is low-order (0) or high-order (1) (for WAMIT)
        ISurfPan;       % Indicates whether the geometry has interior surface panels
        WriteFileMeth;  % Additional methods to support the writing of a geometry file
        WriteParams;    % Parameters of the WriteFileMeth
    end
 
    methods

        function [fb] = FloatingBody(varargin)
            % Constructor
            if (nargin == 0)
                fb.handle = 'newFB';
                fb.geoFile = [];
                fb.panelGeo = [];
                fb.cg = zeros(1, 3);
                fb.m = zeros(6, 6);
                fb.d = zeros(6, 6);
                fb.k = zeros(6, 6);
                fb.position = zeros(1, 3);
                fb.angle = 0;
                fb.wpSec = [];
                fb.modes = ModesOfMotion;
                fb.nGen = 0;
                fb.iGenMds = 0;
                fb.iLowHi = 0;
                fb.iSurfPan = 0;
                fb.writeFileMeth = [];
                fb.writeParams = [];
            elseif (isa(varargin{1},'FloatingBody'))
                fbin = varargin{1};
                fb.geoFile = fbin.geoFile;
                fb.panelGeo = fbin.panelGeo;
                fb.cg = fbin.cg;
                fb.m = fbin.m;
                fb.d = fbin.d;
                fb.k = fbin.k;
                fb.position = fbin.position;
                fb.angle = fbin.angle;
                fb.wpSec = fbin.wpSec;
                fb.modes = fbin.modes;
                fb.nGen = fbin.nGen;
                fb.iGenMds = fbin.iGenMds;
                fb.iLowHi = fbin.iLowHi;
                fb.iSurfPan = fbin.iSurfPan;
                fb.writeFileMeth = fbin.writeFileMeth;
                fb.writeParams = fbin.writeParams;
            else
                error('Input not allowed');
            end
        end
        
        function [han] = get.Handle(fb)
            % Get the handle (descriptive name) of the floating body
            han = fb.handle;
        end
        function [fb] = set.Handle(fb, han)
            % Set the handle (descriptive name) of the floating body
            if (ischar(han))
                fb.handle = han;
            else
                error('The Handle must be a string');
            end
        end

        function [fn] = get.GeoFile(fb)
            % Get the name of geometry file associated with this floating body
            fn = fb.geoFile;
        end
        function [fb] = set.GeoFile(fb, fn)
            % Set the name of geometry file associated with this floating body
            if (~isempty(fb.panelGeo))
                error('The FloatingBody contains a panel geometry that defines the geometry, and so an external geometry file cannot be used');
            else
                if (ischar(fn))
                    fb.geoFile = fn;
                    lfn = length(fn);
                    if(lfn > 4)
                        ending = fn((lfn-3):lfn);
                        if (ending == '.gdf')
                            error('The GeoFile should not include .gdf');
                        end
                    end
                else
                    error('The GeoFile must be a string');
                end
            end
        end
        
        function [pg] = get.PanelGeo(fb)
            % Get the geometry (PanelGeo) object associated with this floating body
            pg = fb.panelGeo;
        end
        function [fb] = set.PanelGeo(fb, panGeo)
            % Set the geometry (PanelGeo) object associated with this floating body
            if (~isa(panGeo, 'PanelGeo'))
                error('Input must be a PanelGeo');
            end
            fb.panelGeo = panGeo;
        end
        
        function [C_g] = get.Cg(fb)
            % Get the center of gravity of the floating body
            C_g = fb.cg;
        end
        function [fb] = set.Cg(fb, cg)
            % Set the center of gravity of the floating body
            fb.checkSizeNx1(cg,3);           
            fb.cg = cg;
        end
        
        function [m_] = get.M(fb)
            % Get the mass matrix
            inds = 6 + fb.nGen;
            m_ = fb.m(1:inds, 1:inds);
        end
        function [fb] = set.M(fb, m_)
            % Set the mass matrix
            fb.checkSize(m_);            
            fb.m = m_;
        end
       
        function [d_] = get.Dpto(fb)
            % Get the PTO damping matrix
            inds = 6 + fb.nGen;
            d_ = fb.d(1:inds, 1:inds);
        end
        function [fb] = set.Dpto(fb, d_)
            % Set the PTO damping matrix
            fb.checkSize(d_);            
            fb.d = d_;
        end
        
        function [k_] = get.K(fb)
            % Get the mechanical stiffness matrix
            inds = 6 + fb.nGen;
            k_ = fb.k(1:inds,1:inds);
        end
        function [fb] = set.K(fb, k_)
            % Get the mechanical stiffness matrix
            fb.checkSize(k_);            
            fb.k = k_;
        end
        
        function [p] = get.XYpos(fb)
            % Get the x-y position of the body origin in global coordinates 
            p = fb.position(1:2);
        end
        function [fb] = set.XYpos(fb, p)
            % Set the x-y position of the body origin in global coordinates 
            fb.checkSizeNx1(p,2);
            for n = 1:2
                fb.position(n) = p(n);
            end
        end
        
        function [p] = get.Zpos(fb)
            % Get the z position of the body origin in global coordinates 
            p = fb.position(3);
        end
        function [fb] = set.Zpos(fb, p)
            % Set the z position of the body origin in global coordinates 
            fb.position(3) = p;
        end
        
        function [wp] = get.WaterPlaneSec(fb)
            % Get the (x, y) coordinates of the cross-section of the body at the water plane in the body's coordinate system
            wp = fb.wpSec;
        end
        function [fb] = set.WaterPlaneSec(fb, wp)
            % Set the (x, y) coordinates of the cross-section of the body at the water plane in the body's coordinate system
            [row, col] = size(wp);
            if (row < 3)
                error('Water plane section must have at least 3 points');
            end
            
            if (col ~= 2)
                error('Water plane section must be a Nx2 array of x-y points');
            end
            
            if ((wp(1,1) ~= wp(row, 1)) || (wp(1,2) ~= wp(row,2)))
                error('First point in section must be the same at the last to close the section');
            end
            
            fb.wpSec = wp;
        end
        
        function [r] = get.Rcir(fb)
            % Get the radius of the circumscribing circle.
            R = sqrt(fb.wpSec(:,1).^2 + fb.wpSec(:,2).^2);
            r = max(R);
        end
        
        function [a] = get.Angle(fb)
            % Get the angle of body rotation from body coordinates to
            % global coordinates
            a = fb.angle;
        end
        function [fb] = set.Angle(fb, a)
            % Set the angle of body rotation from body coordinates to
            % global coordinates
            if (a < 0 || a > 360)
                error('The rotation angle must be between 0 and 360 degrees');
            end
            fb.angle = a;
        end
        
        function [mo] = get.Modes(fb)
            % Get the modes of motion of the floating body
            mo = fb.modes;
        end
        function [fb] = set.Modes(fb, mo)
            % Set the modes of motion of the floating body
            if (isa(mo, 'ModesOfMotion'))
                fb.modes = mo;
                fb.nGen = mo.NGen;
                [n n] = size(fb.m);
                if (n < (6 + fb.nGen))
                    m_ = zeros(6 + fb.nGen, 6 + fb.nGen);
                    m_(1:n, 1:n) = fb.m;
                    fb.m = m_;
                    d_ = zeros(6 + fb.nGen, 6 + fb.nGen);
                    d_(1:n, 1:n) = fb.d;
                    fb.d = d_;
                    k_ = zeros(6 + fb.nGen, 6 + fb.nGen);
                    k_(1:n, 1:n) = fb.k;
                    fb.k = k_;
                end
            else
                error('The Modes must be of type ModesOfMotion');
            end
        end
        
        function [ngen] = get.Ngen(fb)
            % Get the number of generalized modes of the floating body
            ngen = fb.nGen;
        end
        
        function [igm] = get.WamIGenMds(fb)
            % Get the Wamit IGENMDS value, which is used for determining
            % which generalized modes fortran routine is run in Wamit.
            igm = fb.iGenMds;
        end
        function [fb] = set.WamIGenMds(fb, igm)
            % Set the Wamit IGENMDS value, which is used for determining
            % which generalized modes fortran routine is run in Wamit.
            if (~isInt(igm) || (igm < 1))
                error('value must be a positive integer');
            end
            fb.iGenMds = igm;
        end

        function [ilh] = get.WamILowHi(fb)
            % Get the Wamit ILOWHI value, which is used for indicating
            % whether the body geometry file is low-order or higher order.
            ilh = fb.iLowHi;
        end
        function [fb] = set.WamILowHi(fb, ilh)
            % Set the Wamit ILOWHI value, which is used for indicating
            % whether the body geometry file is low-order or higher order.
            if (~isBool(ilh))
                error('value must be boolean, 1 or 0');
            end
            fb.iLowHi = ilh;
        end
        
        function [isp] = get.ISurfPan(fb)
            % Indicates whether the geometry has interior surface panels
            isp = fb.iSurfPan;
        end
        function [fb] = set.ISurfPan(fb, isp)
            % Indicates whether the geometry has interior surface panels
            if (~isBool(isp))
                error('value must be boolean, 1 or 0');
            end
            fb.iSurfPan = isp;
        end
        
        function [wfm] = get.WriteFileMeth(fb)
            % Additional methods to support the writing of a geometry file
            wfm = fb.writeFileMeth;
        end
        
        function [wpar] = get.WriteParams(fb)
            % Additional methods to support the writing of a geometry file
            wpar = fb.writeParams;
        end
    end
   
    methods (Access = protected)
        function [] = checkSize(fb, arg)
            fail = 1;
            if (ndims(arg) == 2)
                [ro co] = size(arg);
                if (ro == co)
                    fail = 0;
                end
                N = 6 + fb.nGen;
                if (ro == N)
                    fail = 0;
                end
            end
            if (fail)
                error(['Value must be a ' num2str(N) 'x' num2str(N) ' matrix, '...
                    'where ' num2str(N) ' = 6 + ' num2str(fb.nGen) ', and ' num2str(fb.nGen) ' is the number of generalized modes.'...
                    'If you intend generalized modes to be used, the Modes property with the generalized modes must be set before this matrix.']);
            end
        end
        
        function [] = checkSizeNx1(fb, arg, N)
            fail = 1;
            if (ndims(arg) == 2)
                [ro co] = size(arg);
                if (((ro == N) && (co == 1)) || ((ro == 1) && (co == N)))
                    fail = 0;
                end
            end
            if (fail)
                error(['Value must be a ' num2str(N) 'x1 or a 1x' num2str(N) ' vector']);
            end
        end
    end 
end