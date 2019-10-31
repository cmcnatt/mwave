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
        compGeoFile;
        panelGeo;
        stlGeoFile;
        stlGeoFileCg;
        stlGeo;
        cg;
        cb;
        wetVol;
        totVol;
        surfArea;
        massBall;
        m;
        dpto;
        dpar;
        c;
        hasc;
        k;
        visc;
        position;
        centRot;
        wpSec;
        angle;
        modes;
        nGen;
        iGenMds;
        iLowHi;
        iSurfPan;
        makeSurfPan;
        wdipole;
        panSize;
        surfAboveZ0;
        writeFileMeth;
        writeParams;
    end

    properties (Dependent)
        Handle;         % Descriptive name of floating body
        GeoFile;        % Name of .gdf geometry file (do not include ".gdf")
        CompGeoFile;    % Name of .gdf geometry file associated with the computation of pessures and velocities on the body surface
        PanelGeo;       % The actual panel geometry  
        StlGeoFile;     % Name and file location of .stl file of geomety
        StlGeoFileCg;   % The CG of the Stl geometry file
        StlGeo;         % The actualy stl geometry object
        Cg;             % Center of Gravity in body coordinates
        CgGlobal;       % Center of Gravity in Global coordinates
        Cb;             % Center of Buoyancy in body coordinates
        WetVolume;      % Submerged volume
        TotalVolume;    % Total volume
        SurfArea;       % Total surface area
        MassBallast;    % Estimated mass of the ballast required for a given submergence
        M;              % Mass matrix
        Dpto;           % PTO Damping matrix
        Dpar;           % Parasitic damping matrix 
        C;              % Externally computed hydrostatic matrix
        K;              % Stiffness matrix
        ViscDampCoefs;  % Array of the viscous damping coefficient
        XYpos;          % XY Position of body origin in global coordinates 
        Zpos;           % Z position of body origin in global coordinates - seperate from the XYPosition because it controls the point about which the body pitches
        CenterRot;      % Positon about which body rotated in body coordinates - if not specified, then is CG.
        WaterPlaneSec;  % (x, y) coordinates of the cross-section of the body at the water plane in the bodies coordinate system
        WaterPlaneSecGlobal; % water plane section in global coords
        Rcir;           % Radius of circumscribed cirle
        Angle;          % Rotation angle (deg) of body x-axis relative to global x-axis
        Modes;          % Modes to be evaluated on this floating body
        Ngen;           % The number of generalized modes;
        WamIGenMds;     % Integer that specifies which NEWMODES subroutine to call.  All bodies in a given run must have the same value.
        WamILowHi;      % Indicates whether the geometry file is low-order (0) or high-order (1) (for WAMIT)
        WamPanelSize;   % See WAMIT PANEL_SIZE .cfg parameter
        WamDipoles;     % Indices of panels or patches that are thin members (dipoles)
        ISurfPan;       % Indicates whether the geometry has interior surface panels
        MakeSurfPan;    % Indicates whether WAMIT should make interior surface panels - overrides ISurfPan
        SurfAboveZ0;    % Indicates that there are surfaces of the geometry above the z = 0 water plane
        WriteFileMeth;  % Additional methods to support the writing of a geometry file
        WriteParams;    % Parameters of the WriteFileMeth
    end
 
    methods

        function [fb] = FloatingBody(varargin)
            % Constructor
            if (nargin == 0)
                fb.handle = 'newFB';
                fb.geoFile = [];
                fb.compGeoFile = [];
                fb.panelGeo = [];
                fb.stlGeoFile = [];
                fb.stlGeoFileCg = [0 0 0];
                fb.stlGeo = [];
                fb.cg = zeros(1, 3);
                fb.cb = zeros(1, 3);
                fb.m = zeros(6, 6);
                fb.dpto = zeros(6, 6);
                fb.dpar = zeros(6, 6);
                fb.k = zeros(6, 6);
                fb.c = zeros(6, 6);
                fb.visc = [];
                fb.hasc = false;
                fb.position = zeros(1, 3);
                fb.centRot = [];
                fb.angle = 0;
                fb.wpSec = [];
                fb.modes = ModesOfMotion;
                fb.nGen = 0;
                fb.iGenMds = 0;
                fb.iLowHi = 0;
                fb.iSurfPan = 0;
                fb.makeSurfPan = 0;
                fb.wdipole = [];
                fb.surfAboveZ0 = false;
                fb.writeFileMeth = [];
                fb.writeParams = [];
                fb.massBall = [];
                fb.wetVol = [];
                fb.totVol = [];
                fb.surfArea = [];
            elseif (isa(varargin{1},'FloatingBody'))
                fbin = varargin{1};
                fb.geoFile = fbin.geoFile;
                fb.compGeoFile = fbin.compGeoFile;
                fb.panelGeo = fbin.panelGeo;
                fb.stlGeoFile = fbin.stlGeoFile;
                fb.stlGeoFileCg = fbin.stlGeoFileCg;
                fb.stlGeo = fbin.stlGeo;
                fb.cg = fbin.cg;
                fb.cb = fbin.cb;
                fb.m = fbin.m;
                fb.dpto = fbin.dpto;
                fb.dpar = fbin.dpar;
                fb.k = fbin.k;
                fb.c = fbin.c;
                fb.visc = fbin.visc;
                fb.hasc = fbin.hasc;
                fb.position = fbin.position;
                fb.centRot = fbin.centRot;
                fb.angle = fbin.angle;
                fb.wpSec = fbin.wpSec;
                fb.modes = fbin.modes;
                fb.nGen = fbin.nGen;
                fb.iGenMds = fbin.iGenMds;
                fb.iLowHi = fbin.iLowHi;
                fb.panSize = fbin.panSize;
                fb.iSurfPan = fbin.iSurfPan;
                fb.makeSurfPan = fbin.makeSurfPan;
                fb.wdipole = fbin.wdipole;
                fb.surfAboveZ0 = fbin.surfAboveZ0;
                fb.writeFileMeth = fbin.writeFileMeth;
                fb.writeParams = fbin.writeParams;
                fb.massBall = fbin.massBall;
                fb.wetVol = fbin.wetVol;
                fb.totVol = fbin.totVol;
                fb.surfArea = fbin.surfArea;
            else
                error('Input not allowed');
            end
        end
        
        function [han] = get.Handle(fb)
            % Get the handle (descriptive name) of the floating body
            han = fb.handle;
        end
        function [] = set.Handle(fb, han)
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
        function [] = set.GeoFile(fb, fn)
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
        
        function [fn] = get.CompGeoFile(fb)
            % Get the name of geometry file associated with the compuation of pointns on the floating body
            fn = fb.compGeoFile;
        end
        function [] = set.CompGeoFile(fb, fn)
            % Get the name of geometry file associated with the compuation of pointns on the floating body
            if (ischar(fn))
                fb.compGeoFile = fn;
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
                
        function [pg] = get.PanelGeo(fb)
            % Get the geometry (PanelGeo) object associated with this floating body
            pg = fb.panelGeo;
        end
        function [] = set.PanelGeo(fb, panGeo)
            % Set the geometry (PanelGeo) object associated with this floating body
            if (~isa(panGeo, 'PanelGeo'))
                error('Input must be a PanelGeo');
            end
            fb.panelGeo = panGeo;
        end
        
        function [fn] = get.StlGeoFile(fb)
            % Name and file location of .stl file of geomety
            fn = fb.stlGeoFile;
        end
        function [] = set.StlGeoFile(fb, val)
            % Name and file location of .stl file of geomety
            if (ischar(val))
                fb.stlGeoFile = val;
            else
                error('The StlGeoFile must be a string');
            end
        end
        
        function [fn] = get.StlGeoFileCg(fb)
            % The CG of the Stl Geometry file
            fn = fb.stlGeoFileCg;
        end
        function [] = set.StlGeoFileCg(fb, val)
            % The CG of the Stl Geometry file
            fb.stlGeoFileCg = val;
        end
        
        function [val] = get.StlGeo(fb)
            % The actualy stl geometry object
            if isempty(fb.stlGeo)
                if ~isempty(fb.stlGeoFile)
                    stl = Stl6DOFGeo;
                    stl.Read(fb.stlGeoFile, fb.stlGeoFileCg);
                    stl.Cg = fb.cg;
                    fb.stlGeo = stl;
                end
            end

            val = fb.stlGeo;
        end
        function [] = set.StlGeo(fb, val)
            % The actualy stl geometry object
            if (~isa(val, 'Stl6DOFGeo'))
                error('Input must be a Stl6DOFGeo');
            end
            fb.stlGeo = val;
        end
        
        function [C_g] = get.Cg(fb)
            % Get the center of gravity of the floating body
            C_g = fb.cg;
        end
        function [] = set.Cg(fb, cg)
            % Set the center of gravity of the floating body
            fb.checkSizeNx1(cg,3);        
            fb.onModifyCg(cg);
            fb.cg = cg;
            fb.setCgOnModes;
        end
        
        function [cgg] = get.CgGlobal(fb)
            % Get the center of gravity of the floating body in Global
            % coords
            cgb = fb.cg;
            
            v = fb.position;
            cgg = cgb + v;
        end
        
        function [C_b] = get.Cb(fb)
            % Get the center of buoyancy of the floating body
            C_b = fb.cb;
        end
        function [] = set.Cb(fb, cb)
            % Set the center of buoyancy of the floating body
            fb.checkSizeNx1(cb,3);           
            fb.cb = cb;
        end
        
        function [val] = get.WetVolume(fb)
            % Get the submerged volume of the floating body
            val = fb.wetVol;
        end
        function [] = set.WetVolume(fb, val)
            % Set the submerged volume of the floating body
            fb.checkSizeNx1(val, 1);        
            fb.wetVol = val;
        end
        
        function [val] = get.TotalVolume(fb)
            % Get the total volume of the floating body
            val = fb.totVol;
        end
        function [] = set.TotalVolume(fb, val)
            % Set the submerged volume of the floating body
            fb.checkSizeNx1(val, 1);        
            fb.totVol = val;
        end
        
        function [val] = get.SurfArea(fb)
            % Get the total surface area of the floating body
            val = fb.surfArea;
        end
        function [] = set.SurfArea(fb, val)
            % Set the total surface area of the floating body
            fb.checkSizeNx1(val, 1);        
            fb.surfArea = val;
        end
        
        function [val] = get.MassBallast(fb)
            % Get the estimated mass of the ballast required for a given submergence
            val = fb.massBall;
        end
        function [] = set.MassBallast(fb, val)
            % Set the estimated mass of the ballast required for a given submergence
            fb.checkSizeNx1(val, 1);        
            fb.massBall = val;
        end
        
        function [m_] = get.M(fb)
            % Get the mass matrix
            inds = 6 + fb.nGen;
            m_ = fb.m(1:inds, 1:inds);
        end
        function [] = set.M(fb, m_)
            % Set the mass matrix
            fb.checkSize(m_);            
            fb.m = m_;
        end
       
        function [d_] = get.Dpto(fb)
            % Get the PTO damping matrix
            inds = 6 + fb.nGen;
            d_ = fb.dpto(1:inds, 1:inds);
        end
        function [] = set.Dpto(fb, d_)
            % Set the PTO damping matrix
            fb.checkSize(d_);            
            fb.dpto = d_;
        end
        
        function [d_] = get.Dpar(fb)
            % Get the parasitic damping matrix
            inds = 6 + fb.nGen;
            d_ = fb.dpar(1:inds, 1:inds);
        end
        function [] = set.Dpar(fb, d_)
            % Set the parasitic damping matrix
            fb.checkSize(d_);            
            fb.dpar = d_;
        end
        
        function [c_] = get.C(fb)
            % Get the hydrostatic stiffness matrix
            if (fb.hasc)
%                 inds = 6 + fb.nGen;
%                 c_ = fb.c(1:inds,1:inds);
                c_ = fb.c;
            else
                c_ = [];
            end
        end
        function [] = set.C(fb, c_)
            % Set the hydrostatic stiffness matrix
            fb.checkSize(c_);            
            fb.c = c_;
            fb.hasc = true;
        end
        
        function [k_] = get.K(fb)
            % Get the mechanical stiffness matrix
            inds = 6 + fb.nGen;
            k_ = fb.k(1:inds,1:inds);
        end
        function [] = set.K(fb, k_)
            % Get the mechanical stiffness matrix
            fb.checkSize(k_);            
            fb.k = k_;
        end
        
        function [val] = get.ViscDampCoefs(fb)
            % The array of viscous damping coefficients
            val = fb.visc;
        end
        function [] = set.ViscDampCoefs(fb, val)
            % The array of viscous damping coefficients
            if ~isa(val, 'ViscDampCoef')
                error('ViscDampCoefs must be of type ViscDampCoef');
            end
            fb.visc = val;
        end
        
        function [p] = get.XYpos(fb)
            % Get the x-y position of the body origin in global coordinates 
            p = fb.position(1:2);
        end
        function [fb] = set.XYpos(fb, p)
            % Set the x-y position of the body origin in global coordinates 
            fb.checkSizeNx1(p,2);
            v = [0 0 0];
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
            v = [0 0 p];
            fb.onModifyPos(v);
            fb.position(3) = p;
        end
        
        function [p] = get.CenterRot(fb)
            % Positon about which body rotated in body coordinates - if not specified, then is CG.
            if (isempty(fb.centRot))
                p = fb.cg;
            else
                p = fb.centRot;
            end
        end
        function [fb] = set.CenterRot(fb, p)
            % Positon about which body rotated in body coordinates - if not specified, then is CG.
            fb.checkSizeNx1(p,3);
            fb.centRot = p;
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
        
        function [wpSecG] = get.WaterPlaneSecGlobal(fb)
            % move the waterplane section
            vxy = fb.position(1:2);
            
            wpSecG = fb.wpSec;
            N = size(fb.wpSec, 1);
            for n = 1:N
                wpSecG = wpSecG + vxy;
            end
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
                    d_(1:n, 1:n) = fb.dpto;
                    fb.dpto = d_;
                    k_ = zeros(6 + fb.nGen, 6 + fb.nGen);
                    k_(1:n, 1:n) = fb.k;
                    fb.k = k_;
                end
                fb.setCgOnModes;
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
        
        function [ilh] = get.WamPanelSize(fb)
            ilh = fb.panSize;
        end
        function [fb] = set.WamPanelSize(fb, ps)
            fb.panSize = ps;
        end
        
        
        function [isp] = get.ISurfPan(fb)
            % Indicates whether the geometry has interior surface panels
            isp = fb.iSurfPan;
        end
        function [] = set.ISurfPan(fb, isp)
            % Indicates whether the geometry has interior surface panels
            if (~isBool(isp))
                error('value must be boolean, 1 or 0');
            end
            fb.iSurfPan = isp;
        end
        
        function [val] = get.MakeSurfPan(fb)
            % Indicates whether WAMIT should make interior surface panels
            % Overrides ISurfPan
            val = fb.makeSurfPan;
        end
        function [] = set.MakeSurfPan(fb, val)
            % Indicates whether WAMIT should make interior surface panels
            % Overrides ISurfPan
            if (~isBool(val))
                error('value must be boolean, 1 or 0');
            end
            fb.makeSurfPan = val;
        end
                
        function [val] = get.WamDipoles(fb)
            % Indices of the panels or patches that are thin members (i.e.
            % dipoles)
            val = fb.wdipole;
        end
        function [] = set.WamDipoles(fb, val)
           % Indices of the panels or patches that are thin members (i.e.
            % dipoles)
            if (~isInt(val))
                error('values must be integers');
            end
            fb.wdipole = val;
        end
        
        function [val] = get.SurfAboveZ0(fb)
            val = fb.surfAboveZ0;
        end
        function [] = set.SurfAboveZ0(fb, val)
            if ~isBool(val)
                error('SurfAboveZ0 must be boolean');
            end
            fb.surfAboveZ0 = val;
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
    
    methods (Static)
        function [] = WriteBodyMetaData(folder, name, bodies, deltaCg, addParam)
            
            if nargin < 4
                deltaCg = [0 0 0];
            end
            if nargin < 5
                addParam = [];
            end
            
            fileName = [folder '\' name '.csv'];
            Nbod = length(bodies);
            
            fid = fopen(fileName, 'w+');
            fprintf(fid, 'parameter,units,');
            for n = 1:Nbod
                fprintf(fid, 'body %i,', n);
            end
            fprintf(fid, '\n');
            
            fprintf(fid, 'name,,');
            for n = 1:Nbod
                fprintf(fid, '%s,', bodies(n).Handle);
            end
            fprintf(fid, '\n');
            
            fprintf(fid, 'mass,kg,');
            for n = 1:Nbod
                fprintf(fid, '%6.2f,', bodies(n).M(1,1));
            end
            fprintf(fid, '\n');
            
            names = {'Ixx', 'Iyy', 'Izz', 'Ixy', 'Ixz', 'Iyz'};
            ind1 = [4, 5, 6, 4, 4, 5];
            ind2 = [4, 5, 6, 5, 6, 6];
            
            for m = 1:6
                fprintf(fid, '%s,kg*m2,', names{m});
                for n = 1:Nbod
                    fprintf(fid, '%6.2f,', bodies(n).M(ind1(m), ind2(m)));
                end
                fprintf(fid, '\n');
            end
            
            names = {'Cg-x', 'Cg-y', 'Cg-z'};
            
            for m = 1:3
                fprintf(fid, '%s,m,', names{m});
                for n = 1:Nbod
                    fprintf(fid, '%6.2f,', bodies(n).Cg(m) + deltaCg(m));
                end
                fprintf(fid, '\n');
            end
            
            if ~isempty(addParam)
                for m = 1:length(addParam)
                    N = length(addParam(m).names);
                    for n = 1:N
                        fprintf(fid, '%s,%s,', addParam(m).names{n}, addParam(m).units{n});
                        for o = 1:Nbod
                            val = addParam(m).values(n, o);
                            if addParam(m).incDelta(n)
                                val = val + deltaCg(n);
                            end
                            fprintf(fid, '%6.2f,', val);
                        end
                        fprintf(fid, '\n');
                    end
                end
            end
            
            fclose(fid);
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
        
        function [] = onModifyCg(fb, cg)
            % do nothing in the superclass
        end
        
        function [] = setCgOnModes(fb)
            if ~isempty(fb.modes)
                fb.modes.Cg = fb.cg;
            end
        end
        
        function [] = onModifyPos(fb, v)
            % do nothing in the superclass
        end
    end 
end