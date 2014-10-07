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
classdef BodySurfWaveField < FBWaveField
    
    properties (Access = private)
        bodyGeo;
    end
    
    properties (Dependent)
        BodyGeo;
    end
    
    methods
        function [wf] = BodySurfWaveField(body, iwave, swave, varargin)
            
            if (~isa(body, 'PanelGeo'))
                error('Body must be a PanelGeo');
            end
            if (iwave.IsArray)
                error('BodySurfWaveField must be a collection of points');
            end
            
            bodyPnts = body.Centroids;
            wfPnts = iwave.FieldPoints;
            
            Npnts = size(bodyPnts, 1);
            
            if (Npnts ~= size(wfPnts, 1))
                error('The body points and wave field points must be the same');
            end
            
            for n = 1:Npnts
                if any(abs(bodyPnts(n,:) - wfPnts(n,:)) > 1e-6*[1 1 1])
   %                  error('The body points and wave field points must be the same');
                end
            end
            
            wf = wf@FBWaveField(iwave, swave, varargin{:});
            
            wf.bodyGeo = body;
        end
        
        function [bod] = get.BodyGeo(wf)
            bod = wf.bodyGeo;
        end
        
        % Pressure
        function [pgeo] = Pressure(wf, type)
            p = wf.pressure(type);
            pgeo = cell(size(p));
            N = numel(p);
            
            for n = 1:N
                pgeon = PanelGeo(wf.bodyGeo);
                pgeon.Values = p{n};
                pgeo{n} = pgeon;
            end
        end
                
        % Elevation
        function [eta] = Elevation(wf, type)
            error('Elevation not implemented for BodySurfWaveField');
        end
        
        % Velocity
        function [velgeo] = Velocity(wf, type)
            vel = wf.velocity(type);
            velgeo = cell(size(vel));
            N = numel(vel);
            
            for n = 1:N
                velgeon = PanelGeo(wf.bodyGeo);
                velgeon.Values = vel{n};
                velgeo{n} = velgeon;
            end
        end
        
        % SignificantWaveHeight
        function [hs] = SigWaveHeight(wf, type)
            error('SigWaveHeight not implemented for BodySurfWaveField');
        end
        
        % Gets spectrums at the points closest to the desired points
        function [specs, actPoints] = Spectra(wf, varargin)
            error('Spectra not implemented for BodySurfWaveField');
        end
        
        % EnergyFlux
        function [fluxgeo] = EnergyFlux(wf, surf, varargin)
            fluxgeo = wf.bodyGeo;
            
            pnts = fluxgeo.Centroids;
            nrms = fluxgeo.Normals;
            ars = fluxgeo.Areas;
            surf = ControlSurface(pnts, nrms, ars);
            
            flux = wf.energyFlux(surf, varargin{:}, 'PerPoint');

            fluxgeo = cell(size(flux));
            N = numel(flux);
            
            for n = 1:N
                fluxgeon = PanelGeo(wf.bodyGeo);
                fluxgeon.Values = flux{n};
                fluxgeo{n} = fluxgeon;
            end
        end
    end
end