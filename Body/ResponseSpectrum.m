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
classdef ResponseSpectrum < handle
    
    properties (Access = private)
        rao;
        spec;
        nf;
        nd;
    end
    
    properties (Dependent)
        WaveSpectrum;
        RAO;
        RMSResponse;
    end
    
    methods
        
        function [rspec] = ResponseSpectrum(varargin)
            if ~isempty(varargin)
                rspec.setCheckSize(varargin{1});
                rspec.rao = varargin{1};
                if length(varargin) > 1
                    if ~isa(varargin{2}, 'WaveSpectrum')
                        error('Second input must be of type WaveSpectrum');
                    end
                    rspec.setCheckSize(varargin{2}.Spectrum);
                    rspec.spec = varargin{2};
                end
            end
        end
        
        function [val] = get.WaveSpectrum(rspec)
            % the wave spectrum used in the response
            val = rspec.rao;
        end
        function [] = set.WaveSpectrum(rspec, val)
            % the wave spectrum used in the response
            if ~isa(val, 'WaveSpectrum')
                error('WaveSpectrum must be of type WaveSpectrum');
            end

            rspec.setCheckSize(val.Spectrum);
            rspec.spec = val;
        end
        
        function [val] = get.RAO(rspec)
            % The response amplitude operator
            val = rspec.rao;
        end
        function [] = set.RAO(rspec, val)
            % The response amplitude operator
            rspec.setCheckSize(val);
            rspec.rao = val;
        end
        
        function [val] = Response(rspec, varargin)
            % Get the response spectrum
            val = [];
            if ~isempty(rspec.spec) && ~isempty(rspec.rao)
                val = abs(rspec.rao).^2.*rspec.spec.Spectrum(varargin{:});
            end
        end
        
        function [val] = get.RMSResponse(rspec)
            % Get the RMS response spectrum
            val = [];
            if ~isempty(rspec.spec) && ~isempty(rspec.rao)
                res = abs(rspec.rao).^2.*rspec.spec.Spectrum;
                
                [df, ddir] = rspec.spec.Deltas('Both');
                [Ddir, Df] = meshgrid(ddir, df);
                
                val = sqrt(abs(sum(sum(res.*Ddir.*Df))));
            end
        end
        
        function [p] = ProbabilityOfExceeding(rspec, val)
            % compute the probability of exceeding a particular value
            sigma = rspec.RMSResponse;
            p = [];
            if ~isempty(sigma)
                p = exp(-val^2/(2*sigma^2));
            end
        end
    end
    
    methods (Access = private)
        function [] = setCheckSize(rspec, res)

            if isvector(res)
                Nf = length(res);
                Nd = 1;
            else
                [Nf, Nd] = size(res);
            end

            if isempty(rspec.nf)
                rspec.nf = Nf;
                rspec.nd = Nd;
            else
                if (rspec.nf ~= Nf) || (rspec.nd ~= Nd)
                    error('The Spectrum or RAO is not the correct size');
                end
            end
        end
    end
end