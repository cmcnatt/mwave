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
classdef VideoTimeTrans
    properties (Access = protected)
        tvid;
        texp;
    end
    
    methods
         
        function [vtt] = VideoTimeTrans(vidObj, expdt)
            Nfr = vidObj.Duration*vidObj.FrameRate;
            vtt.tvid = 0:1/vidObj.FrameRate:vidObj.Duration;
            vtt.texp = 0:expdt:Nfr*expdt;
        end
        
        function [val] = VidTime(vtt, timeIn)
            iframe = indexOf(vtt.texp, timeIn);
            val = vtt.tvid(iframe);
        end
     end
end