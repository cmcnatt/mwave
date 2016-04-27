classdef MColor
    properties (Constant)
        Blue = [0 0.4470 0.7410];
        Red = [0.8500 0.3250 0.0980];
        Yellow = [0.9290 0.6940 0.1250];
        Purple = [0.4940 0.1840 0.5560];
        Green = [0.4660 0.6740 0.1880];
        LightBlue = [0.3010 0.7450 0.9330];
        Maroon = [0.6350 0.0780 0.1840];
        Black = [0 0 0];
        White = [1 1 1];
    end
    
    methods (Static)
        function [cols] = Colors()
            cols = cell(7,1);
            cols{1} = MColor.Blue;
            cols{2} = MColor.Red;
            cols{3} = MColor.Yellow;
            cols{4} = MColor.Purple;
            cols{5} = MColor.Green;
            cols{6} = MColor.LightBlue;
            cols{7} = MColor.Maroon;
        end
    end
end