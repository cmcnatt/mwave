classdef ResponseSpectrumUT < matlab.unittest.TestCase
    
    methods (Test)
        function test1(testCase)
            
            comp = floatingCylComp;
            Hs = 2;
            Tp = 7.5;
            spec = Bretschneider(Hs, Tp, comp.T);
            
            xi = comp.Motions;
            
            heave = ResponseSpectrum(squeeze(xi(:,1,3)), spec);
            
            figure;
            subplot(2,1,1);
            yyaxis right;
            plot(spec.Frequencies, spec.Spectrum);
            ylabel('Wave Spectrum (m^2/Hz)')
            
            yyaxis left;
            plot(1./comp.T, heave.Response);
            ylabel({'Heave Response' 'Spectrum (m^2/Hz)'});
            xlabel('Frequency (Hz)');
            
            title('ResponseSpectrumUT');
            
            subplot(2,1,2);
            yyaxis right;
            plot(1./spec.Frequencies, spec.Spectrum('Period'));
            ylabel('Wave Spectrum (m^2/s)')
            
            yyaxis left;
            plot(comp.T, heave.Response('Period'));
            ylabel({'Heave Response' 'Spectrum (m^2/s)'});
            xlabel('Period (s)');
            
            heaveLim = 1.0;
            prob = heave.ProbabilityOfExceeding(heaveLim);
            
            fprintf('\nRMS Response: %4.2f m\n', heave.RMSResponse);
            fprintf('Probability of exceeding %3.1f m motion: %4.2f\n\n', heaveLim, prob); 
            
        end
    end
end