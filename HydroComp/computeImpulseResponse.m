function [K, t] = computeImpulseResponse(B, omega, Nt, tmax, varargin)

Nomega = [];
omegaMin = [];
omegaMax = [];
    
if ~isempty(varargin)
    
    Nomega = varargin{1};

    if length(varargin) == 3 
        omegaMin = varargin{2};
        omegaMax = varargin{3};
    elseif length(varargin) == 2
        error('Three optional inputs are required: Nomega, omegaMin, omegaMax');
    end
end

waitB = waitbar(0,'Calculating radiation IRFs...');  % Progress bar

if size(B,1) ~= length(omega)
    error('The size of the damping matrix does not match the number of frequencies');
end

% Set defaults if empty
if isempty(tmax)
    tmax = 100;          
end

if isempty(Nt)  
    Nt = 1001;            
end

omegaB = omega;

interpB = false;
if ~isempty(Nomega)
    interpB = true;
    if isempty(omegaMin)
        omegaMin = min(omega);  
    end
    if isempty(omegaMax)
        omegaMax = max(omega);  
    end

    omega = linspace(omegaMin, omegaMax, Nomega);
end

t = linspace(0, tmax, Nt);
dof = size(B, 3);
K = zeros(Nt, dof, dof);

NN = length(t)*dof^2;

% Calculate the impulse response function for radiation
nn = 0;

for m = 1:dof;
    for n = 1:dof;
        ra_B = squeeze(B(:,m,n)).';
        if interpB
            ra_B = interp1(omegaB, ra_B, omega, 'spline');
        end
        for o = 1:length(t);
            %K(o,m,n) = (2/pi)*trapz(omega, ra_B.*(cos(omega*t(o)).*omega));
            K(o,m,n) = (2/pi)*trapz(omega, ra_B.*(cos(omega*t(o))));
            nn = nn+1;
        end
    end
    waitbar(nn/NN)
end


close(waitB)

end