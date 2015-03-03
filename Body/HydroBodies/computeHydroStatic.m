function [C, C0, Km, cb] = computeHydroStatic(rho, panelGeo, zpos, modes)

panelGeo = PanelGeo(panelGeo);
panelGeo.Translate([0 0 zpos]);

g = IWaves.G;

dof = modes.DoF;

motFuncs = modes.MotionFuncs;
for n = 1:length(motFuncs)
    motFuncs(n).Cg = motFuncs(n).Cg + [0 0 zpos];
    if (isa(motFuncs(n), 'HingeYFunc'))
        motFuncs(n).HingePos = motFuncs(n).HingePos + [0 zpos];
    end
end

cents = panelGeo.Centroids;
norms = panelGeo.Normals;
areas = panelGeo.Areas;
nPan = panelGeo.Count;
isWets = panelGeo.IsWets;
isIns = panelGeo.IsInteriors;
isBods = panelGeo.IsBodies;

% Compute the wetted volume, center of buoyancy
VolWet = [0 0 0];
VolBod = [0 0 0];
cb = [0 0 0];
cg = [0 0 0];

for n = 1:nPan
    
    pnt = cents(n,:);
    nrm = norms(n,:);
    if (isWets(n))
        VolWet = VolWet - nrm.*pnt*areas(n);
        cb = cb - nrm.*pnt.^2*areas(n);
    end
    
    if (isBods(n))
        VolBod = VolBod - nrm.*pnt*areas(n);
        cg = cg - nrm.*pnt.^2*areas(n);
    end
end

VolWet = VolWet(1);
VolBod = VolBod(1);

cb = cb./(2*VolWet);
cg = cg./(2*VolBod);

% Compute the mass 'static matrix'
rhoBod = rho*VolWet/VolBod;

Km = zeros(dof, dof);

for l = 1:dof
    for m = 1:dof
        kmlm = 0;
        
        for n = 1:nPan
            if (isBods(n))
                pnt = cents(n,:);        
                nrm = norms(n,:);
                
                fgl = motFuncs(l).GravityForce(pnt);
                
                nm = motFuncs(m).Evaluate(pnt);
                nmS = nm(1)*nrm(1) + nm(2)*nrm(2) + nm(3)*nrm(3);
                
                kmlm = kmlm + nmS*fgl*areas(n);
            end
        end
        Km(l, m) = kmlm;
    end
end

% znorm = [0 0 -1];
% 
% for m = 1:6
%     Fm = [0 0 0];
%     Mm = [0 0 0];
%     for n = 1:nPan
%         if (isBods(n))
%             pnt = cents(n,:);        
%             nrm = norms(n,:);
%             nm = motFuncs(m).Evaluate(pnt);
% 
%             nmS = nm(1)*nrm(1) + nm(2)*nrm(2) + nm(3)*nrm(3);
% 
%             Fm = Fm + nmS*znorm*areas(n);
%             Mm = Mm + nmS*cross(pnt-cg, znorm)*areas(n);
%         end
%     end
%     Km(1:3,m) = Fm';
%     Km(4:6,m) = Mm';
% end

Km = rhoBod*g*Km;
mKm = max(max(abs(Km)));
izero = (abs(Km) < 1e-8*mKm);
Km(izero) = 0;

% compute the hydrostatic matrix, not including the mass
C0 = zeros(dof, dof);

for l = 1:dof
    for m = 1:dof
        c0lm = 0;
        
        for n = 1:nPan
            if (isWets(n) && ~isIns(n))
                pnt = cents(n,:);
                nrm = norms(n,:);

                nl = motFuncs(l).Evaluate(pnt);
                divl = motFuncs(l).Divergence(pnt);

                nm = motFuncs(m).Evaluate(pnt);
                nmS = nm(1)*nrm(1) + nm(2)*nrm(2) + nm(3)*nrm(3);

                c0lm = c0lm + nmS*(nl(3) + pnt(3)*divl)*areas(n);
            end
        end
        C0(l, m) = c0lm;
    end    
end

C0 = -rho*g*C0;
mC0 = max(max(C0));
izero = (abs(C0) < 1e-8*mC0);
C0(izero) = 0;

% add the hydrostatic and mass static to get the final hydrostatic matrix
C = C0 - Km;

mC = max(max(C));
izero = (abs(C) < 1e-8*mC);
C(izero) = 0;


end