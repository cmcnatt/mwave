function [vol, cg] = computeVolume(panelGeo)

cents = panelGeo.Centroids;
norms = panelGeo.Normals;
areas = panelGeo.Areas;
nPan = panelGeo.Count;
isBods = panelGeo.IsBodies;

% Compute the wetted volume, center of buoyancy
vol = [0 0 0];
cg = [0 0 0];

for n = 1:nPan
    
    pnt = cents(n,:);
    nrm = norms(n,:);
    area = areas(n);

    if (isBods(n))
        vol = vol - nrm.*pnt*area;
        cg = cg - nrm.*pnt.^2*area;
    end
end

vol = vol(2);
cg = cg./(2*vol);
vol = -vol;
