function erni = calculateErrorRelative2NoIntervention(eNormal2Cortex,tmap,E0,tMin)

w = abs(tmap);
if ~isempty(tMin)
tmap(abs(tmap)<=tMin) = 0;
w(tmap==0) = tMin;
end

y = E0 * tmap;
ew = w .* eNormal2Cortex;

erni = -2 * y .* ew + ew .* ew;

