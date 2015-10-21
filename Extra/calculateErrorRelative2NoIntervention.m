function erni = calculateErrorRelative2NoIntervention(eNormal2Cortex,tmap,E0,tMin)


tmap_thresholded = tmap;
tmap_thresholded(abs(tmap_thresholded)<=2) = 0;
w = abs(tmap_thresholded);
w(w==0) = 2;

y = E0 * tmap_thresholded;
ew = w .* eNormal2Cortex;

erni = -2 * y .* ew + ew .* ew;

