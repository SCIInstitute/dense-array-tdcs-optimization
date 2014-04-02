function x = lcs(d,i,j)

ind = evalin('base','ind');
tot = evalin('base','tot');
p = evalin('base','p');
Te = evalin('base','Te');

indi = ind(i);
pi = p(j);

aa = load(['roi1/elecCurrent' num2str(d) '1' num2str(i) num2str(j) '.mat']);

nun = findNumUnknown(aa.currentArray);

if (nun > 5)
    disp('number of unknowns greater than 7, exiting....');
else
    disp(['numbers of unknowns: ' num2str(nun)]);
    wq = load(['wQMatrix' num2str(d) '.mat']);
    w = wq.w;
    Q = wq.Q;
    x = oneSourceOptimization(w,Q,tot,indi,pi,Te);
end

end

function nun = findNumUnknown(x)

minX = min(x);
maxX = max(x);

unknownNegatives = find(x > 0.95*minX & x < -0.001* min(maxX,-minX));
unknownPositives = find(x < 0.95*maxX & x > 0.001* min(maxX,-minX));

nN = numel(unknownNegatives);
nP = numel(unknownPositives);

nun = nP+ nN;

end





