function centers = elemCenter(node,elem)
%%Finds the centers of the elements by taking average of the
%nodes belonging to same element.
%Written by Seyhmus Guler, 3/8/14
%INPUTS:
%node: Node positions. Size: #dimension x #nodes
%elem: Elements. Size: #nodes per element x # elements
%OUTPUTS:
%centers: centers of the elements. Size: #dimension x #elements

[nD,~] = size(node); %node dimension (nD) and number of nodes (~)
[eD,M] = size(elem); %element dimension (eD) and number of elements (M)

expandedNode = zeros(eD,nD*M);
for i=1:eD
    expandedNode(i,:) = reshape(node(:,elem(i,:)),1,[]);
end

centers = reshape(sum(expandedNode)/eD,nD,[]);