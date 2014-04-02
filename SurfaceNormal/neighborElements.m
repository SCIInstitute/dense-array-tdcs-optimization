function neighborList = neighborElements(elem,numOfCommNodes4Neigh)
switch numOfCommNodes4Neigh
    case 1
        disp('1')
    case 2
        disp('2')
    case 3
        disp('3')
    case 4
        disp('4')
    otherwise
        disp('No neighborhood relation for this definition')
end
neighborList = [];
end