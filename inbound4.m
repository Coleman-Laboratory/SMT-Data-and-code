function [in_x, in_y ,in_id] = inbound4(mol_x, mol_y, mol_id, bound_x, bound_y, enfac) 
% mol_x = nMol_x; mol_y = nMol_y; mol_id = nucl_i;
% % p = 1; 
% enfac = 10;
% % zbd_x = periZone{p}(:,1); zbd_y = periZone{p}(:,2);
% % [newzb_x, newzb_y] = membound(zbd_y, zbd_x, 1);
% bound_x = new_x; bound_y = new_y;
% function "inbound"
%   
% Updates:
%   08/12/2015 V. 4.0
%   Use inpolygon to test the track position one by one.
%
%   08/06/2015 V. 3.0
%   No continuouse boundary is required when using poly2mask.  
%
%   07/17/2015 v. 2.1
%   Use original bound_x and bound_y for input, magnification is calculated 
%   within the function. 
%   Use poly2mask to create the testing matrix.   
%
%   07/07/2015 v. 2.0
%   The input bound_x & bound_y are already magnified by the factor, enfac.
%   The update adjusted the coordinates from input, mol_x and mol_y by the
%   factor, enfac.
%
%   04/17/2015 v. 1.0
%   molx and moly are the coordinates for (x, y) of the objects for test.
%   mol_id are the id of the objects. 
%   bound_x and bound_y are the coordinates of the boundary where
%       objects are tested for whether they are inside it. 
%   The function returns the object coordinates in_x and in_y of those that
%       are inside the boundary. The indices of the objects are also
%       returned as in_id which refer to the input matrices. 
%   
%


test_i = zeros(length(mol_id),1);

cor_x = round(mol_x * enfac);
cor_y = round(mol_y * enfac);

% rbound = zeros(max(cor_y), max(cor_x));

[eb_x, eb_y] = membound(bound_y, bound_x, enfac);

% ebound_x = bound_x * enfac;
% ebound_y = bound_y * enfac;

% rbound_x = round(eb_x);
% rbound_y = round(eb_y);

% for i = 1:length(rbound_x)
%     rbound(rbound_y(i),rbound_x(i)) = 1;
% end
% 
% %filter = imfill(rbound,'holes');

% filter = poly2mask(ebound_y, ebound_x, max(cor_y), max(cor_x));

for i = 1:length(cor_x);
    inb = inpolygon(cor_x(i), cor_y(i), eb_y, eb_x);
    if inb == 1
%     if filter(cor_y(i), cor_x(i)) == 1
        test_i(i) = 1;
    else
        test_i(i) = 0;
    end
end

in_x = mol_x(logical(test_i));
in_y = mol_y(logical(test_i));
in_id = mol_id(logical(test_i));


end

    
    
    





