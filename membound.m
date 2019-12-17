function [new_x, new_y] = membound(mem_x, mem_y, enfac) 
%   
% Updates:
%   07/07/2015 V. 2.0
%   Create continuous points pixel by pixel without extra redundant points.
%   Magnify the map by the factor set with enfac.
%
%   04/20/2015 v. 1.0
%   membrane connects manually defined spots with additional points.  
%   The mem_x and mem_y are the coordinates of the difined spots.
%   enfac is the factor for enlarging the matrix size, which improves the
%       accruacy when difining the boundary.
%   intfac is the number or queries for interpolating between the defined spots.
%   The function exports new coordinates of the linked membrane [mew_x, new_y]. 
%
%   * tested values for enfac: 10
%   * tested values for intfac: 100 for small region 
%                               1000 for the nuclear envelop
%
memc_x = round([mem_x; mem_x(1)] * enfac);
memc_y = round([mem_y; mem_y(1)] * enfac);

new_x = []; new_y = [];

for i = 1: length(memc_x)-1
    no_pt = max([abs(memc_x(i+1) - memc_x(i)) + 1, abs(memc_y(i+1) - memc_y(i)) + 1]);
    
    xline = linspace(memc_x(i), memc_x(i+1), no_pt);
    xint = xline(2:no_pt-1);
    
    yline = linspace(memc_y(i), memc_y(i+1), no_pt);
    yint = yline(2:no_pt-1);  
    
    new_x = cat(1, new_x, memc_x(i), xint(:));
    new_y = cat(1, new_y, memc_y(i), yint(:));
end

end





    
    
    





