function [dmap] = mdensity4(tx, ty, dpar, w, E_x, E_y, mapx, mapy)
% Written by Yu-Jen Chen, Department of Anatomy and Structural Biology,
% Albert Einstein College of Medince, Coleman laboratory

% The function 'mdensity' calculate heat map for the given parameter 'dpar'.
% 'dpar' is associated with the coordinates [x, y]. If only density (or the 
% appearance of a 'signal' at [x, y]) is wanted, put the string 'nopar' for
% dpar. 
% The heat map is weighted and only calculated signals within the boundary
% defined by [xbound, ybound].
% Output dmap has dimension [mapy, mapx].
%
% Updates:
% 7/7/2015  V4.0
% Use circular window for density calculation.
%
% 7/03/2015 V3.0
% use a matrix (Mdata) for counting tracks in the window
% create a cicular boundary to calculate the density
%
% 7/02/2015 V2.0   
% Note: use original 25 mem defining points for xbound and ybound
%
% 5/14/2015 v1.0
% Note: When dpar is provided, the heat map is not normalized to the number
% of weighted pixels counted (the surface area of the window). 
% 
% 

Mx = round(tx); My = round(ty);
tf = strcmp(dpar,'nopar');

ths = 0:1:8;
th = (pi/4)*ths + pi/8;

dmap = zeros(mapy, mapx); 

[NE_x, NE_y] = membound(E_x, E_y, 1);

avweight = poly2mask(NE_x, NE_y, mapy, mapx);

if tf == 1
    Mdata = zeros(mapy, mapx);
    
    for k = 1:length(Mx)
        Mdata(My(k), Mx(k)) = Mdata(My(k), Mx(k)) + 1;
    end
    
    x0 = w + 1; y0 = w + 1;
    R = w / cos(pi/16);
    xi = R * cos(th) + x0; yi = R * sin(th) + y0;

    inW = poly2mask(xi, yi, 2*w+1, 2*w+1);
    
    [in_r, in_c] = find(avweight == 1);

    for j = 1:length(in_r)
        
        miniW = avweight(in_r(j)-w:in_r(j)+w, in_c(j)-w:in_c(j)+w);
        miniM = Mdata(in_r(j)-w:in_r(j)+w, in_c(j)-w:in_c(j)+w);
        
        pxnum = sum(sum(inW .* miniW));
        in_num = sum(sum(miniM .* inW));

        dmap(in_r(j), in_c(j)) = in_num / pxnum;
        
    end
    
else
    Mdata = zeros(mapy, mapx); Pdata = zeros(mapy, mapx);
    
    for k = 1:length(Mx)
        Mdata(My(k), Mx(k)) = Mdata(My(k), Mx(k)) + 1;
        Pdata(My(k), Mx(k)) = Pdata(My(k), Mx(k)) + dpar(k);
    end
    
    x0 = w + 1; y0 = w + 1;
    R = w / cos(pi/16);
    xi = R * cos(th) + x0; yi = R * sin(th) + y0;

    inW = poly2mask(xi, yi, 2*w+1, 2*w+1);
    
    [in_r, in_c] = find(avweight == 1);

    for j = 1:length(in_r)
        
        miniW = avweight(in_r(j)-w:in_r(j)+w, in_c(j)-w:in_c(j)+w);
        miniM = Mdata(in_r(j)-w:in_r(j)+w, in_c(j)-w:in_c(j)+w);
        miniP = Pdata(in_r(j)-w:in_r(j)+w, in_c(j)-w:in_c(j)+w);
        
        %pxnum = sum(sum(inW .* miniW));
        in_par = sum(sum(miniP .* inW .* miniW));
        in_num = sum(sum(miniM .* inW .* miniW));

        dmap(in_r(j),in_c(j)) = in_par / in_num;
        
    end
  
end

end