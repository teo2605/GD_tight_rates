function cons=Hypoconvex(pt1,pt2,L,mu)
%
% This routine implements the interpolation conditions for non-convex
% smooth functions. Two parameter may be provided: 
%       - smoothness constant L (L is nonnegative and possibly infinite).
%       - hypoconvex constant mu, |mu|<=L.
%
% Note: not defining the value of L automatically corresponds to set L=Inf.
%
% To generate a non-convex smooth function 'f' with smoothness paramater
% L=1 from an instance of PEP called P:
%  >> P=pep();
%  >> param.L=1;
%  >> f=P.AddObjective('Hypoconvex',param);
%
% For details about interpolation conditions, we refer to the following
% references:
%
% (1) Taylor, Adrien B., Julien M. Hendrickx, and François Glineur. 
%     "Smooth strongly convex interpolation and exact worst-case 
%     performance of first-order methods." 
%     Mathematical Programming 161.1-2 (2017): 307-345.
%
% (2) Taylor, Adrien B., Julien M. Hendrickx, and François Glineur.
%     "Exact Worst-case Performance of First-order Methods for Composite
%     Convex Optimization." SIAM Journal on Optimization (2017)
%
% assert(L>=0 & L>=mu,'Constants provided to the functional class are not valid');
if ~(pt1.x.isEqual(pt2.x) && pt1.g.isEqual(pt2.g) && pt1.f.isEqual(pt2.f))
    if mu == L
       cons = ( (pt1.f - pt2.f + pt1.g * (pt2.x - pt1.x)) + L/2 * (pt2.x-pt1.x)^2  == 0);        
    else
        if L == 0
            cons=( (pt1.f - pt2.f + pt1.g * (pt2.x - pt1.x) + ...
                    -1/(2*mu) * (pt1.g-pt2.g)^2 + ...
                     (pt1.x-pt2.x)*(pt1.g-pt2.g)) <= 0);
        else
            if L ~= Inf % L < Inf, smooth        
                % [Theorem 3.3, 2] Smooth - hypoconvex
                cons=((pt1.f - pt2.f + pt1.g * (pt2.x - pt1.x) + ...
                        1/(2*(1-mu/L)) * ...
                        (1/L*(pt1.g-pt2.g)^2 + ...
                        mu*(pt1.x-pt2.x)^2 - ...
                        2*mu/L*(pt1.x-pt2.x)*(pt1.g-pt2.g))) <= 0);      
            else % non-smooth hypoconvex
                cons=((pt1.f - pt2.f + pt1.g * (pt2.x - pt1.x) + ...
                        mu/2 * (pt1.x - pt2.x)^2) <= 0);             
            end 
        end       
    end
else
    cons=[];
end
end
