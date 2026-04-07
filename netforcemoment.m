function output = netforcemoment(guess,buoyantdistfunction,selfWeightInterp,xInterp,xUnique,paddlers)
% NETFORCEMOMENT Determines the net force and moment from a trim/draft of a
% buoyantdist func.
% 
%   Detailed explanation goes here
    buoyantdist = buoyantdistfunction(guess);
    buoyantdistInterp = interp1(xUnique,buoyantdist,xInterp,"linear");
    netForceDist = selfWeightInterp + buoyantdistInterp;
    distributedForce = trapz(xInterp,netForceDist);
    Fnet = distributedForce + sum(paddlers(2,:));
    Mnet = momentabout(0,netForceDist,xInterp) + sum(paddlers(1,:) .* paddlers(2,:));
    output = [Fnet, Mnet];
end

