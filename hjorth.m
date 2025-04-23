function [activity,mobility,complexity] = hjorth(sig)
    dx  = diff(sig);
    ddx = diff(dx);
    activity = var(sig);
    mobility = sqrt(var(dx)/activity);
    complexity = sqrt(var(ddx)/var(dx)) / mobility;
end
