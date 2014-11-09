function G=MakeEllipse(a,b,cx,cy)
% Makes a structure for each ellipse containing: 
% parametric curve: (rx ry)
% normal vector: (nx,ny) (not normalized)
% curvature:     (curv)
% length:        (norm of tangent/normal vector)

G.rx=@(t) cx+a*cos(t);  %parametric curve
G.ry=@(t) cy+b*sin(t);

rxp=@(t) -a*sin(t);     %derivative of curve
ryp=@(t) b*cos(t);

G.length=@(t) (norm([rxp(t) ryp(t)]));  % norm of derivative

sxp=@(t) -a*cos(t)*b^2./((b*cos(t)).^2+a^2-(a*cos(t)).^2).^(3/2);   % derivative of tangent
syp=@(t) -b*sin(t)*a^2./((b*cos(t)).^2+a^2-(a*cos(t)).^2).^(3/2);

G.nx=@(t) ryp(t);    % normal vector (not normalized)
G.ny=@(t) (-rxp(t));

G.curv=@(t) norm([sxp(t) syp(t)])./norm([rxp(t) ryp(t)]);  % curvature