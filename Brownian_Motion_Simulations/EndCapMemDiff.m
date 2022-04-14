function [xPos2,xcap2,yPos2,zPos2] = EndCapMemDiff(r,l,phi,d,zPrevPos,xPrevPos,yPrevPos)
%Written by Alecia Achimovich
%Define the radius of the circle of the circle-circle intersection (of the
%sphere defining the cell cap & the sphere defining the movement vector).

%Use the small angle approximating to find the angle
theta_sa = d/r;
a = r*sin(theta_sa);
%a = (1/(2*r))*sqrt(4*r^4-(2*r^2-d^2)^2);
xDelta = a*cos(phi);
yDelta = a*sin(phi);
%Define vector normal to starting position
if xPrevPos <0
    xcap = xPrevPos+l;
elseif xPrevPos >0
    xcap=xPrevPos-l;
end
u = [xcap; yPrevPos; zPrevPos];
%Define vector perpendicular to the normal vector. This vector will be
%parallel to the plane of the circle around starting position.
if abs(xcap) <= 0.05 && abs(yPrevPos) <= 0.025
    v = [(u(3));0;(-u(1))];
else
    v = [-(u(2));u(1);0];
end
%Define vector perpendicular to both u & v by finding the cross product.
w = cross(u,v);
%Normalize each vector component (so that the magnitude is equal to one)
%such that it defines a unit vector of the circle.
mag_u = sqrt(u(1)^2+u(2)^2+u(3)^2);
norm_u = u./mag_u;

mag_v = sqrt(v(1)^2+v(2)^2+v(3)^2);
norm_v = v./mag_v;

mag_w = sqrt(w(1)^2+w(2)^2+w(3)^2);
norm_w = w./mag_w;

%Define vector orginating from previous position
p = [(xDelta.*norm_v)+(yDelta.*norm_w)];
%Find origin of intersecting circle
%Find where these imaginary numbers are coming from.

test = d^2-(sqrt(p(1)^2+p(2)^2+p(3)^2))^2;
psi = acos(a/d);
if test > 10^(-29)
    delta = d*sin(psi);
    scaleu= (mag_u-delta)/mag_u;
    oCircle = scaleu.*u;
    %Find final vector
    newPos = oCircle+p;
else
    newPos = u+p;
end

yPos = newPos(2);
zPos = newPos(3);
xcappre = newPos(1);
%Make sure that new position is face of cylinder
scalor = r/(sqrt(yPos^2+zPos^2+xcappre^2));
newPosScaled = scalor.*newPos;

yPos2 = newPosScaled(2);
zPos2 = newPosScaled(3);
xcap2 = newPosScaled(1);
% yPos2 = newPos(2,1);
% zPos2 = newPos(3,1);
% xcap2 = newPos(1,1);
%Redefine x in terms of spherocylinder coordinates
% check = isreal(zPos2);
% check = check;
% errPos = [0,0,0,0,0,0,0];
if xPrevPos < 0 && xcap2 > 0 || xPrevPos > 0 && xcap2 < 0 
%     xcap2 = -xcap2;
%     if xcap2 > 0
%         xPos2 = xcap2+l;
%     elseif xcap2 < 0
%         xPos2 = xcap2 - l;
%     end
    dxover = abs(xcap2);
    sphereD = 1-(dxover/(abs(xcap2-xcap)));
    dy = yPos2-yPrevPos;
    dySphere = dy*sphereD;
    ySpherePos = yPrevPos + dySphere;
    dz = zPos2-zPrevPos;
    dzSphere = dz*sphereD;
    zSpherePos = zPrevPos + dzSphere;
    dyzCyl = sqrt(((1-sphereD)*dy)^2+((1-sphereD)*dz)^2);
        if dz < 0
            dyzCyl = -dyzCyl;
        end
    if xcap >= 0 && xcap2 < 0
        dxover=-dxover;
        xSpherePos = l;
    elseif xcap <= 0 && xcap2 > 0
        xSpherePos=-l;
    end
[xPos2, yPos2, zPos2] = CylMemDiff(r,l,dxover,dyzCyl,xSpherePos,ySpherePos,zSpherePos)  
% check = isreal(zPos2);
% check = check;
% if abs(xPos2) > l
%  errPos = [xPrevPos,xPos2,yPrevPos,yPos2,zPrevPos,zPos2,phi];
% end
else
    if xcap2 > 0
        xPos2 = xcap2+l;
    elseif xcap2 < 0
        xPos2 = xcap2 - l;
    end
    
end
end