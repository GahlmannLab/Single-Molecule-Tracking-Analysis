function [xPos2, yPos2, zPos2,errPos] = CylMemDiff(r,l,dx,dyz,xPrevPos,yPrevPos,zPrevPos)
%Written by Alecia Achimovich
%Use small angle theorem to find angle corresponding to z-y diffusion
phi = dyz/r;
%Find elevation angle of previous position
rota=[1 0 0 ;0 cos(phi) sin(phi);0 -sin(phi) cos(phi)];
rotap=rota*[xPrevPos;yPrevPos;zPrevPos];
x2=rotap(1);
yPos2=rotap(2);
zPos2=rotap(3);
xPos2 = x2+dx;
crossSph = 0;
errPos = [0,0,0,0,0,0];
if abs(xPos2) > l
 dover = abs(xPos2)-l;  
 if xPos2 < -l
     dover = -dover;
 end
    magV = sqrt(dover^2 + yPos2^2 + zPos2^2);
    Pos = [dover;yPos2;zPos2];
    scalor = 0.4/magV;
    newPos = scalor.*Pos;
    xcap2 = newPos(1);
    yPos2 = newPos(2);
    zPos2 = newPos(3);
    if xPos2 <= -l
        xPos2 = -l+xcap2
    elseif xPos2 >= l
        xPos2 = l+xcap2;
    end
        
%     if dx < 0
%         dover = dover*(-1);
%     end
%     dcyl = (abs(dx)-abs(dover))/abs(dx);
%     dyzcyl = dyz*dcyl;
%     phicyl = dyzcyl/r;
%     rota=[1 0 0 ;0 cos(phi) sin(phi);0 -sin(phi) cos(phi)];
%     rotap=rota*[xPrevPos;yPrevPos;zPrevPos];
%     yPosCyl = rotap(2);
%     zPosCyl = rotap(3);
%     if xPos2 > 0
%         xPosCyl = l+(10^(-30));
%     elseif xPos2 < 0
%         xPosCyl = -l-(10^(-30));
%     end
%     dyzsphere =dyz*(1-dcyl);
%     dsphereTotal = sqrt(dyzsphere^2 + dover^2);    
%     
%         if yPrevPos > 0  && xPos2 > 0  || yPrevPos < 0  && xPos2 < 0
%             dover = -dover;      
%            
%         end        
%     thetaSphere = atan(dyzsphere/dover);
%         if thetaSphere > 0
%             if dover > 0
%                 thetaSphere = thetaSphere;
%             elseif dover < 0 
%                 thetaSphere = thetaSphere + pi;
%             end
%         elseif thetaSphere < 0
%             if dover < 0
%                 thetaSphere = thetaSphere + pi;
%             elseif dover > 0
%                 thetaSphere = thetaSphere + 2*pi;
%             end
%         end
%     [xPos2,xcap2,yPos2,zPos2] = EndCapMemDiff(r,l,thetaSphere,dsphereTotal,zPosCyl,xPosCyl,yPosCyl);
%     
    if abs(xPos2) < abs(l)        
        errPos = [xPrevPos,xPos2,yPrevPos,yPos2,zPrevPos,zPos2];
%     
    end
end
end
