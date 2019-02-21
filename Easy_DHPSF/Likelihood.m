function y = Likelihood(par,Zdata,bkgrnd,datavar, ii,jj)
%PDF Summary of this function goes here
%   Detailed explanation goes here
%f =a(1)*exp(-((x-a(3)).^2+(x-a(4)).^2)/(2*a(7)^2))+a(2)*exp(-((x-a(5)).^2+(x-a(6)).^2)/(2*a(8).^2));
Zmodel = (par(1).*exp( -((ii-par(3)).^2+(jj-par(4)).^2) / (2*par(7).^2)) ...
          +par(2).*exp( -((ii-par(5)).^2+(jj-par(6)).^2) / (2*par(8).^2)) ...
          +bkgrnd+par(9));
lowVal = min(min(min(Zmodel))-1e-6, 0);
par(9) = par(9) + abs(lowVal);
Zmodel= Zmodel+ datavar + abs(lowVal);
%  if abs(lowVal) > max(max(Zdata))*.05
%     warning('background fitting is in error'); 
%     par
%  end
%Zmodel= Zmodel+ datavar;
model = reshape(Zmodel, 1, []);
%%model = (reshape(Zmodel, 1, []))./10^4;
%Just for Bright beads, adjust because of bad fit
% for m =1:17
% for n =1:17
%     if (Zdata(m, n) < .2*max(max(Zdata)))
%     Zdata(m, n) = bkgrnd + datavar(m, n);
%     else
%     Zdata(m, n) = (Zdata(m,n));
%     end
% end
% end
data1 = reshape(Zdata,1, []);
%%data1 = (reshape(Zdata,1, []))./10^4;
y = sum(-((data1.*log(model)-model)));
%%y = sum(-(log((model).^(data1).*(exp(-model))./(gamma(data1+1)))));