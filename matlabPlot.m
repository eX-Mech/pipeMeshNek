%close all
clear
clc

r = load('plotElements');

R  = 0.5;
nR = 5;

nVe = size(r,1);
nEl = nVe/8;

% plot front view
figure(1)
hold off
% outer circle
%
x = linspace(-R,R,100);
plot(x, (R^2 - x.^2).^0.5,'-k','linewidth',2)
hold on
plot(x,-(R^2 - x.^2).^0.5,'-k','linewidth',2)
% inner round circles
%
%for i = 1 : nR
%	x = linspace(-R/nR*i,R/nR*i,100);
%	plot(x, ((R/nR*i)^2 - x.^2).^0.5,'-k','linewidth',1)
%	plot(x,-((R/nR*i)^2 - x.^2).^0.5,'-k','linewidth',1)
%end
% mesh
%
for i = 1:8:nVe

   plot([r(i,  2) r(i+1,2)], [r(i,  1) r(i+1,1)], '-b')
   plot([r(i+1,2) r(i+2,2)], [r(i+1,1) r(i+2,1)], '-b')
   plot([r(i+2,2) r(i+3,2)], [r(i+2,1) r(i+3,1)], '-b')
   plot([r(i+3,2) r(i,  2)], [r(i+3,1) r(i,  1)], '-b')

   plot(r(i,  2), r(i,  1), '.k')
   plot(r(i+1,2), r(i+1,1), '.k')
   plot(r(i+2,2), r(i+2,1), '.k')
   plot(r(i+3,2), r(i+3,1), '.k')

end
grid on
axis equal
%axis([-0.05 .55 -0.05 .55])
xlabel('y')
ylabel('x')
