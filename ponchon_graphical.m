clear;
close;
clc;

xe = [0 0.08 0.18 0.25 0.49 0.65 0.79 0.91 1];
ye = [0 0.28 0.43 0.51 0.73 0.83 0.9 0.96 1];
Hl = [24.3 24.1 23.2 22.8 22.05 21.75 21.7 21.6 21.4];
Hv = [61.2 59.6 58.5 58.1 56.5 55.2 54.4 53.8 53.3];

E_curve_y = polyfit(xe,ye,8);
E_curve_x = polyfit(ye,xe,8);

xe_new = 0:0.01:1;
for k = 1:length(xe_new) 
    ye_new(k) = polyval(E_curve_y,xe_new(k));
end

F = 1000/103.48;
zf = 0.42;
xf = 0.42;
xd = 0.97;
xw = 0.01;
R = 2.5;
D = 4.2681;
W = 5.7256;

figure(2)
plot(xe,Hl,'b');
hold on
plot(ye,Hv,'r');
ylim([-70 150]);
title('Enthalpy Concentration Diagram');
ylabel('H 10^-3 KJ/Kmol');
xlabel('x-y');


Hl_p = polyfit(xe,Hl,1);
Hv_p = polyfit(ye,Hv,1);
Hw = polyval(Hl_p,xw);
Hv2 = polyval(Hv_p,xw);
Hd = polyval(Hl_p,xd);
Hv1 = polyval(Hv_p,xd);
Hf = polyval(Hl_p,xf);

qd_dash = R*(Hv1-Hd) + Hv1;
slope_F = (qd_dash - Hf)/(xd-xf);
intercept_F = -slope_F*xd + qd_dash;
qw_dash = slope_F*xw + intercept_F;

plot(xd,qd_dash,'or','MarkerSize',5,'MarkerFaceColor','red');
plot(xw,qw_dash,'or','MarkerSize',5,'MarkerFaceColor','red');
plot(xf,Hf,'sr','MarkerSize',7,'MarkerFaceColor','red');
plot(xd,Hd,'^r','MarkerSize',6,'MarkerFaceColor','red');
plot(xw,Hw,'^r','MarkerSize',6,'MarkerFaceColor','red');
set(line([xd,xw],[qd_dash,qw_dash]),'Color',[1,0,1]);
set(line([xd,xd],[qd_dash,Hd]),'Color',[0,1,0]);
set(line([xw,xw],[qw_dash,Hw]),'Color',[0,1,0]);

xc(1) = xd;
yc(1) = xd;

x_int = xd;
y_int = xd;

i = 2;

while(x_int > xf)
     xc(i) = interp1(ye_new,xe_new,y_int);
   if(xc(i) >= xf)
       Hv_old = polyval(Hv_p,y_int);
       Hl_new = polyval(Hl_p,xc(i));
   
       plot([xc(i),y_int],[Hl_new, Hv_old],'--b');
       plot([xd,xc(i)],[qd_dash,Hl_new],'g');
   
       line = polyfit([xd,xc(i)],[qd_dash,Hl_new],1);
       y_line = polyval(line,xe);
   
       yc(i) = interp1(y_line - Hv,xe,0);
       
       y_int = yc(i);
       x_int = xc(i);
       
       i = i+1;
       
   else
       Hv_old = polyval(Hv_p,y_int);
%        plot([xf,y_int],[Hf, Hv_old],'--b');
       x_int = xc(i); 
   end
end

xstr(1) = xw;
ystr(1) = xw;

x_int = xw;
y_int = xw;


j = 2;
yf = 0;
while(x_int < xf )
   ystr(j) = polyval(E_curve_y,x_int);
   Hl_old = polyval(Hl_p,x_int);
   Hv_new = polyval(Hv_p,ystr(j));
   
   line = polyfit([xw,ystr(j)],[qw_dash,Hv_new],1);
   x_line = polyval(line,xe);
   
   xstr(j) = interp1(x_line - Hl,xe,0);
   if(xstr(j) < xf) 
       x_data = [ystr(j),x_int];
       y_data = [Hv_new, Hl_old];
       
       plot([ystr(j),x_int],[Hv_new, Hl_old],'--b');
       plot([xw,ystr(j)],[qw_dash,Hv_new],'g')
       
       x_int = xstr(j);
       y_int = ystr(j);
       
       j = j+1;
   else
       
       line = polyfit([xd,xf],[qd_dash,Hf],1);
       y_line = polyval(line,xe);
   
       yf = interp1(y_line - Hv,xe,0);
       Hfv = polyval(Hv_p,yf);
       
%        plot([yf,x_int],[Hfv, Hl_old],'--b');
       x_int = xstr(j);
       
       y_int = ystr(j);
   end
   
end

n1 = 1;
n2 = j-1;
while (n1 <= i + j - 1)
    
    if(n1 < i)
        x_points(n1) = xc(n1);
        y_points(n1) = yc(n1);
    elseif (n1 == i)
        x_points(n1) = xf;
        y_points(n1) = yf;
    else
        x_points(n1) = xstr(n2);
        y_points(n1) = ystr(n2);
        
        n2 = n2 - 1;
    end
    n1 = n1 + 1;
end

figure(1)
plot(xe,xe,'g');
hold on
plot(xe_new,ye_new,'r');
xlim([0 1]);
ylim([0 1]);

plot(x_points,y_points,'b');
plot(x_points,y_points,'or');
title('Equilibrium curve with stage construction');
xlabel('x (molefraction of n-heptane in liquid phase)');
ylabel('y (molefraction of n-heptane in vapour phase)');

xstage(1) = xd;
ystage(1) = xd;
k =1;

while (xstage(k) > xw)
    ydummy = ystage(k);
    xstage(k+1) = interp1(ye_new,xe_new,ydummy);
    ystage(k+1) = interp1(x_points,y_points,xstage(k+1));
    plot([xstage(k),xstage(k+1)],[ystage(k),ystage(k)],'g');
    if (xstage(k+1) > xw) 
        plot([xstage(k+1),xstage(k+1)],[ystage(k),ystage(k+1)],'g');
    end
    k = k+1;
end

hold off

n = 1;
for i=2:length(xstage)
    Xn(n) = xstage(i);
    Yn(n) = ystage(i);
    n= n+1;
end

Xn = Xn';
Yn = Yn';


for i=1:length(Xn)
    if(i <= 7) 
        x = Xn(i);
%         y = Yn(i);
        y = interp1(x_points,y_points,x);
        check(i) = y;
        V(i) = (D*xd - D*x)/(y-x);
        L(i) = V(i) - D;
    else
        x = Xn(i);
%         y = Yn(i);
        y = interp1(x_points,y_points,x);
        check(i) = y;
        L(i) = (W*xw - W*y)/(x-y);
        V(i) = L(i) - W;       
    end
end

HLn = interp1(xe,Hl,Xn);
HVn = interp1(xe,Hv,Yn);
Ln = L';
Vn = V';

figure(3)
x_names = {'x', 'y','L','V','HL','HV'};
y_names = 1:1:11;
h = heatmap(x_names,y_names,[Xn,Yn,Ln,Vn,HLn,HVn]);

h.Title = 'By using ponchon savarit method'; 
h.YLabel = 'x,y,L,HL,HV in each stage';
h.XLabel = 'All variables in each stage';

