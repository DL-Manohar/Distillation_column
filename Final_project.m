clear;
close;
clc;

global X

% information given in the question
F = 1000/103.48; % Kmol/hr
zf = 0.42;  % feed conc.
xf = 0.42; % feed conc.
xd = 0.97; % distillate composition
xw = 0.01; % bottoms composition
R = 2.5; % reflux ratio
D = 4.2681; % distillate in kmol/hr, calculated from overall material balance
W = 5.7256; % bottoms in kmol/hr, calculated from overall material balance

% no. of stages under observation in trial and error method
n=11;

% feed tray at no. 6
feed_tray=6;

% equilibrium data
xe = [0 0.08 0.18 0.25 0.49 0.65 0.79 0.91 1];
ye = [0 0.28 0.43 0.51 0.73 0.83 0.9 0.96 1];

E_curve_y = polyfit(xe,ye,8);

xe_new = 0:0.001:1;
for k = 1:length(xe_new) 
    ye_new(k) = polyval(E_curve_y,xe_new(k));
end

% recifying section functions
fun1 = @Rectifying_equations;

% calculation of V1 and other varibles for first tray from given xd 
X = xd;
x0 = [0.5002,10.82,15.09,21.63,53.85];
y = fsolve(fun1,x0);

% variables to store the data at each tray
x_tray(1)=X;
y_tray(1)=y(1);
L_tray(1)=y(2);
V_tray(1)=y(3);
Hl_tray(1)=y(4);
Hv_tray(1)=y(5);

% for all the tray in rectifying section we use the 
% same strategy, calculate x from previous iteration 
% from equilibrium data, and then in the current iteration 
% use this x value to calculate variables namely
% y,L,V,Hl,Hv

for i= 1:1:feed_tray
    X = interp1(ye_new,xe_new,y(1));
    x_tray(i)=X;
    fun1 = @Rectifying_equations;
    x0 = [0.5002,10.82,15.09,21.63,53.85];
    y = fsolve(fun1,x0);
    y_tray(i)=y(1);
    L_tray(i)=y(2);
    V_tray(i)=y(3);
    Hl_tray(i)=y(4);
    Hv_tray(i)=y(5);
end

% for all the tray in stripping section we use the 
% same strategy yet again, calculate x from previous iteration 
% from equilibrium data, and then in the current iteration 
% use this x value to calculate variables namely
% y,L,V,Hl,Hv

for i= (feed_tray+1):1:n
    X = interp1(ye_new,xe_new,y(1));
    x_tray(i)=X;
    fun2 = @Stripping_equations;
    x0 = [0.5002,10.82,15.09,21.63,53.85];
    y = fsolve(fun2,x0);
    y_tray(i)=y(1);
    L_tray(i)=y(2);
    V_tray(i)=y(3);
    Hl_tray(i)=y(4);
    Hv_tray(i)=y(5);
end


x_tray = x_tray';
y_tray = y_tray';
L_tray = L_tray';
V_tray = V_tray';
Hl_tray = Hl_tray';
Hv_tray = Hv_tray';

figure(1)
x_names = {'x', 'y','L','V','HL','HV'};
y_names = 1:1:11;
h = heatmap(x_names,y_names,[x_tray,y_tray,L_tray,V_tray,Hl_tray,Hv_tray]);

h.Title = 'using fsolve to solve equations at every stage'; 
h.YLabel = 'x,y,L,HL,HV in each stage';
h.XLabel = 'All variables in each stage';




% rectifying section equations
function F = Rectifying_equations(x)

    global X
    
    % information given in the question
    zf = 0.42;
    xf = 0.42;
    xd = 0.97;
    xw = 0.01;
    R = 2.5;
    D = 4.2681;
    W = 5.7256;

    % data given in question
    xe = [0 0.08 0.18 0.25 0.49 0.65 0.79 0.91 1];
    ye = [0 0.28 0.43 0.51 0.73 0.83 0.9 0.96 1];
    Hl = [24.3 24.1 23.2 22.8 22.05 21.75 21.7 21.6 21.4];
    Hv = [61.2 59.6 58.5 58.1 56.5 55.2 54.4 53.8 53.3];


    % fitting a polynomial function for equilibrium curve
    E_curve_y = polyfit(xe,ye,8);

    xe_new = 0:0.001:1;
    for k = 1:length(xe_new) 
        ye_new(k) = polyval(E_curve_y,xe_new(k));
    end

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

    %x(1) = y
    %x(2) = Ln
    %x(3) = Vn
    %x(4) = HLn
    %x(5) = HVn
    
    % material balance equation
    F(1) = x(3) - x(2) - D;
    
    % component balance equation
    F(2) = x(3)*x(1) - x(2)*X - D*xd;
    
    % energybalance equation
    F(3) = x(3)*x(5) - x(2)*x(4) - qd_dash*D;
    
    % interpolation of Hl
    F(4) = x(4) - interp1(xe,Hl,X);
    
    % interpolation of Hv
    F(5) = x(5) - interp1(ye,Hv,x(1));
end


% stripping section equation
function F = Stripping_equations(x)

    global X

    % information given in the question
    zf = 0.42;
    xf = 0.42;
    xd = 0.97;
    xw = 0.01;
    R = 2.5;
    D = 4.2681;
    W = 5.7256;

    % data given in question
    xe = [0 0.08 0.18 0.25 0.49 0.65 0.79 0.91 1];
    ye = [0 0.28 0.43 0.51 0.73 0.83 0.9 0.96 1];
    Hl = [24.3 24.1 23.2 22.8 22.05 21.75 21.7 21.6 21.4];
    Hv = [61.2 59.6 58.5 58.1 56.5 55.2 54.4 53.8 53.3];

    
    % polynomial fit for equilibrium
    E_curve_y = polyfit(xe,ye,8);

    xe_new = 0:0.001:1;
    for k = 1:length(xe_new) 
        ye_new(k) = polyval(E_curve_y,xe_new(k));
    end

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

    %x(1) = y
    %x(2) = Ln
    %x(3) = Vn
    %x(4) = HLn
    %x(5) = HVn

    % material balance equation
    F(1) = x(2) - x(3) - W;
    
    % component balance equation
    F(2) = x(2)*X - x(3)*x(1) - W*xw;
    
    % energy balance equation
    F(3) = x(2)*x(4) - x(3)*x(5) - qw_dash*W;
    
    % interpolation of Hl from x
    F(4) = x(4) - interp1(xe,Hl,X);
    
    % interpolation of Hv from y 
    F(5) = x(5) - interp1(ye,Hv,x(1));
end