clc; clear all; close all; format short g; tic
try
    load('results.mat')
catch
    error('No results file.')
end

% Set n to the SCvx iteration you want to plot.
n = numel(x);

% GIF parameters
filename = 'SCvx.gif';
fps = 12;

% Plot parameters
az = -145;
el = 30;
%% Plot Trajectory
set_canvas(az, el)
annotate_plot(n,nu,sigma)

plot3(x{n}(:,3),x{n}(:,4),x{n}(:,2),':k')
thrust = u{n} / T_max / 4;
for k = 1:K
    C_I_B = DCM(x{n}(k,8:11))';
    attitude(k,:) = C_I_B*[1 0 0]';
    attitude(k,:) = .5 * attitude(k,:) / norm(attitude(k,:));
    point(k,:) = x{n}(k,2:4) + (C_I_B*r_T_B')';

    quiver3(point(k,2),point(k,3),point(k,1),attitude(k,2),attitude(k,3),attitude(k,1),'ShowArrowHead','off','Color','b','LineWidth',1.5)
    quiver3(point(k,2),point(k,3),point(k,1),-thrust(k,2),-thrust(k,3),-thrust(k,1),'ShowArrowHead','off','Color','r','LineWidth',1.5)
end
gif_axis = axis;

%% Plot 
figure
time = linspace(0,sigma{n},K);

subplot(6,1,1)
hold on
plot(time, [x{n}(:,3),x{n}(:,4),x{n}(:,2)],'.-')
grid on
title("Position")
xlim([0,sigma{n}])

subplot(6,1,2)
hold on
plot(time, [x{n}(:,6),x{n}(:,7),x{n}(:,5)],'.-')
grid on
title("Velocity")
xlim([0,sigma{n}])

subplot(6,1,3)
hold on
plot(time, acosd(1-2*(x{n}(:,9).^2 + x{n}(:,10).^2)),'.-')
plot(time, rad2deg(theta_max)*ones(1,K),'--k')
grid on
title("Tilt angle")
xlim([0,sigma{n}])

subplot(6,1,4)
hold on
plot(time, rad2deg(sqrt(x{n}(:,12).^2 + x{n}(:,13).^2 + x{n}(:,14).^2)),'.-')
plot(time, rad2deg(omega_max)*ones(1,K),'--k')
grid on
title("Angular rate")
xlim([0,sigma{n}])

subplot(6,1,5)
hold on
plot(time, sqrt(u{n}(:,2).^2 + u{n}(:,3).^2 + u{n}(:,1).^2),'.-')
plot(time, T_max*ones(1,K),'--k')
plot(time, T_min*ones(1,K),'--k')
grid on
title("Thrust")
xlim([0,sigma{n}])


subplot(6,1,6)
hold on
plot(time, x{n}(:,1),'.-')
plot(time, m_wet*ones(1,K),'--k')
plot(time, m_dry*ones(1,K),'--k')
grid on
title("Mass")
xlim([0,sigma{n}])
%% GIF maker
figure
set_canvas(az, el, gif_axis)
annotate_plot(1,nu,sigma)
gif(filename,'DelayTime',1/fps,'overwrite',true)
thrust = u{n} / T_max / 4;
for k = 1:K
    set_canvas(az, el, gif_axis)
    annotate_plot(n,nu,sigma)

    C_I_B = DCM(x{n}(k,8:11))';
    attitude(k,:) = C_I_B*[1 0 0]';
    attitude(k,:) = .5 * attitude(k,:) / norm(attitude(k,:));
    point(k,:) = x{n}(k,2:4) + (C_I_B*r_T_B')';

    quiver3(point(k,2),point(k,3),point(k,1),attitude(k,2),attitude(k,3),attitude(k,1),'ShowArrowHead','off','Color','b','LineWidth',1.5)
    quiver3(point(k,2),point(k,3),point(k,1),-thrust(k,2),-thrust(k,3),-thrust(k,1),'ShowArrowHead','off','Color','r','LineWidth',1.5)
    plot3(x{n}(1:k,3),x{n}(1:k,4),x{n}(1:k,2),':k')
    gif
end
for i = 1:fps
    set_canvas(az, el, gif_axis)
    annotate_plot(n,nu,sigma)

    quiver3(point(k,2),point(k,3),point(k,1),attitude(k,2),attitude(k,3),attitude(k,1),'ShowArrowHead','off','Color','b','LineWidth',1.5)
    plot3(x{n}(1:k,3),x{n}(1:k,4),x{n}(1:k,2),':k')
    gif
end

%% Functions
function [] = set_canvas(varargin)
    switch numel(varargin)
        case 2
            [az, el] = varargin{:};
        case 3
            [az, el, gif_axis] = varargin{:};
    end

    clf
    hold on
    view(az,el)
    axis equal
    if exist('gif_axis')
        axis(gif_axis)
    end
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
    set(gcf,'color','w');
end
function [] = annotate_plot(n,nu,sigma)
        title("SCvx Iteration: "+(n-1))
        annotation('textbox','Position',[0.75, 0.75, 0.1, 0.1],'String', ...
        "||\nu||_1: "+norm(nu{n},1)+newline ...
        +"\sigma: "+sigma{n})
end
function C = DCM(q)
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    C = [1-2*(q2^2+q3^2), 2*(q1*q2+q0*q3), 2*(q1*q3-q0*q2);
         2*(q1*q2-q0*q3), 1-2*(q1^2+q3^2), 2*(q2*q3+q0*q1);
         2*(q1*q3+q0*q2), 2*(q2*q3-q0*q1), 1-2*(q1^2+q2^2)];
end