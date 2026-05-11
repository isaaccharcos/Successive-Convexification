clc; clear all; close all; format short g; tic
try
    load('results.mat')
catch
    error('No results file.')
end

% GIF parameters
filename = 'SCvx Iterations.gif';
fps = 12;

% Plot parameters
az = -145;
el = 30;
%% Plotting
hold on
for n = 1:N_iter_max+1
    plot3(x{n}(:,3),x{n}(:,4),x{n}(:,2))
    thrust = u{n} / T_max / 4;
    for k = 1:K
        C_I_B = DCM(x{n}(k,8:11))';
        attitude(k,:) = C_I_B*[1 0 0]';
        attitude(k,:) = .5 * attitude(k,:) / norm(attitude(k,:));
        point(k,:) = x{n}(k,2:4) + (C_I_B*r_T_B')';
    
        quiver3(point(k,2),point(k,3),point(k,1),attitude(k,2),attitude(k,3),attitude(k,1),'ShowArrowHead','off','Color','b','LineWidth',1.5)
        quiver3(point(k,2),point(k,3),point(k,1),-thrust(k,2),-thrust(k,3),-thrust(k,1),'ShowArrowHead','off','Color','r','LineWidth',1.5)
    end
end
axis equal
view(az,el)
gif_axis = axis;

%% GIF maker
set_canvas(az, el, gif_axis)
annotate_plot(1,nu,sigma)
gif(filename,'DelayTime',1/fps,'overwrite',true)
for n = 1:N_iter_max+1
    thrust = u{n} / T_max / 4;
    for k = 1:K
        set_canvas(az, el, gif_axis)
        annotate_plot(n,nu,sigma)

        C_I_B = DCM(x{n}(k,8:11))';
        attitude(k,:) = C_I_B*[1 0 0]';
        attitude(k,:) = .5 * attitude(k,:) / norm(attitude(k,:));
        point(k,:) = x{n}(k,2:4) + (C_I_B*r_T_B')';

        for j = 1:n-1
            plot3(x{j}(:,3),x{j}(:,4),x{j}(:,2),':','Color',.75*ones(1,3))
        end
        quiver3(point(k,2),point(k,3),point(k,1),attitude(k,2),attitude(k,3),attitude(k,1),'ShowArrowHead','off','Color','b','LineWidth',1.5)
        quiver3(point(k,2),point(k,3),point(k,1),-thrust(k,2),-thrust(k,3),-thrust(k,1),'ShowArrowHead','off','Color','r','LineWidth',1.5)
        plot3(x{n}(1:k,3),x{n}(1:k,4),x{n}(1:k,2),':k')
        gif
    end
    for i = 1:fps
        set_canvas(az, el, gif_axis)
        annotate_plot(n,nu,sigma)

        for j = 1:n-1
            plot3(x{j}(:,3),x{j}(:,4),x{j}(:,2),':','Color',.75*ones(1,3))
        end
        quiver3(point(k,2),point(k,3),point(k,1),attitude(k,2),attitude(k,3),attitude(k,1),'ShowArrowHead','off','Color','b','LineWidth',1.5)
        plot3(x{n}(1:k,3),x{n}(1:k,4),x{n}(1:k,2),':k')
        gif
    end
end
%% Functions
function [] = set_canvas(az, el, gif_axis)
    clf
    hold on
    view(az,el)
    axis equal
    axis(gif_axis)
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