%% Initialization
clc; clear all; close all; format short g; tic

% load_results = {results file};

% BCs and Algorithm parameters
w_nu = 1e5;
w_Delta_bar = 1e-3;
w_Delta_sigma = 1e-1;
nu_tol = 1e-10;
Delta_tol = 1e-3;
N_iter_max = 15;
K = 50;
% Initial               % Final
r_i = [4 4 0]';         r_f = [0 0 0]';
v_i = [0 0 2]';         v_f = [-1e-1 0 0]';
w_i = [0 0 0]';         w_f = [0 0 0]';
q_i = [1 0 0 0]';       q_f = [1 0 0 0]';
tf_guess = 1;
alpha_mdot = 0.05; % Not provided in paper.

% Simulation Parameters
g = [-1 0 0]';
m_wet = 2;
m_dry = 1;
T_min = 0.3;
T_max = 5;
delta_max = deg2rad(20);
theta_max = deg2rad(90);
gamma_gs = deg2rad(0);
omega_max = deg2rad(60);
J_B = 1e-2*eye(3);
r_T_B = -1e-2*[1 0 0];
H_23 = [0 1 0; 0 0 1]';
H_q = [0 0 1 0; 0 0 0 1];

% Initial state
x_0 = [m_wet r_i' v_i' q_i' w_i'];

% Symbolic variables for dynamics
m = sym('m', 'real');               % Mass
r = sym('r', [3 1], 'real');        % Position
v = sym('v', [3 1], 'real');        % Velocity
q = sym('q', [4 1], 'real');        % Attitude (quaternion)
w = sym('w', [3 1], 'real');        % Angular velocity
x = [m r' v' q' w']';               % State vector
u = sym('T', [3 1], 'real');        % Thrust
sigma = sym('sigma', 'real');       % Time dilation coefficient

C_B_I = DCM(q);
C_I_B = C_B_I';

mdot = -alpha_mdot * norm(u);
rdot = v;
vdot = 1/m * C_I_B * u + g;
qdot = 1/2 * Omega(w) * q;
wdot = inv(J_B) * (skew(r_T_B)*u - skew(w)*J_B*w);
xdot = simplify([mdot rdot' vdot' qdot' wdot']');

A = simplify(sigma*jacobian(xdot,x));
B = simplify(sigma*jacobian(xdot,u));
z = simplify(-A*x - B*u);

A_func = matlabFunction(A,'Vars',{[m r' v' q' w'], u', sigma});
B_func = matlabFunction(B,'Vars',{[m r' v' q' w'], u', sigma});
Sigma_func = matlabFunction(xdot,'Vars',{[m r' v' q' w'], u', sigma});
z_func = matlabFunction(z,'Vars',{[m r' v' q' w'], u', sigma});

clear m r v q w x u sigma A B z
%% SCvx Loop

% Define variables as cell{SCvx iteration}(time)
K_set = 0:(K-1);
dt = 1/(K-1);
% Initialization
try 
    load(load_results,'x','u','sigma')
    n_start = numel(x) + 1;
catch
    for i = 1:K
        k = K_set(i);
        tau(i) = k/(K-1);
        alpha_1 = (K-k)/K;
        alpha_2 = k/K;
    
        % Guess reference state
        m = alpha_1*m_wet + alpha_2*m_dry;
        r = alpha_1*r_i;
        v = alpha_1*v_i + alpha_2*v_f;
        q = q_i;
        w = w_i;
        x{1}(i,:) = [m r' v' q' w'];
    
        % Guess reference control
        u{1}(i,:) = -m * g;
        sigma{1} = tf_guess;

        n_start = 2;
    end
end
%%%%%%%%%% SCVX LOOP %%%%%%%%%%
for n = n_start:n_start+N_iter_max
    % Set up transition matrices using trajectory from n-1
    for i = 1:K-1
        A{i} = A_func(x{n-1}(i,:), u{n-1}(i,:), sigma{n-1});
        B{i} = B_func(x{n-1}(i,:), u{n-1}(i,:), sigma{n-1});
        Sigma{i} = Sigma_func(x{n-1}(i,:), u{n-1}(i,:), sigma{n-1});
        z{i} = z_func(x{n-1}(i,:), u{n-1}(i,:), sigma{n-1});

        alpha_k = @(xi) (tau(i+1) - xi)/(tau(i+1) - tau(i));
        beta_k = @(xi) (xi - tau(i))/(tau(i+1) - tau(i));

        Phi_func = @(xi) expm(A{i}*(tau(i+1)-xi));
        Abar{i} = Phi_func(tau(i));
        Bbar_integrand = @(xi) Phi_func(xi) * B{i} * alpha_k(xi);
        Bbar{i} = integral(Bbar_integrand, tau(i), tau(i+1), 'ArrayValued', 1);
        Cbar_integrand = @(xi) Phi_func(xi) * B{i} * beta_k(xi);
        Cbar{i} = integral(Cbar_integrand, tau(i), tau(i+1), 'ArrayValued', 1);
        Sigmabar_integrand = @(xi) Phi_func(xi) * Sigma{i};
        Sigmabar{i} = integral(Sigmabar_integrand, tau(i), tau(i+1), 'ArrayValued', 1);
        zbar_integrand = @(xi) Phi_func(xi) * z{i};
        zbar{i} = integral(zbar_integrand, tau(i), tau(i+1), 'ArrayValued', 1);
    end
    
    %%%%%%%%%% CONVEX PROBLEM START %%%%%%%%%%
    cvx_begin
        %%% Initialization %%%
        variable u_cvx(K,3)                         % Thrust control
        variable sigma_cvx(1) nonnegative           % Time dilation control
        variable nu_cvx(K,14)                       % Virtual state control
        variable Delta_bar_cvx(K) nonnegative       % State/control trust region
        variable Delta_sigma_cvx(1) nonnegative     % Time dilation trust region
        expression x_cvx(K,14)                      % State expression
        expression delta_x_cvx(K,14)                % State trust region distance
        expression delta_u_cvx(K,3)                 % Thrust control trust region distance
        expression delta_sigma_cvx(1)               % Time dilation control trust region distance

        %%% Expressions %%%
        % State
        x_cvx(1,:) = x_0;
        for k = 1:K-1
            x_cvx(k+1,:) = Abar{k}*x_cvx(k,:)' + Bbar{k}*u_cvx(k,:)' + Cbar{k}*u_cvx(k+1,:)' + Sigmabar{k}*sigma_cvx + zbar{k} + nu_cvx(k,:)';
        end

        % Control
        for k = 1:K
            B_g(k,:) = u{n-1}(k,:) / norm(u{n-1}(k,:),2);
        end
        
        % Trust regions
        for k = 1:K
            delta_x_cvx(k,:) = x_cvx(k,:) - x{n-1}(k,:);
            delta_u_cvx(k,:) = u_cvx(k,:) - u{n-1}(k,:);
        end
        delta_sigma_cvx = sigma_cvx - sigma{n-1};
        
        %%% Constraints %%%
        for k = 1:K
            % State constraints
            m_dry <= x_cvx(k,1)
            tan(gamma_gs)*norm(H_23'*x_cvx(k,2:4)',2) <= dot([1 0 0]',x_cvx(k,2:4))
            cos(theta_max) <= 1 - 2*sum_square_abs(H_q*x_cvx(k,8:11)')
            norm(x_cvx(k,12:14),2) <= omega_max
            
            % Control constraints
            T_min <= B_g(k,:)*u_cvx(k,:)'
            norm(u_cvx(k,:),2) <= T_max
            cos(delta_max)*norm(u_cvx(k,:),2) <= dot([1 0 0], u_cvx(k,:))

            % Trust region constraint
            dot(delta_x_cvx(k,:),delta_x_cvx(k,:)) + dot(delta_u_cvx(k,:),delta_u_cvx(k,:)) <= Delta_bar_cvx(k)
        end
        norm(delta_sigma_cvx,1) <= Delta_sigma_cvx
        
        % Boundary Conditions
        x_cvx(end,2:4) == r_f';
        x_cvx(end,5:7) == v_f';
        x_cvx(end,8:11) == q_f';
        x_cvx(end,12:14) == w_f';

        %%% Objective %%%
        minimize(sigma_cvx + w_nu*norm(nu_cvx,1) + w_Delta_bar*norm(Delta_bar_cvx,2) + w_Delta_sigma*norm(Delta_sigma_cvx,1))

    cvx_end
    %%%%%%%%%% CONVEX PROBLEM END %%%%%%%%%%

    % Store results for next iteration
    x{n} = x_cvx;
    u{n} = u_cvx;
    sigma{n} = sigma_cvx;
    nu{n} = nu_cvx;
    
    disp("Iteration:")
    disp(n-1)
    disp("Objectives:")
    disp([sigma_cvx, norm(nu_cvx,1), norm(Delta_bar_cvx,2), norm(Delta_sigma_cvx,1)])

    % Tolerance breaks
    if norm(Delta_bar_cvx,2) <= Delta_tol && norm(nu_cvx,1) <= nu_tol
        disp("Tolerance condition met.")
        break
    end

    clear Abar Bbar Cbar Sigmabar zbar
end
toc
toc_SCvx = toc;
save results
%%%%%%%%%% SCVX LOOP %%%%%%%%%%
%% Plotting
try
    load('results.mat')
catch
    disp('No results file.')
end
close all
clf
hold on
az = -145;
el = 30;
view(az,el)
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
grid on
set(gcf,'color','w');

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
axis equal
grid on
gif_axis = axis;
view(az,el)
annotation('textbox','Position',[0.75, 0.75, 0.1, 0.1],'String', ...
    "||\nu||_1: "+norm(nu{end},1)+newline ...
    +"\sigma: "+sigma{end})

%% GIF maker
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
annotation('textbox','Position',[0.75, 0.75, 0.1, 0.1],'String', ...
    "||\nu||_1: "+norm(nu{end},1)+newline ...
    +"\sigma: "+sigma{end})
% gif('SCvx.gif','DelayTime',1/12,'overwrite',true,'resolution',300)
fps = 1.5 * K/sigma{n};
gif('SCvx.gif','DelayTime',1/fps,'overwrite',true)
for k = 1:K
    clf
    hold on
    C_I_B = DCM(x{n}(k,8:11))';
    attitude(k,:) = C_I_B*[1 0 0]';
    attitude(k,:) = .5 * attitude(k,:) / norm(attitude(k,:));
    point(k,:) = x{n}(k,2:4) + (C_I_B*r_T_B')';

    quiver3(point(k,2),point(k,3),point(k,1),attitude(k,2),attitude(k,3),attitude(k,1),'ShowArrowHead','off','Color','b','LineWidth',1.5)
    quiver3(point(k,2),point(k,3),point(k,1),-thrust(k,2),-thrust(k,3),-thrust(k,1),'ShowArrowHead','off','Color','r','LineWidth',1.5)
    plot3(x{n}(1:k,3),x{n}(1:k,4),x{n}(1:k,2),':k')
    view(az,el)
    axis equal    
    axis(gif_axis)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
    set(gcf,'color','w');
    annotation('textbox','Position',[0.75, 0.75, 0.1, 0.1],'String', ...
    "||\nu||_1: "+norm(nu{end},1)+newline ...
    +"\sigma: "+sigma{end})
    gif
end
for i = 1:fps
    clf
    hold on
    quiver3(point(k,2),point(k,3),point(k,1),attitude(k,2),attitude(k,3),attitude(k,1),'ShowArrowHead','off','Color','b','LineWidth',1.5)
    plot3(x{n}(:,3),x{n}(:,4),x{n}(:,2),':k')
    view(az,el)
    axis equal    
    axis(gif_axis)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    grid on
    set(gcf,'color','w');
    annotation('textbox','Position',[0.75, 0.75, 0.1, 0.1],'String', ...
    "||\nu||_1: "+norm(nu{end},1)+newline ...
    +"\sigma: "+sigma{end})
    gif
end
%% Functions
function C = DCM(q)
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    C = [1-2*(q2^2+q3^2), 2*(q1*q2+q0*q3), 2*(q1*q3-q0*q2);
         2*(q1*q2-q0*q3), 1-2*(q1^2+q3^2), 2*(q2*q3+q0*q1);
         2*(q1*q3+q0*q2), 2*(q2*q3-q0*q1), 1-2*(q1^2+q2^2)];
end

function xi_skew = skew(xi)
    xi_x = xi(1);
    xi_y = xi(2);
    xi_z = xi(3);
    
    xi_skew = [0,       -xi_z,  xi_y;
               xi_z,    0,      -xi_x;
              -xi_y,    xi_x,   0];
end

function Omega = Omega(xi)
    xi_x = xi(1);
    xi_y = xi(2);
    xi_z = xi(3);
    
    Omega = [0,     -xi_x,  -xi_y,  -xi_z;
             xi_x,  0,      xi_z,   -xi_y;
             xi_y,  -xi_z,  0,      xi_x;
             xi_z,  xi_y,   -xi_x,  0];
end