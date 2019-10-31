armParamHWB;  % load parameters

% instantiate arm, and reference input classes 
arm = armDynamics(P);  
addpath('../hw_a'); reference = signalGenerator(0.01, 0.02);  
addpath('../hw_a'); torque = signalGenerator(0.2, 0.05);

% instantiate the data plots and animation
addpath('../hw_a'); dataPlot = dataPlotter(P);
addpath('../hw_a'); animation = armAnimation(P);

% main simulation loop
t = P.t_start;  % time starts at t_start
while t < P.t_end  
    % Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot;
    while t < t_next_plot % updates control and dynamics at faster simulation rate
        r = reference.square(t);
        u = torque.square(t);  % Calculate the input force
        y = arm.update(u);  % Propagate the dynamics
        t = t + P.Ts; % advance time by Ts
    end
    % update animation and data plots
    animation.update(arm.state);
    dataPlot.update(t, r, arm.state, u);
end
