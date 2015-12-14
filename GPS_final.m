%% Final Variables:

% Errors in b
%%% error_gauss_b 
%%% error_steepest_b 

% Errors in position
%%% error_gauss_pos
%%% error_steepest_pos

% Loss Functions per iteration
%%% loss_gauss
%%% loss_steepest

% xopt_gauss and xopt_steepest are the final results of position and b
% given by each algorithm

%% Act 1
alpha = 1;
ER = 6370000;
epsilon = 0.0000001;

% Genie Information
receiver_pos = [1; 0; 0];
sat1 = [3.5852; 2.07;        0];
sat2 = [2.9274; 2.9274;      0];
sat3 = [2.6612; 0;      3.1712];
sat4 = [1.4159; 0;      3.8904];
b_actual = 2.354788068e-3;

% Generating Pseudoranges
yl = pseudorange([receiver_pos;b_actual]);

% Initial Conditions
s = [0.9331; 0.25; 0.258819];
b = 0;

%% Step 1: Initial Linearization
H  = [(s - sat1)'/sqrt((s-sat1)'*(s-sat1));...
    (s - sat2)'/sqrt((s-sat2)'*(s-sat2));...
    (s - sat3)'/sqrt((s-sat3)'*(s-sat3));...
    (s - sat4)'/sqrt((s-sat4)'*(s-sat4))];
H = [H ones(4,1)];

k = 1;
i = 100;
xopt = [zeros(4,100) [s;b]];

%% Steps 2 and 3: Gauss-Newton Simulation
% Termination Criteria
while(sqrt((xopt(:,end)-xopt(:,end-10))'*(xopt(:,end)-xopt(:,end-10))) > epsilon)  
    
    % Calculating the loss function for plotting
    loss_gauss(k) = 0.5*(pseudorange([s;b])-yl)'*(pseudorange([s;b])-yl); 

    % Linearization about the current estimate
    H  = [(s - sat1)'/sqrt((s-sat1)'*(s-sat1));...
        (s - sat2)'/sqrt((s-sat2)'*(s-sat2));...
        (s - sat3)'/sqrt((s-sat3)'*(s-sat3));...
        (s - sat4)'/sqrt((s-sat4)'*(s-sat4))];
    H = [H ones(4,1)];
    
    % Pseudorange at current estimate
    hl = [sqrt((s - sat1)'*(s - sat1));...
        sqrt((s - sat2)'*(s - sat2));...
        sqrt((s - sat3)'*(s - sat3));...
        sqrt((s - sat4)'*(s - sat4));] + b;
    
    delta_x  = alpha*H\(yl - hl);
    
    % Updating parameters
    xopt(:,i) = xopt(:,end) + delta_x;
    
    s = xopt(1:3,end);
    b = xopt(end,end);
    
    i = i+1;
    
    % Step 4: Error calculations per iteration
    error_gauss_pos(k) = sqrt((receiver_pos-s)'*(receiver_pos-s))*ER;    
    error_gauss_b(k) = abs(b_actual-b)*ER;
    k = k+1;
end
xopt_gauss = xopt(:,end)

%% Steepest Descent 
epsilon = 0.00000000000001;
alpha = 0.1;

% Initial Conditions
s = [0.9331; 0.25; 0.258819];
b = 0;

% Initial Linearization
H  = [(s - sat1)'/sqrt((s-sat1)'*(s-sat1));...
    (s - sat2)'/sqrt((s-sat2)'*(s-sat2));...
    (s - sat3)'/sqrt((s-sat3)'*(s-sat3));...
    (s - sat4)'/sqrt((s-sat4)'*(s-sat4))];

H = [H ones(4,1)];

k = 1;
i = 100;
xopt = [zeros(4,100) [s;b]];

%% Steps 2 and 3: Steepest Descent Simulation
% Termination Criteria
while(sqrt((xopt(:,end)-xopt(:,end-100))'*(xopt(:,end)-xopt(:,end-100))) > epsilon)
   
    % Calculating loss function for plotting
    loss_steepest(k) = 0.5*(pseudorange([s;b])-yl)'*(pseudorange([s;b])-yl);
   
    % Linearization about the current estimate
    H  = [(s - sat1)'/sqrt((s-sat1)'*(s-sat1));...
        (s - sat2)'/sqrt((s-sat2)'*(s-sat2));...
        (s - sat3)'/sqrt((s-sat3)'*(s-sat3));...
        (s - sat4)'/sqrt((s-sat4)'*(s-sat4))];
    H = [H ones(4,1)];
    
    % Calculating pseudorange at current estimate
    hl = [sqrt((s - sat1)'*(s - sat1));...
        sqrt((s - sat2)'*(s - sat2));...
        sqrt((s - sat3)'*(s - sat3));...
        sqrt((s - sat4)'*(s - sat4));] + b;
    
    delta_x  = alpha*H'*(yl - hl);
   
    xopt(:,i) = xopt(:,end) + delta_x;
    
    s = xopt(1:3,end);
    b = xopt(end,end);
    
    i = i+1;
    
    % Step 4: Error calculations per iteration
    error_steepest_pos(k) = sqrt((receiver_pos-s)'*(receiver_pos-s))*ER;    
    error_steepest_b(k) = abs(b_actual-b)*ER;
    k = k+1;
end
xopt_steepest = xopt(:,end)

%% Figures
plot(error_gauss_pos,'--');
hold on;
plot(error_gauss_b);
legend({'Position Error S','Clock Bias Error b'},'Interpreter','latex');
xlabel('Number of Iterations','Interpreter','latex');
ylabel('Error (meters)','Interpreter','latex');
title('GPS Error vs Gauss-Newton Iteration','Interpreter','latex');
figure;
plot(loss_gauss);
xlabel('Number of Iterations','Interpreter','latex');
ylabel('Loss', 'Interpreter', 'latex');
title('Loss Function Per Iteration','Interpreter','latex');

figure;
plot(error_steepest_pos,'--');
hold on;
plot(error_steepest_b);
legend({'Position Error S','Clock Bias Error b'},'Interpreter','latex');
xlabel('Number of Iterations','Interpreter','latex');
ylabel('Error (meters)','Interpreter','latex');
title('GPS Error vs Gauss-Newton Iteration','Interpreter','latex');
figure;
semilogy(loss_steepest);
xlabel('Number of Iterations','Interpreter','latex');
ylabel('Loss', 'Interpreter', 'latex');
title('Loss Function Per Iteration','Interpreter','latex');