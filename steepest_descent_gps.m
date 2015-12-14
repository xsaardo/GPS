clc; clear; close all;

%% Act 1
alpha = .1;
ER = 6370000;
epsilon = 0.00000000000001;

% Genie Information
receiver_pos = [1; 0; 0];
sat1 = [3.5852; 2.07;        0];
sat2 = [2.9274; 2.9274;      0];
sat3 = [2.6612; 0;      3.1712];
sat4 = [1.4159; 0;      3.8904];
b_actual = 2.354788068e-3;

% Pseudoranges
yl = pseudorange([receiver_pos;b_actual]);

%% Act 2
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

% Steepest Descent
while(sqrt((xopt(:,end)-xopt(:,end-100))'*(xopt(:,end)-xopt(:,end-100))) > epsilon)
    H  = [(s - sat1)'/sqrt((s-sat1)'*(s-sat1));...
        (s - sat2)'/sqrt((s-sat2)'*(s-sat2));...
        (s - sat3)'/sqrt((s-sat3)'*(s-sat3));...
        (s - sat4)'/sqrt((s-sat4)'*(s-sat4))];
    
    H = [H ones(4,1)];
    
    hl = [sqrt((s - sat1)'*(s - sat1));...
        sqrt((s - sat2)'*(s - sat2));...
        sqrt((s - sat3)'*(s - sat3));...
        sqrt((s - sat4)'*(s - sat4));] + b;
    
    delta_x  = alpha*H'*(yl - hl);
   
    xopt(:,i) = xopt(:,end) + delta_x;
    
    s = xopt(1:3,end);
    b = xopt(end,end);
    
    i = i+1;
    
    % Error calculations per iteration
    error_pos(k) = sqrt((receiver_pos-s)'*(receiver_pos-s))*ER;    
    error_b(k) = abs(b_actual-b)*ER;
    k = k+1;
    loss(k) = 0.5*(pseudorange([s;b])-yl)'*(pseudorange([s;b])-yl);
end

plot(error_pos,'--');
hold on;
plot(error_b);
legend({'Position Error S','Clock Bias Error b'},'Interpreter','latex');
xlabel('Number of Iterations','Interpreter','latex');
ylabel('Error (meters)','Interpreter','latex');
title('GPS Error vs Steepest Descent Iteration','Interpreter','latex');
figure;
semilogy(loss);
xlabel('Number of Iterations','Interpreter','latex');
ylabel('Loss', 'Interpreter', 'latex');
title('Loss Function Per Iteration','Interpreter','latex');