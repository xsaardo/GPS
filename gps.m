clc; clear; close all;

%% Act 1
alpha = 1;
sigma = 0.000400;
m = 256;

% Genie Information
receiver_pos = [1; 0; 0];
sat1 = [3.5852; 2.07;        0];
sat2 = [2.9274; 2.9274;      0];
sat3 = [2.6612; 0;      3.1712];
sat4 = [1.4159; 0;      3.8904];

b0 = 2.354788068e-3;

% Pseudoranges
y = zeros(4,1);
for i = 1:m 
y = y + [sqrt((receiver_pos - sat1)'*(receiver_pos - sat1));...
    sqrt((receiver_pos - sat2)'*(receiver_pos - sat2));...
    sqrt((receiver_pos - sat3)'*(receiver_pos - sat3));...
    sqrt((receiver_pos - sat4)'*(receiver_pos - sat4));] + b0 + normrnd(0,sigma);
end
yavg = y/m;
yl = yavg;

% Initial Conditions
x0 = [0.9331; 0.25; 0.258819];
b0 = 0;

% Initial Linearization
H  = [(x0 - sat1)'/sqrt((x0-sat1)'*(x0-sat1));...
    (x0 - sat2)'/sqrt((x0-sat2)'*(x0-sat2));...
    (x0 - sat3)'/sqrt((x0-sat3)'*(x0-sat3));...
    (x0 - sat4)'/sqrt((x0-sat4)'*(x0-sat4))];

H = [H ones(4,1)];

xopt = [x0;b0];

% Gauss-Newton
for i = 1:100000
    s = xopt(1:3);
    b = xopt(end);
    
    H  = [(s - sat1)'/sqrt((s-sat1)'*(s-sat1));...
        (s - sat2)'/sqrt((s-sat2)'*(s-sat2));...
        (s - sat3)'/sqrt((s-sat3)'*(s-sat3))
        (s - sat4)'/sqrt((s-sat4)'*(s-sat4))];
    
    H = [H ones(4,1)];
    
    hl = [sqrt((s - sat1)'*(s - sat1));...
        sqrt((s - sat2)'*(s - sat2));...
        sqrt((s - sat3)'*(s - sat3));...
        sqrt((s - sat4)'*(s - sat4));] + b;
    
    delta_x  = alpha*inv(H)*(yl - hl);
    
    xopt = xopt + delta_x;
end