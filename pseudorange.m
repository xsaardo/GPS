function [h] = pseudorange(X)
sat1 = [3.5852; 2.07;        0];
sat2 = [2.9274; 2.9274;      0];
sat3 = [2.6612; 0;      3.1712];
sat4 = [1.4159; 0;      3.8904];

receiver_pos = X(1:3);
b = X(end);

h = [sqrt((receiver_pos - sat1)'*(receiver_pos - sat1));...
    sqrt((receiver_pos - sat2)'*(receiver_pos - sat2));...
    sqrt((receiver_pos - sat3)'*(receiver_pos - sat3));...
    sqrt((receiver_pos - sat4)'*(receiver_pos - sat4));] + b;

end