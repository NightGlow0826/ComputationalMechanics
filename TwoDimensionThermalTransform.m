clear
clc

% MESH = 8;


Problem = 2



meshes = [8, 24, 96];

% 第一个subplot，绘制温度曲线
figure;
subplot(1,2,1);
for i = 1:3
    [r, temp] = get_temperature(meshes(i), Problem);
    plot(r, temp);
    hold on;
    if i == 3
        plot(r, get_ana_temp(r, Problem));
        hold on;
    end
end
legend('8', '24', '96', 'ana');
title = sprintf('Situation %d, T~R', Problem);
xlabel('r');
ylabel('T');
hold off;

% 第二个subplot，绘制残差曲线
subplot(1,2,2);
for i = 1:3
    [r, temp] = get_temperature(meshes(i), Problem);
    ana = get_ana_temp(r, Problem)';
    plot(r, ana - temp);
    hold on;
end
hold off;
legend('8', '24', '96');
title = sprintf('Situation %d, Res T~R', Problem);
xlabel('r');
ylabel('Res T');

% FUNCTIONS
% FUNCTIONS
% FUNCTIONS
% FUNCTIONS
function [r, temp] = get_temperature(MESH, Problem)
a = 0.5; b =1.5; k = 1.;

node_coords = get_coords(MESH, a, b);
triangle_node = get_triangles_node(MESH);
nodes_num = (MESH + 1) .^2;
trans_mat = zeros(nodes_num, nodes_num);
for t = 1: length(triangle_node)
    nodes = triangle_node(t, :);
    node1 = nodes(1);
    node2 = nodes(2);
    node3 = nodes(3);
    x1 = node_coords(node1, 1);
    y1 = node_coords(node1, 2);
    x2 = node_coords(node2, 1);
    y2 = node_coords(node2, 2);
    x3 = node_coords(node3, 1);
    y3 = node_coords(node3, 2);
    b1 = y2 - y3;
    b2 = y3 - y1;
    b3 = y1 - y2;
    c1 = -(x2 - x3);
    c2 = -(x3 - x1);
    c3 = -(x1 - x2);
    S = abs(0.5 * det([1, 1, 1; x1, x2, x3;y1, y2, y3]));
    B = 1/ (2*S) * [b1, b2, b3; c1, c2, c3];
    trans_mat([node1, node2, node3], [node1, node2, node3]) = trans_mat([node1, node2, node3], [node1, node2, node3]) + k*S*(B'*B);
end

Q = zeros(nodes_num, 1);
T = zeros(nodes_num, 1);
inner_constraints = 1:MESH + 1;
outer_constraints = 1 + (MESH + 1)*MESH:(MESH + 1).^2;
constraints = [inner_constraints, outer_constraints];
unconstraints = setdiff(1:nodes_num, constraints);
un_inner_constraints = setdiff(1:nodes_num, inner_constraints);
un_outter_constraints = setdiff(1:nodes_num, outer_constraints);
% For Problem 1:
if Problem == 1
    T(inner_constraints) = 400;
    T(outer_constraints) = 300;
    Q(constraints) = trans_mat(constraints, constraints) * T(constraints);
    % since Q of the middle nodes == 0
    size(trans_mat(constraints, unconstraints));
    T(unconstraints) = -inv(trans_mat(unconstraints, unconstraints)) * trans_mat(unconstraints, constraints) * T(constraints);
    
else
    
    % For Problem 2:
    Q_side = -100 * 0.5 * b * pi / MESH;
    T(inner_constraints) = 400;
    Q([outer_constraints(1), outer_constraints(end)]) = [Q_side / 2, Q_side / 2];
    Q(outer_constraints(2:end-1)) = Q_side;
    size(Q(un_inner_constraints));
    T(un_inner_constraints) = inv(trans_mat(un_inner_constraints, un_inner_constraints)) * (Q(un_inner_constraints) - trans_mat(un_inner_constraints, inner_constraints) * T(inner_constraints));
end
temp = T(1:MESH+1:end);
delta_r = (b - a) / MESH;
r = a + delta_r * [0:MESH];

end



function node_coords = get_coords(MESH, a, b);
node_coords = [];
delta_r = (b - a) / MESH;
delta_theta = pi / 2 / MESH;
for i = 0:MESH
    r = a + delta_r * i;
    for j = 0:MESH
        theta = 0 + delta_theta * j;
        node_coords = [node_coords; r * cos(theta), r * sin(theta)];
    end
end
end


function triangle_node = get_triangles_node(MESH)
triangle_node = [];
for i = 1:MESH
    for j = 1:MESH
        for k = 1:2
            node = position_2_node(i, j, MESH);
            if k ==1
                % triangle_node = [triangle_node; [i, i+1, ]];
                triangle_node = [triangle_node; [node, node + 1, node + 1 + MESH + 1]];
            else
                triangle_node = [triangle_node; [node, node + MESH + 1, node + MESH + 2]];
                
            end
        end
    end
end
end

function node = position_2_node(i, j, MESH)
node = (MESH + 1) * (i - 1) + j;
end

function ana_temp = get_ana_temp(r, Problem)
a = 0.5;
b = 1.5;
if Problem == 1
    ana_temp = 400 - (400-300)*(log(r)-log(a))/(log(b)-log(a));
else
    ana_temp = 400 - 100*b*log(r/a);
end
end
