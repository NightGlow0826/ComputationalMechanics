clc;
% Some global variables
E = 100; G = 38.5;
A = 0.2;
L = 1;
Iz = 0.003; Iy = 0.003; Ix = 1; J = 0.006;
Fx = -3194.0; Fy = -856.0;
% Create the local stiffness matrix
% k_local = get_local_stiffness(E, G, A, J, L, Iz, Iy, Ix);
% Create the transformation matrix
% T stands for a 12^2 matrix, which is composed of 4 3*3 mat in the diagonal

% r = get_r(1, 2, 3, 4, 5, 6)

node_coords = [58.0, 38.0, 0.0;
    48.0, 38., 0.0;
    31.0, 38., 0.0;
    17.0, 38., 22.0;
    0.0 ,38.0, 24.0;
    58.0, 38.0, 42.0;
    48.0, 38.0, 42.0;
    36.0, 38.0, 70.0;
    0.0 ,38.0, 75.0;
    58.0, 17.0, 42.0;
    58.0, 17.0, 0.0;
    0.0 ,17.0, 0.0;
    0.0 ,17.0, 24.0;
    18.0, 0.0, 72.0;
    0.0 ,0.0 37.5];


% Create the coordinates and the connections
symmetry_node = [node_coords(:, 1), -node_coords(:, 2), node_coords(:, 3)]; % lets k' them with k + 15
node_coords = [node_coords; symmetry_node];
node_coords = node_coords(1:28, :); % Since we do not want to duplicate point 14,15

node_link = [1, 2; 1, 11; 2, 3; 3, 4; 4,5; 5, 9; 8, 9; 8, 14; 9, 14; 2, 7; 6, 7; 6, 10;
    10, 11; 11, 12; 7, 8; 12, 13; 5, 15; 9, 15; 13, 5];
sym_node = node_link + 15;
self_connect = [11, 10, 8, 13, 9];
self_connect = transpose([self_connect; self_connect + 15]);
node_link = [node_link; sym_node; self_connect];
node_link(node_link == 30) = 15;
node_link(node_link == 29) = 14;
length(node_link);
node_link;


global_force_dof = 6 * length(node_coords);
force = zeros([global_force_dof, 1]);
force(1) = Fx; force(2) = Fy;


stiffness = zeros([global_force_dof, global_force_dof]);
% Compose the stiffness
for i = 1:length(node_link)
    node1 = node_link(i, 1);
    node2 = node_link(i, 2);
    [node1, node2];
    x1 = node_coords(node1, 1);
    y1 = node_coords(node1, 2);
    z1 = node_coords(node1, 3);
    x2 = node_coords(node2, 1);
    y2 = node_coords(node2, 2);
    z2 = node_coords(node2, 3);
    L = sqrt((x2 - x1).^2 + (y2 - y1).^2 + (z2 - z1).^2);
    k_local = get_local_stiffness(E, G, A, J, L, Iz, Iy, Ix);
    % find the k_th dof coordinate
    a = -5:0;
    dofs = [6*node1 + a, 6*node2 + a];
    r = get_r(x1, y1, z1, x2, y2, z2);
    % a1 = [x1, y1, z1]
    % z2 = [x2, y2, z2]
    k_universal = r'* k_local * r;
    stiffness(dofs, dofs) = stiffness(dofs, dofs) + k_universal;
end
bounded_dofs = [];
boundaries = [11, 12, 26, 27];
for i = 1:length(boundaries)
    % bounded_dofs = [bounded_dofs, 6*boundaries(i) - 5:6*boundaries(i) - 0];
    bounded_dofs = [bounded_dofs, 6*boundaries(i) - 5:6*boundaries(i) - 3];
end
bounded_dofs;

displacements = solution(global_force_dof, bounded_dofs, stiffness, force);

nodes = [1, 2, 6, 7, 10, 11];
a = -5:0;
node_seq = [];
% displacements(1:168)
for i = 1:length(nodes)
    node_seq = [node_seq, nodes(i)*6 + a];
end
node_seq;
ans = reshape(displacements(node_seq), 6,  6)


forces = [];
for i = 1:length(node_link)
    node1 = node_link(i, 1);
    node2 = node_link(i, 2);
    x1 = node_coords(node1, 1);
    y1 = node_coords(node1, 2);
    z1 = node_coords(node1, 3);
    x2 = node_coords(node2, 1);
    y2 = node_coords(node2, 2);
    z2 = node_coords(node2, 3);
    L = sqrt((x2 - x1).^2 + (y2 - y1).^2 + (z2 - z1).^2);
    k_local = get_local_stiffness(E, G, A, J, L, Iz, Iy, Ix);
    r = get_r(x1, y1, z1, x2, y2, z2);
    k_universal = r'*k_local*r;
    a = -5:0;
    dofs = [6*node1 + a, 6*node2 + a];
    force_on_beam = k_universal * displacements(dofs);
    m1 = sqrt(sum(force_on_beam(4:6).^2));
    m2 = sqrt(sum(force_on_beam(10:12).^2));
    
    forces = [forces, m1, m2];
end
[max_m, idx] = max(forces)
node_max = node_link(ceil(idx / 2), mod(idx, 2) + 1)



% After get the displacements, we could use the transformed universal stiffness matrix for a beam to calculate the 3 decomponets of M







function k = get_local_stiffness(E, G, A, J, L, Iz, Iy, Ix)
k1 = E*A/L;
k2 = 12*E*Iz/L.^3;
k3 = 6*E*Iz/L.^2;
k4 = 4*E*Iz/L;
k5 = 2*E*Iz/L;
k6 = 12*E*Iy/L.^3;
k7 = 6*E*Iy/L.^2;
k8 = 4*E*Iy/L;
k9 = 2*E*Iy/L;
k10 = G*J/L;

k = [k1 0 0 0 0 0 -k1 0 0 0 0 0;
    0 k2 0 0 0 k3 0 -k2 0 0 0 k3;
    0 0 k6 0 -k7 0 0 0 -k6 0 -k7 0;
    0 0 0 k10 0 0 0 0 0 -k10 0 0;
    0 0 -k7 0 k8 0 0 0 k7 0 k9 0;
    0 k3 0 0 0 k4 0 -k3 0 0 0 k5;
    -k1 0 0 0 0 0 k1 0 0 0 0 0;
    0 -k2 0 0 0 -k3 0 k2 0 0 0 -k3;
    0 0 -k6 0 k7 0 0 0 k6 0 k7 0;
    0 0 0 -k10 0 0 0 0 0 k10 0 0;
    0 0 -k7 0 k9 0 0 0 k7 0 k8 0;
    0 k3 0 0 0 k5 0 -k3 0 0 0 k4];
end
function lambda_3 = get_lambda(x1, y1, z1, x2, y2, z2)
% If x, y are equal, nan may appear in the following matrix
if x1 == x2 && y1 == y2
    % fprintf('vertical\n')
    if z2 > z1
        lambda_3 = [0, 0, 1; 0, 1, 0; -1, 0, 0];
    else
        lambda_3 = [0, 0, -1; 0, 1, 0; 1, 0, 0];
    end
else
    % fprintf('nonvertical\n')
    % We need to calculate the cosine beteen axis of local and universal
    L = sqrt((x2 - x1).^2 + (y2 - y1).^2 + (z2 - z1).^2);
    % Since x of local coor is along the truss
    % lambda could be decomposed of matrixs of rotation of 3 axis
    % Suppose the local z doesnt change
    
    CXx = (x2 - x1) / L;
    CYx = (y2 - y1) / L;
    CZx = (z2 - z1) / L;
    rho = sqrt(CXx.^2 + CYx.^2);
    CXy = -CYx /rho;
    CYy = CXx / rho;
    CZy = 0;
    CXz = -CXx*CZx / rho;
    CYz = -CYx*CZx / rho;
    CZz = rho;
    lambda_3 = [CXx, CYx, CZx; CXy, CYy, CZy;CXz, CYz, CZz];
end
end

function r = get_r(x1, y1, z1, x2, y2, z2)


r = zeros([12, 12]);
r(1: 3, 1: 3) = get_lambda(x1, y1, z1, x2, y2, z2);
r(4:6, 4:6) = get_lambda(x1, y1, z1, x2, y2, z2);
r(7:9, 7:9) = get_lambda(x1, y1, z1, x2, y2, z2);
r(10:12, 10:12) = get_lambda(x1, y1, z1, x2, y2, z2);
end

function displacements = solution(GDof,boundary_dof,stiffness,force)

activeDof = setdiff((1:GDof)', boundary_dof); % Make the matrix full rank, remove the bounded dofs
stiffness;
U = stiffness(activeDof,activeDof)\force(activeDof);
displacements = zeros(GDof,1);
displacements(activeDof) = U;
end




