import time
import numpy as np
from numpy import linalg
from typing import Literal
import copy
import sympy
with open('Elastic/mesh_1.txt', 'r', encoding='utf-8') as f:
    info = f.readlines()

    for (i, line) in enumerate(info):
        if r'% Coordinates' in line:
            cor_idx_found = True
            cor_idx = i
        if r'% Elements (triangles)' in line:
            elem_idx_found = True
            elem_idx = i
    assert cor_idx_found and elem_idx_found


def position_to_edge(coords, xlow=-np.inf, xhigh=np.inf, ylow=-np.inf, yhigh=np.inf, type='node'):
    tol = 1e-3
    condition_x = (xlow-tol <= coords[:, 0]) & (coords[:, 0] <= xhigh+tol)
    condition_y = (ylow-tol <= coords[:, 1]) & (coords[:, 1] <= yhigh+tol)

    # 使用 np.where 找到满足条件的索引
    indices = np.where(condition_x & condition_y)[0]
    if type == 'node':
        return indices
    elif type == 'x':
        return indices * 2
    elif type == 'y':
        return indices * 2 + 1


def elems_contains_coord(coords, elems, node_coord: list | np.ndarray):
    # node_coord = coords[node].reshape(1, -1)
    if isinstance(node_coord, list):
        node_coord = np.array(node_coord)
    elem_container = []
    for (i, elem) in enumerate(elems):
        elem_coords = coords[elem]
        vector = elem_coords - np.tile(node_coord, [3, 1])
        in_tri = np.cross(vector[0], vector[1]) >= 0 \
            and np.cross(vector[1], vector[2]) >= 0\
            and np.cross(vector[2], vector[0]) >= 0
        if in_tri:
            elem_container.append(i)
    return elem_container


def to_float(x): return [float(i) for i in x]
def to_int(x): return [int(i) for i in x]


def find_third_point(elems, n1, n2):
    for elem in elems:
        temp = copy.deepcopy(elem)
        if n1 in elem and n2 in elem:
            return temp[~np.isin(temp, [n1, n2])]


def get_N(node_i, node_j, node_k):
    # 对顺序有要求
    x_i, y_i = coords[node_i]
    x_j, y_j = coords[node_j]
    x_k, y_k = coords[node_k]
    a_i = x_j*y_k - x_k * y_j
    b_i = y_j - y_k
    c_i = x_k - x_j
    Delta = linalg.det(
        np.hstack([np.ones([3, 1]), coords[[node_i, node_j, node_k]]])) / 2

    x = sympy.Symbol('x')
    y = sympy.Symbol('y')
    N = 1 / (2*Delta) * (a_i + b_i * x + c_i * y)
    return N


def get_Delta(x1, y1, x2, y2, x3, y3) -> np.ndarray:
    G = np.array([[1, x1, y1],
                  [1, x2, y2],
                  [1, x3, y3],])
    Delta = linalg.det(G) / 2
    return Delta


def get_B(x1, y1, x2, y2, x3, y3) -> np.ndarray:
    Delta = get_Delta(x1, y1, x2, y2, x3, y3)
    return 1 / (2 * Delta) * np.array([[y2-y3, 0, y3-y1, 0, y1-y2, 0],
                                       [0, x3-x2, 0, x1-x3, 0, x2-x1],
                                       [x3-x2, y2-y3, x1-x3, y3-y1, x2-x1, y1-y2]])


def node_to_dof(nodes: np.array, type):
    if isinstance(elem, list):
        nodes = np.array(nodes)
    if type == 'x':
        return 2 * nodes
    elif type == 'y':
        return 2 * nodes + 1
    elif type == 'all':
        dofs = np.zeros(2*len(nodes), dtype=int)
        dof_x = 2 * nodes
        dof_y = 2 * nodes + 1
        dofs[::2] = dof_x
        dofs[1::2] = dof_y
        return dofs


coords = np.array([to_float(i.split()) for i in info[cor_idx+1:elem_idx]])
elems = np.array([to_int(i.split()) for i in info[elem_idx+1:]]) - 1


# 平面应力
E = 200000000000
nu = 0.3
t = 1
# 平面应变
pmyb = 1

if pmyb:
    E = E / (1 - nu**2)
    nu = nu / (1-nu)
D = E / (1 - nu ** 2) * np.array([[1, nu, 0],
                                  [nu, 1, 0],
                                  [0, 0, (1 - nu)/2]])


total_dofs = len(coords) * 2

K_total = np.zeros([total_dofs, total_dofs])

displacements_total = np.zeros(total_dofs)
f_total = np.zeros(total_dofs, dtype=float)

# Get the global K
for node in elems:
    # nodes = elems[0] # 如 [0, 1, 2]
    # print(nodes)
    node_coord: np.ndarray = coords[node]

    G = np.hstack([np.ones([3, 1]), node_coord])
    # Delta = linalg.det(G) / 2

    B = get_B(*node_coord.flatten())
    Delta = get_Delta(*node_coord.flatten())
    K_local = Delta * t * B.T @ D @ B
    dof_x = 2 * node
    dof_y = 2 * node + 1

    dof_activated = np.zeros([6], dtype=int)
    dof_activated[::2] = dof_x
    dof_activated[1::2] = dof_y

    K_total[np.ix_(dof_activated, dof_activated)] += K_local


# print(displacements_total)
######################################
# 力
force_edge_node = position_to_edge(coords, 0, 4, 0, 0)
node_dct = dict(zip(force_edge_node, coords[force_edge_node, 0]))
sorted_node = sorted(node_dct, key=lambda k: node_dct[k])


x = sympy.Symbol('x')
y = sympy.Symbol('y')
fy = -500
for idx in range(len(sorted_node) - 1):
    node_i, node_j = sorted_node[idx:idx+2]
    node_k = find_third_point(elems, node_i, node_j)[0]
    a = float(sympy.integrate(get_N(node_i, node_j, node_k) * fy,
              (x, coords[node_i, 0], coords[node_j, 0])).subs(y, 0))
    b = float(sympy.integrate(get_N(node_j, node_i, node_k) * fy,
              (x, coords[node_i, 0], coords[node_j, 0])).subs(y, 0))

    f_total[[node_i*2 + 1, node_j*2 + 1]] += np.array([a, b])

node_jz_y = position_to_edge(coords, 2, 2, 0, 0, type='y')

f_total[node_jz_y] = -300

for idx in range(len(sorted_node) - 1):
    node_i, node_j = sorted_node[idx:idx+2]
    f_inte = float(sympy.integrate(
        fy, (x, coords[node_i, 0], coords[node_j, 0])).subs(y, 0))
    f_total[[node_i*2 + 1, node_j*2 + 1]] += np.array([f_inte/2, f_inte/2])


###################################
# 约束
deactivated_dof = []
deactivated_dof.extend(position_to_edge(coords, 0, 4, 2, 2, type='x'))
deactivated_dof.extend(position_to_edge(coords, 0, 4, 2, 2, type='y'))

activated_dof = np.arange(total_dofs)[~np.isin(
    np.arange(total_dofs), deactivated_dof)]


displacements_total[activated_dof] = (linalg.inv(
    K_total[np.ix_(activated_dof, activated_dof)])@f_total[activated_dof])

##################################################
n = position_to_edge(coords, 2, 2, 0, 0)
print(f'Displacement: \n{displacements_total[[n*2, n*2 + 1]]}')
print('*******************************')
####################
# 应变

stress_lst = []
for elem in elems[elems_contains_coord(coords, elems, [2, 0])]:
    elem_dof = node_to_dof(elem, 'all')
    # print(coords[elem])
    # print(elem)
    elem_displacement = displacements_total[elem_dof]
    # print(elem_displacement)
    B = get_B(*coords[elem].flatten())
    # print(B)
    strain_elem = get_B(*coords[elem].flatten()) @ elem_displacement
    stress = D @ strain_elem
    sx, sy, tau = stress
    print(stress)
    stress_lst.append(stress)
    ms = np.sqrt(0.25*(sx - sy)**2 + (tau) ** 2) + 0.5 * (sx + sy)
    # print(ms)
    # break
mean_stress = np.mean(np.array(stress_lst).reshape(-1, 3), axis=0)
print(f'mean stress: \n{mean_stress}')
print('********************************')
# quit()
main_stress = []
for (i, elem) in enumerate(elems):

    elem_dof = node_to_dof(elem, 'all')
    elem_displacement = displacements_total[elem_dof]
    Delta = get_Delta(*coords[elem].flatten())

    strain_elem = get_B(*coords[elem].flatten()) @ elem_displacement
    stress = D @ strain_elem
    sx, sy, tau = stress
    # ms = np.sqrt(0.25*(sx - sy)**2 + (tau) ** 2) + 0.5 * (sx + sy)
    ms = np.sqrt(0.5*((sx - sy)**2 + sx**2 + sy ** 2) + 3*tau**2)
    main_stress.append(
        ms
    )

idx = np.argmax(main_stress)
print(f'max stress: {max(main_stress)} at {idx} elem')
weight_coord = np.mean(coords[elems[idx]], axis=0)
print(weight_coord)
