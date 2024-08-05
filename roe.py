# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from numpy import *
import pandas as pd
from tqdm import *


def compute_flux(q, gamma, seps):
    intervals = seps - 1
    dx = 1 / intervals
    seps = intervals + 1
    x = np.linspace(0 + dx / 2., 1, seps)
    areas = 4.2 - 2 * np.sqrt(4 - (x - 0.5) ** 2)

    rho = q[0] / areas
    u = q[1] / rho / areas
    E = q[2] / rho / areas
    p = (gamma - 1.) * rho * (E - 0.5 * u ** 2)

    # Flux vector
    F0 = np.array(rho * u) * areas
    F1 = np.array(rho * u ** 2 + p) * areas
    F2 = np.array(u * (rho * E + p)) * areas
    flux = np.array([F0, F1, F2])

    return flux


def compute_roe_flux(q, dx, gamma, a, nx):
    rho = q[0]
    u = q[1] / rho
    E = q[2] / rho
    p = (gamma - 1.) * rho * (E - 0.5 * u ** 2)
    htot = gamma / (gamma - 1) * p / rho + 0.5 * u ** 2

    Phi = np.zeros((3, nx - 1))

    for j in range(0, nx - 1):
        R = sqrt(rho[j + 1] / rho[j])
        rmoy = R * rho[j]
        u_hat = (R * u[j + 1] + u[j]) / (R + 1)
        h_hat = (R * htot[j + 1] + htot[j]) / (R + 1)
        a_hat = sqrt((gamma - 1.0) * (h_hat - 0.5 * u_hat * u_hat))

        alpha1 = (gamma - 1) * u_hat * u_hat / (2 * a_hat * a_hat)
        alpha2 = (gamma - 1) / (a_hat * a_hat)

        wdif = q[:, j + 1] - q[:, j]

        P = np.array([[1, 1, 1],
                      [u_hat - a_hat, u_hat, u_hat + a_hat],
                      [h_hat - a_hat * u_hat, 0.5 * u_hat * u_hat, h_hat + a_hat * u_hat]])

        Pinv = np.array([[0.5 * (alpha1 + u_hat / a_hat), -0.5 * (alpha2 * u_hat + 1 / a_hat), alpha2 / 2],
                         [1 - alpha1, alpha2 *
                          u_hat, -alpha2],
                         [0.5 * (alpha1 - u_hat / a_hat), -0.5 * (alpha2 * u_hat - 1 / a_hat), alpha2 / 2]])

        # Compute matrix Lambda_{j+1/2}
        lamb = np.array([[abs(u_hat - a_hat), 0, 0],
                         [0, abs(u_hat), 0],
                         [0, 0, abs(u_hat + a_hat)]])

        # Compute Roe matrix |A_{j+1/2}|
        A = np.dot(P, lamb)
        A = np.dot(A, Pinv)

        # Compute |A_{j+1/2}| (W_{j+1}-W_j)
        Phi[:, j] = np.dot(A, wdif)

    # Compute Phi=(F(W_{j+1}+F(W_j))/2-|A_{j+1/2}| (W_{j+1}-W_j)/2
    F = compute_flux(q, gamma, nx)
    Phi = 0.5 * (F[:, 0:nx - 1] + F[:, 1:nx]) - 0.5 * Phi

    dflux = (Phi[:, 1:-1] - Phi[:, 0:-2])

    return dflux


def theory_plot():
    df = pd.read_excel(
        r'D:\Python Projects\24Spring\Computation Mechanics\ROE\参考解.xlsx', engine='openpyxl')

    x, rho, u, p = df.to_numpy().T
    plt.subplot(3, 1, 1)
    plt.plot(x, rho, label='theory', alpha=0.5)
    plt.ylabel('$\\rho$', fontsize=16)
    plt.tick_params(axis='x', bottom=False, labelbottom=False)
    plt.grid(True)

    plt.subplot(3, 1, 2)
    plt.plot(x, u, label='theory', alpha=0.5)
    plt.ylabel('$U$', fontsize=16)
    plt.tick_params(axis='x', bottom=False, labelbottom=False)
    plt.grid(True)

    plt.subplot(3, 1, 3)
    plt.plot(x, p, label='theory', alpha=0.5)
    plt.ylabel('$p$', fontsize=16)
    plt.tick_params(axis='x', bottom=False, labelbottom=False)
    plt.grid(True)

    plt.subplots_adjust(left=0.2)
    plt.subplots_adjust(bottom=0.15)
    plt.subplots_adjust(top=0.95)


def interpolator(x0, y0, x1=None):
    if x1 is None:
        x1 = np.linspace(1 / 2 * 1 / 1600, 1, 1601)
    y1 = np.interp(x1, x0, y0)
    return x1, y1


def solver(ncells, show=False, tEnd=0.18, rho_only=False):
    time_step = 0.50
    gamma = 1.4
    tEnd = tEnd
    ncells = ncells
    dx = 1 / ncells
    nx = ncells + 1
    x = np.linspace(0 + dx / 2., 1, nx)  # Mesh

    rho0 = np.zeros(nx)
    u0 = np.zeros(nx)
    p0 = np.zeros(nx)
    mid = int(ncells / 2)
    areas = 4.2 - 2 * np.sqrt(4 - (x - 0.5) ** 2)

    areas_x = 2 * ((x - 0.5) / np.sqrt(4 - (x - 0.5) ** 2))

    p0[:mid] = 1.0
    p0[mid:] = 0.1
    u0[:mid] = 0.0
    u0[mid:] = 0.0
    rho0[:mid] = 1.0
    rho0[mid:] = 0.125

    E0 = p0 / ((gamma - 1.) * rho0) + 0.5 * u0 ** 2  # Total Energy density
    a0 = sqrt(gamma * p0 / rho0)  # Speed of sound
    # Vector of conserved variables
    q = np.array([rho0 * areas, rho0 * u0 * areas, rho0 * E0 * areas])

    t = 0
    it = 0
    a = a0
    dt = time_step * dx / max(abs(u0) + a0)

    p = p0
    while True:
        extra = np.zeros([3, nx])
        extra[1, :] = p * areas_x

        q0 = q.copy()
        dF = compute_roe_flux(q0, dx, gamma, a, nx)

        q[:, 1:-2] = q0[:, 1:-2] - dt / dx * dF + dt * extra[:, 1:-2]

        q[:, 0] = q0[:, 0]
        q[:, -1] = q0[:, -1]  # Dirichlet

        rho = q[0] / areas
        u = q[1] / (rho * areas)
        E = q[2] / (rho * areas)
        p = (gamma - 1.) * rho * (E - 0.5 * u ** 2)
        a = sqrt(gamma * p / rho)

        dt = time_step * dx / max(abs(u) + a)
        delta_t = tEnd - t
        if abs(delta_t) <= abs(dt):
            dt = delta_t

        t = t + dt
        it = it + 1

        if t >= (tEnd - 1e-5):
            print(it, t)
            if rho_only:
                return x, rho

            plt.subplot(3, 1, 1)
            plt.title(f'T = {tEnd}')
            plt.plot(x, rho, label=f'mesh={ncells}', alpha=0.5)
            plt.ylabel('$\\rho$', fontsize=16)
            plt.tick_params(axis='x', bottom=False, labelbottom=False)
            plt.grid(True)

            plt.subplot(3, 1, 2)
            plt.plot(x, u, label=f'mesh={ncells}', alpha=0.5)
            plt.ylabel('$U$', fontsize=16)
            plt.tick_params(axis='x', bottom=False, labelbottom=False)
            plt.grid(True)

            plt.subplot(3, 1, 3)
            plt.plot(x, p, label=f'mesh={ncells}', alpha=0.5)
            plt.ylabel('$p$', fontsize=16)
            plt.tick_params(axis='x', bottom=False, labelbottom=False)
            plt.grid(True)

            plt.subplots_adjust(left=0.2)
            plt.subplots_adjust(bottom=0.15)
            plt.subplots_adjust(top=0.95)
            if show:
                plt.show()

            break


fig, axes = plt.subplots(nrows=3, ncols=1)
tEnd = 0.18
solver(200, tEnd=tEnd)
solver(400, tEnd=tEnd)
solver(800, tEnd=tEnd)
theory_plot()
plt.legend()
plt.show()

x_std, y_std = solver(300, rho_only=True)

eps = []
# mesh_lst = [10, 20, 50, 100, 200, 400, 600, 800, 1200]
mesh_lst = [10, 20, 50, ]
interval_lst = [1 / i for i in mesh_lst]
for (mesh, interval) in tqdm(zip(mesh_lst, interval_lst)):
    x_i, y_i = interpolator(*solver(mesh, rho_only=True), x_std)
    e = np.sqrt(interval * np.sum((y_std - y_i) ** 2))
    eps.append(e)
plt.title('$\\epsilon_{\\rho} \sim interval$')
plt.xlabel('Interval')
plt.ylabel('Eps')

plt.loglog(interval_lst, eps, '-*')
plt.show()
