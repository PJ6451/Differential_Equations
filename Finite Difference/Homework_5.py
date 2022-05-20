import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import scipy.linalg as scilin
from numba import vectorize
import sympy as sym
import pandas as pd
import os
import imageio

def u0(x, y):
    cosh = lambda ex: .5*(np.exp(ex) + np.exp(-ex))
    return np.sin(1.2*(x - y))*cosh(x + 2*y)

def exact_sol(t, x, y):
    return np.exp(1.68*t)*u0(x,y)

def LF(x_bounds, t_max, u_init, l_bound, r_bound, h=1/20, lamb=.8, a=1):
    xl, xr = x_bounds
    x_points = np.arange(xl, xr+h, h)

    k = lamb*h
    t_points = np.arange(0, t_max+k, k)

    u_grid = np.tile(u_init(x_points), (t_points.size, 1))

    u_grid[:, 0] = l_bound(t_points)
    u_grid[:, -1] = r_bound(t_points)

    for ii in range(1, t_points.size):
        n = ii - 1

        u_grid[ii,1:-1] = u_grid[n,1:-1] - .5*a*lamb*(u_grid[n,2:]-u_grid[n,:-2]) + \
                          .5*a*a*lamb*lamb*(u_grid[n,2:]-2*u_grid[n,1:-1]+u_grid[n,:-2])


    return x_points, t_points, u_grid

def FTBS(x_bounds, t_max, u_init, l_bound, r_bound, h=1/2, lamb=.8, a=1):
    xl, xr = x_bounds
    x_points = np.arange(xl, xr+h, h)

    k = lamb*h
    t_points = np.arange(0, t_max+k, k)

    u_grid = np.tile(u_init(x_points), (t_points.size, 1))

    u_grid[:, 0] = l_bound(t_points)
    u_grid[:, -1] = r_bound(t_points)

    for ii in range(1, t_points.size):
        n = ii - 1

        u_grid[ii,1:-1] = a*lamb*u_grid[n,:-2] + (1 - a*lamb)*u_grid[n,1:-1]


    return x_points, t_points, u_grid

def box_scheme(x_bounds, t_max, u_init, boundary, forcing, h, lamb):
    xl, xr = x_bounds
    x_points = np.arange(xl, xr+h, h)

    k = lamb*h
    t_points = np.arange(0, t_max+k, k)

    u_grid = np.tile(u_init(x_points), (t_points.size, 1))
    u_grid[:, 0] = boundary(t_points)

    T, X = np.meshgrid(t_points, t_points)
    f_grid = forcing(T, X)

    for ii in range(1, t_points.size):
        for jj in range(1, x_points.size):
            #P = .125*np.sum(f_grid[ii-1:ii+1,jj-1:jj+1])*k
            P = (forcing(t_points[ii], x_points[jj]) + forcing(t_points[ii], x_points[jj-1]) +
                      forcing(t_points[ii-1], x_points[jj]) + forcing(t_points[ii-1], x_points[jj-1]))
            
            # u_grid[ii, jj] = (P+(lamb-1)*u_grid[ii,jj-1]+(1+lamb)*u_grid[ii-1,jj-1]+(1-lamb)*u_grid[ii-1,jj])/(1+lamb)

            u_grid[ii, jj] = (2*(k-h)*u_grid[ii,jj-1] + 2*(k+h)*u_grid[ii-1, jj-1] +2*(h-k)*u_grid[ii-1,jj] + P*k*h)/(2*(k+h))            

    return x_points, t_points, u_grid

def testing(n):
    x_points = np.linspace(0,1,n)
    y_points = np.linspace(0,1,n)

    u = lambda ex, why: np.exp(ex*why)

    X, Y = np.meshgrid(x_points, y_points)

    U = u(X, Y)

    main_diag = -2*np.ones(x_points.size)
    sub_diag = np.ones(x_points.size-1)

    A = np.diag(main_diag) + np.diag(sub_diag, k=-1) + np.diag(sub_diag, k=1)
    A[0, :2] = [1,0]
    A[-1, -2:] = [0,1]
    I = np.eye(A.shape[0])


    res1 = np.matmul(A, U)
    res2 = np.matmul(np.kron(I, A), U.T.flatten())
    res2 = np.reshape(res2, A.shape, order='F')

    print(res1)
    print(res2)

    res3 = np.matmul(A, U.T).T
    res4 = np.matmul(np.kron(A, I), U.T.flatten())
    res4 = np.reshape(res4, A.shape, order='F')

    print(res3)
    print(res4)
    
def peaceman_rachford(x_bounds, y_bounds, h1, h2, k, t_max, exact, a1=2, a2=1):
    xl, xr = x_bounds; yl, yr = y_bounds

    x_points = np.arange(xl, xr+h1, h1)
    y_points = np.arange(yl, yr+h2, h2)
    t_points = np.arange(0, t_max+k, k)

    mu_x = k/h1**2
    mu_y = k/h2**2

    halfbmux = .5*a1*mu_x
    halfbmuy = .5*a2*mu_y

    X, Y, T = np.meshgrid(x_points, y_points, t_points)

    U_exact = exact(T, X, Y).T
    U = np.transpose(np.copy(U_exact), axes=[0,2,1])

    for ii in range(1, t_points.size):
        v = np.copy(U[ii-1])
        w = np.copy(v)
        
        for m in range(1, y_points.size-1):
            p = [0]
            q = [U[ii,0,m]]

            for l in range(1, x_points.size-1):
                dd = v[l,m] + halfbmuy*(v[l, m-1] -2*v[l,m] + v[l, m+1])
                denom = 1 + halfbmux*(2 - p[-1])

                p.append(halfbmux/denom)
                q.append((dd + halfbmux*q[-1])/denom)

            w[-1, m] = U[ii,-1,m]

            for l in range(x_points.size-1)[::-1]:
                w[l, m] = p[l]*w[l+1,m]+q[l]


        for l in range(1, x_points.size-1):
            p = [0]
            q = [U[ii,l,0]]

            for m in range(1, y_points.size-1):
                dd = w[l,m] + halfbmux*(w[l-1,m]-2*w[l,m]+w[l+1,m])
                denom = 1 + halfbmuy*(2 - p[-1])
                p.append(halfbmuy/denom)
                q.append((dd + halfbmuy*q[-1])/denom)
            
            v[l,-1] = U[ii,l,-1]

            for m in range(y_points.size-1)[::-1]:
                v[l, m] = p[m]*v[l,m+1]+q[m]
            
        v[0,:] = U[ii,0,:]
        v[-1,:] = U[ii,-1,:]

        U[ii] = v

    return t_points, x_points, y_points, np.transpose(U,axes=[0,2,1]), U_exact


P737 = 1
P312 = 0

# P 7.3.7
# ===========================================================================================
if P737:
    x_bounds = (0,1)
    y_bounds = (0,1)

    t_max = 1.2

    dxs = [.1, .05, .025]

    l2_errors = np.zeros(len(dxs))
    inf_errors = np.zeros(len(dxs))

    for jj, dx in enumerate(dxs):

 
        t, x, y, result, exact = peaceman_rachford(x_bounds, y_bounds, dx, dx, dx, t_max, exact_sol)
        mid = t.size//4

        l2_errors[jj] = np.sqrt(dx)*np.linalg.norm(result[mid] - exact[mid], ord=2)
        # inf_errors[jj] = np.linalg.norm(result[mid] - exact[mid], ord=np.inf) 
        inf_errors[jj] = np.max(np.abs(result[mid] - exact[mid])) 

        X, Y = np.meshgrid(x,y)
        images = [None]*t.size
        
        for ii, (e_sol, n_sol) in enumerate(zip(exact, result)):

            fig = plt.figure(figsize=(16,8))

            ax1 = fig.add_subplot(1,2,1, projection='3d')
            ax1.plot_surface(X, Y, n_sol)
            ax1.set_xlabel("$x$", size=15)
            ax1.set_ylabel("$y$", size=15)
            ax1.set_zlabel("$u(t,x,y)$", size=15)
            ax1.set_zlim(result.min(), result.max())
            ax1.tick_params(axis='both', labelsize=10)
            ax1.set_title("Peaceman-Rachford", size=20)

            ax2 = fig.add_subplot(1,2,2, projection='3d')
            ax2.plot_surface(X, Y, e_sol)
            ax2.set_xlabel("$x$", size=15)
            ax2.set_ylabel("$y$", size=15)
            ax2.set_zlabel("$u(t,x,y)$", size=15)
            ax2.set_zlim(exact.min(), exact.max())
            ax2.tick_params(axis='both', labelsize=10)
            ax2.set_title("Exact", size=20)

            fig.suptitle(f"$dx = dy = dt = {dx},\\ t = {t[ii]:.3f}$", size=25)

            fig.tight_layout()
            fig.savefig(f"HW5_P737_Peaceman_Rachford_{jj:03d}_{ii:03d}.png")        
            plt.close(fig)

            images[ii] = imageio.imread(f"HW5_P737_Peaceman_Rachford_{jj:03d}_{ii:03d}.png")
            os.remove(f"HW5_P737_Peaceman_Rachford_{jj:03d}_{ii:03d}.png")

        imageio.mimsave(f"HW5_P737_Peaceman_Rachford_GIF_dx={dx:.4f}.gif", images, duration=5/t.size)


    info_dict = {"dx": dxs}

    error_measurements = [l2_errors, inf_errors]
    error_labels = ["L2", "Infty"]

    for error_vec, e_label in zip(error_measurements, error_labels):

        info_dict[e_label+" Error"] = error_vec
        RoCs = np.zeros(error_vec.size)
        RoCs[1:] = np.log(error_vec[:-1]/error_vec[1:])/np.log(2)

        info_dict[e_label+" OoA"] = RoCs


    data = pd.DataFrame(info_dict).set_index('dx')

    print("Peaceman-Rachford")
    print(data.round(7),"\n")

# P 3.1.2
# ===========================================================================================

if P312:
    u_exact = lambda tea, ex: np.sin(2*np.pi*(ex - tea))
    u_init = lambda ex: u_exact(0, ex)

    x_bounds = (-1,1)

    l_bc = lambda tea: u_exact(tea, x_bounds[0])
    r_bc = lambda tea: u_exact(tea, x_bounds[-1])

    t_max = 1.2

    hvals = [.1, .05, .025, .0125]

    methods = [FTBS, LF]
    method_labels = ["Forward-time Backward-space", "Lax-Friedrichs"]

    for method, meth_lab in zip(methods,method_labels):

        l2_errors = np.zeros(len(hvals))
        inf_errors = np.zeros(len(hvals))

        for ii, h in enumerate(hvals):

            x, t, u = method(x_bounds, t_max, u_init, l_bc, r_bc, h)

            mid = t.size//2

            X, T = np.meshgrid(x, t)
            U_exact = u_exact(T, X)

            images = [None]*t.size

            for jj, (e_sol, n_sol) in enumerate(zip(U_exact, u)):


                # fig = plt.figure(figsize=(16,8))
                fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,7))


                # ax1 = fig.add_subplot(1,2,1, projection='3d')
                # ax2 = fig.add_subplot(1,2,2, projection='3d')


                # ax1.set_prop_cycle('color', plt.cm.jet(steps))
                # ax1.plot_surface(X, T, u, cmap='jet', edgecolor='none')

                ax1.plot(x, n_sol, '-k')
                ax1.set_xlabel("$x$", size=15)
                # ax1.set_ylabel("$t$", size=15)
                ax1.set_ylabel("$u(t,x)$", size=15)
                ax1.set_ylim(U_exact.min(), U_exact.max())
                ax1.set_title(meth_lab, size=20)
                ax1.grid(True, which='both')
                ax1.tick_params(axis='both', labelsize=10)
                # ax1.view_init(azim=-115, elev=30)


                # ax2.set_prop_cycle('color', plt.cm.jet(steps))
                # ax2.plot_surface(X, T, U_exact, cmap='jet')

                ax2.plot(x, e_sol, '-k')
                ax2.set_xlabel("$x$", size=15)
                # ax2.set_ylabel("$t$", size=15)
                ax2.set_ylabel("$u(t,x)$", size=15)
                ax2.set_ylim(U_exact.min(), U_exact.max())
                ax2.set_title(f"Exact", size=20)
                ax2.grid(True, which='both')
                ax2.tick_params(axis='both', labelsize=10)
                # ax2.view_init(azim=-115, elev=30)

                fig.suptitle(f"$\\lambda = 1.2,\\ h = {h},\\ k = {.8*h:.3f},\\ t = {t[jj]:.3f}$", size=25)

                fig.tight_layout()


                fig.savefig(f"HW5_P312_{meth_lab.replace(' ','-').replace('-','_')}_{ii:03d}_{jj:03d}.png")
                plt.close(fig)

                images[jj] = imageio.imread(f"HW5_P312_{meth_lab.replace(' ','-').replace('-','_')}_{ii:03d}_{jj:03d}.png")
                os.remove(f"HW5_P312_{meth_lab.replace(' ','-').replace('-','_')}_{ii:03d}_{jj:03d}.png")

            imageio.mimsave(f"HW5_P312_{meth_lab}_GIF_h={h:.4f}.gif", images, duration=5/t.size)


            l2_errors[ii] = np.sqrt(h)*np.linalg.norm(u[mid] - U_exact[mid], ord=2)
            inf_errors[ii] = np.linalg.norm(u[mid] - U_exact[mid], ord=np.inf)

        info_dict = {"h": hvals}

        error_measurements = [l2_errors, inf_errors]
        error_labels = ["L2", "Infty"]

        for error_vec, e_label in zip(error_measurements, error_labels):

            info_dict[e_label+" Error"] = error_vec
            RoCs = np.zeros(error_vec.size)
            RoCs[1:] = np.log(error_vec[:-1]/error_vec[1:])/np.log(2)

            info_dict[e_label+" OoA"] = RoCs


        data = pd.DataFrame(info_dict).set_index('h')

        print(meth_lab)
        print(data.round(7),"\n")

