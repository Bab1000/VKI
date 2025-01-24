import numpy as np
import matplotlib.pyplot as plt
import os

def FDM_coeff(stencil, derivative_order):

    # number of points in the stencil
    n = len(stencil)

    # create matrix A
    A = np.zeros([n,n])

    for col in range(n):
        z = 0
        for row in range (n):
            # general formula for each coefficient based on taylor expansion
            A[row][col] = stencil[col]**z/np.math.factorial(z)
            z= z+1

    # create matrix b
    b = np.zeros(n)
    b[derivative_order] = 1

    # compute linear system to obtain the coefficients 
    c = np.linalg.inv(A)@b
    
    return c

def Mesh_generation(L,n_points):
    mesh = np.linspace(0,L,n_points)
    return mesh

def Source_term(k,mesh):
    source = k * 4 * np.pi**2 * np.sin(2 * np.pi * mesh)
    return source

def Source_term_variable_k(k0,mesh):
    source = - 2*np.pi*np.cos(2*np.pi*mesh)*2*np.pi*np.cos(2*np.pi*mesh) + (225 + np.sin(2*np.pi*mesh))*4*np.pi**2*np.sin(2*np.pi*mesh)
    return source

def GP_manager_nonperiodicBC(u0,idx,coeff_Tx):
     # Pour idx-3
    if idx >= 3:
        u0minus3 = u0[idx-3]
    else:
        u0minus3 = 800

    # Pour idx-2
    if idx >= 2:
        u0minus2 = u0[idx-2]
    else:
        u0minus2 = 800

    # Pour idx-1
    if idx >= 1:
        u0minus1 = u0[idx-1]
    else:
        u0minus1 = 800

    # Pour idx+1
    if idx < len(u0)-1:
        u0plus1 = u0[idx+1]
    else:
        #compute the RHS ghost point value
        stencil_GP = [u0[idx - 3], u0[idx - 2], u0[idx - 1], u0[idx - 0]]

        # compute the ghost point
        RHS_ghost_point_value = -sum(coeff_Tx[0:(len(coeff_Tx)-1)] * stencil_GP)
        RHS_ghost_point_value = RHS_ghost_point_value/coeff_Tx[-1]

        # updating the RHS ghost point value
        u0plus1 = RHS_ghost_point_value

    return u0minus3, u0minus2, u0minus1, u0plus1

def correction_periodicBC(Txx,u0,coeff):
    # Periodic boundary condition
    Txx[0] = (
    coeff[0] * u0[-4] +
    coeff[1] * u0[-3] +
    coeff[2] * u0[-2] +
    coeff[3] * u0[0] +
    coeff[4] * u0[1] 
    )
    
    Txx[1] = (
    coeff[0] * u0[-3] +
    coeff[1] * u0[-2] +
    coeff[2] * u0[0] +
    coeff[3] * u0[1] +
    coeff[4] * u0[2] 
    )

    Txx[2] = (
    coeff[0] * u0[-2] +
    coeff[1] * u0[0] +
    coeff[2] * u0[1] +
    coeff[3] * u0[2] +
    coeff[4] * u0[3] 
    )

    Txx[-1] = (
    coeff[0] * u0[-4] +
    coeff[1] * u0[-3] +
    coeff[2] * u0[-2] +
    coeff[3] * u0[-1] +
    coeff[4] * u0[1] 
    )
    
    return Txx

def plot_mesh_convergence(err,test_mesh,convergence_slope,res):

    print("Plotting mesh convergence results ...")

    os.makedirs("Results/Mesh_convergence", exist_ok=True)

    save_path_1 = os.path.join("Results", "Mesh_convergence", "mesh_convergence.pdf")
    plt.figure()
    y1 = np.array(test_mesh, dtype=float)**(-3)
    plt.loglog(test_mesh, y1, '--', label='Slope : 3')
    y2 = np.array(test_mesh, dtype=float)**(-2)
    plt.loglog(test_mesh, y2, '--', label='Slope : 2')
    y3 = np.array(test_mesh, dtype=float)**(-1)
    plt.loglog(test_mesh, y3, '--', label='Slope : 1')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')  
    plt.rc('axes', unicode_minus=False)
    plt.loglog(test_mesh, err, 'o-', label=f'Scheme slope : {convergence_slope:.2f}')
    plt.xlabel('Mesh size', fontsize=12)
    plt.ylabel('Error', fontsize=12)
    plt.legend(fontsize=12)
    label = [f"{mesh}" for mesh in test_mesh]
    plt.xticks(test_mesh, labels=label, fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    plt.savefig(save_path_1, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()

    save_path_2 = os.path.join("Results", "Mesh_convergence", "results.pdf")
    mesh = Mesh_generation(1,test_mesh[-1])
    ref = np.sin(2 * np.pi * mesh)
    plt.figure()
    plt.plot(mesh, ref, linestyle='-', color='b', linewidth=1.5, label="$T_{ref}(x,t)$")
    plt.plot(mesh, res, linestyle='--', color='r', linewidth=1.5, label="T(x,t) (num.)")
    plt.xlabel('x [m]', fontsize=12)
    plt.ylabel('Temperature [K]', fontsize=12)
    plt.legend(fontsize=12)
    plt.grid(True)
    plt.savefig(save_path_2, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"   ---> Graph saved in : {save_path_1} | {save_path_2}")

def plot_mesh_convergence_variable_k(err, test_mesh, convergence_slope,residuals):
    print("Plotting mesh convergence results for variable k...")
    
    os.makedirs("Results/Mesh_convergence_variable_k", exist_ok=True)

    save_path = os.path.join("Results", "Mesh_convergence_variable_k", "mesh_convergence_variable_k.pdf")

    plt.figure()
    y1 = np.array(test_mesh, dtype=float)**(-3)
    plt.loglog(test_mesh, y1, '--', label='Slope : 3')
    y2 = np.array(test_mesh, dtype=float)**(-2)
    plt.loglog(test_mesh, y2, '--', label='Slope : 2')
    y3 = np.array(test_mesh, dtype=float)**(-1)
    plt.loglog(test_mesh, y3, '--', label='Slope : 1')

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')  
    plt.rc('axes', unicode_minus=False)

    plt.loglog(test_mesh, err, 'o-', label=f'Scheme slope : {convergence_slope:.2f}')
    plt.xlabel('Mesh size', fontsize=12)
    plt.ylabel('Error', fontsize=12)
    plt.legend(fontsize=12)
    
    label = [f"{mesh}" for mesh in test_mesh]
    plt.xticks(test_mesh, labels=label, fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)

    plt.savefig(save_path, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"   ---> Graph saved in : {save_path}")

    print("Plotting residuals for variable k simulation ...")
    
    os.makedirs("Results/Residuals_variable_k", exist_ok=True)

    save_path = os.path.join("Results", "Residuals_variable_k", "Residuals_variable_k.pdf")

    nt = np.linspace(0,len(residuals[0][:])-1,len(residuals[0][:]))
    for i in range(len(test_mesh)):
        plt.loglog(nt, residuals[i][:], linestyle='-', linewidth=1.5, label=f"Residuals for mesh {test_mesh[i]} ")
    plt.xlabel('Iterations', fontsize=12)
    plt.ylabel('Residual', fontsize=12)
    plt.legend()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    plt.savefig(save_path, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"   ---> Graph saved in : {save_path}")

def plot_stability(residual):
    print("Plotting stability results...")

    os.makedirs("Results/Stability", exist_ok=True)
    save_path = os.path.join("Results", "Stability", "stability.pdf")

    plt.figure()
    plt.loglog(residual, linestyle='-', color='b', linewidth=1.5, label='Residual Convergence')
    plt.xlabel('Iterations', fontsize=12)
    plt.ylabel('Residual', fontsize=12)
    plt.legend(fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.savefig(save_path, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"   ---> Graph saved in: {save_path}")

def plot_nonperiodicBC(res, mesh,dt):
    print("Plotting non-periodic boundary conditions results...")

    # Create the directory to save the results
    os.makedirs("Results/NonperiodicBC", exist_ok=True)
    save_path = os.path.join("Results", "NonperiodicBC", "nonperiodicBC.pdf")

    # Select 10 equally spaced indices
    res = np.array(res)
    indices = np.linspace(0, res.shape[0] - 1, 7, dtype=int)

    # Plot the 10 selected lines
    plt.figure()
    for i in indices:
        plt.plot(mesh, res[i][:], linestyle='-', linewidth=1.5, label=f'Time : {i*dt:.4}')

    # Add labels and customize the display
    plt.xlabel('Position [m]', fontsize=12)
    plt.ylabel('Temperature [k]', fontsize=12)
    plt.legend()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)

    # Save the plot in PDF format
    plt.savefig(save_path, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"   ---> Graph saved in: {save_path}")
