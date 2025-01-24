import numpy as np
import matplotlib.pyplot as plt
from Utils import *
import pdb

def Mesh_convergence(test_mesh,stencil,coeff,k,beta,nt):

    err = []
    res_vec = []

    for i in range(len(test_mesh)):

        # generating the mesh currently being analyzed
        mesh = Mesh_generation(1,test_mesh[i])

        # computing global inputs
        dx = mesh[1] - mesh[0]
        dt = 2*beta*dx**2/k

        # computing initial condition
        u0 = np.sin(2*np.pi*mesh)

        # computing the result of the mesh analysis
        res,residual = Temperature_diffusion_with_source_periodicBC(u0,stencil,mesh,coeff,k,dx,dt,nt)

        res = res[-1]

        res_vec.append(res)

        # error computation
        err.append(np.linalg.norm(res - u0 ,ord=2) * dx)

    return err,res_vec

def Mesh_convergence_variable_k(test_mesh,stencil,coeff_Txx,coeff_Tx,coeff_k,k,beta,nt):

    err = []

    residual_vec = []

    for i in range(len(test_mesh)):

        # generating the mesh currently being analyzed
        mesh = Mesh_generation(1,test_mesh[i])

        # computing global inputs
        dx = mesh[1] - mesh[0]
        dt = 2*beta*dx**2/k

        # computing initial condition
        u0 = np.sin(2*np.pi*mesh)
        k0 = k+np.sin(2*np.pi*mesh)

        # computing the result of the mesh analysis
        res,residual = Temperature_diffusion_with_source_periodicBC_variableK(u0,k0,stencil,mesh,coeff_Txx,coeff_Tx,coeff_k,dx,dt,nt)

        res = res[-1]

        # error computation
        err.append(np.linalg.norm(res - u0 ,ord=2))
        residual_vec.append(residual)

    return err,residual_vec

def Stability_analysis(test_beta,mesh,stencil,coeff,k,nt):

    # computing global inputs
    dx = mesh[1] - mesh[0]

    for i in range(len(test_beta)):

        print(f"Stability test for \033[34mbeta = {test_beta[i]:.2f}\033[0m")

        dt = 2*test_beta[i]*dx**2/k

        # computing initial condition
        u0 = np.sin(2*np.pi*mesh)

        # computing the result of the mesh analysis
        res,residual_vec = Temperature_diffusion_with_source_periodicBC(u0,stencil,mesh,coeff,k,dx,dt,nt)

        print(f"   ---> residual : {residual_vec[-1]}")

        if (residual_vec[-1] > 1):
            break

    return residual_vec

def Temperature_diffusion_with_source_periodicBC(u0,stencil,mesh,coeff,k,dx,dt,nt):

    # Prepare to store solutions at each time step
    u_res = [u0.copy()]

    # compute the source term (base on the result obtained by injecting u0 in the equation while assuming d/dt = 0)
    source = Source_term(k,mesh)

    residual_vec = []
    residual_0 = []

    for i in range(nt):
        
        # compute the spatial discretization
        Txx = sum(coeff[k] * np.roll(u0,-stencil[k]) for k in range(len(coeff)))

        Txx = correction_periodicBC(Txx,u0,coeff)

        Txx = k*Txx/dx**2
        RHS = Txx + source

        residual_0.append(np.linalg.norm(RHS))

        # compute the time discretization based on Heun's method (improved Euler method --> second order)
        u_bar = u0 + dt * RHS

        Txx_bar = sum(coeff[k] * np.roll(u_bar,-stencil[k]) for k in range(len(coeff)))

        Txx_bar = correction_periodicBC(Txx_bar,u_bar,coeff)    

        Txx_bar = k*Txx_bar/dx**2
        RHS_bar = Txx_bar + source

        # compute the new values for the next time step
        u_new = u0 + 0.5 * dt * (RHS + RHS_bar)

        # storing values
        u_res.append(u_new)

        # checking for convergence
        residual = np.linalg.norm(RHS)       
        residual = residual/residual_0[0]
        residual_vec.append(residual)

        if(residual < 1e-7):
            print(f"   ---> \033[34mMesh {len(mesh)}\033[0m : The solution has converged after {i} iterations")
            break
        elif(residual > 1):
            print(f"   ---> \033[34mMesh {len(mesh)}\033[0m : \033[31m/!\ The solution has diverged after {i} iterations\033[0m")
            break
        elif(i == nt):
            print(f"   ---> \033[34mMesh {len(mesh)}\033[0m : Residuals did not reach 1e-7 convergence tolerance after {i} iterations")

        # updating the initail conditions
        u0 = u_new

    return u_res,residual_vec

def Temperature_diffusion_nonperiodicBC(u0,stencil,mesh,coeff_Txx,coeff_Tx,k,dx,dt,nt,dirichlet_bc):

   
    # Prepare to store solutions at each time step
    u_res = [u0.copy()]

    # applying BC
    print("Applying BC ...")
    print(f"   ---> \033[34mLeft wall BC\033[0m  : Dirichlet with T = {dirichlet_bc}")
    print("   ---> \033[34mRight wall BC\033[0m : Neumann with dT/dx = 0")
    print("Simulation of the problem ...")
    for i in range(nt):

        Txx = np.zeros(len(u0))
        Txx_bar = np.zeros(len(u0))

        for j in range(len(u0)):
                
            u0minus3, u0minus2, u0minus1, u0plus1 = GP_manager_nonperiodicBC(u0,j,coeff_Tx)

            Txx[j] = coeff_Txx[0]*u0minus3 + coeff_Txx[1]*u0minus2 + coeff_Txx[2]*u0minus1 + coeff_Txx[3]*u0[j] + coeff_Txx[4]*u0plus1 

        Txx = k*Txx/dx**2
        RHS = Txx
       
         # compute the time discretization based on Heun's method (improved Euler method --> second order)
        u_bar = u0 + dt * Txx

        for j in range(3, len(u0) - 1):

            u0minus3, u0minus2, u0minus1, u0plus1 = GP_manager_nonperiodicBC(u_bar,j,coeff_Tx)

            Txx_bar[j] = coeff_Txx[0]*u0minus3 + coeff_Txx[1]*u0minus2 + coeff_Txx[2]*u0minus1 + coeff_Txx[3]*u0[j] + coeff_Txx[4]*u0plus1 

        Txx_bar = k*Txx_bar/dx**2
        RHS_bar = Txx_bar

        # compute the new values for the next time step
        u_new = u0 + 0.5 * dt * (RHS + RHS_bar)

        # storing values
        u_res.append(u_new)

        # updating the initail conditions
        u0 = u_new

    return u_res

def Temperature_diffusion_with_source_periodicBC_variableK(u0,k0,stencil,mesh,coeff_Txx,coeff_Tx,coeff_k,dx,dt,nt):

    # Prepare to store solutions at each time step
    u_res = [u0.copy()]

    residual_vec = []
    residual_0 = []

    # compute the source term (base on the result obtained by injecting u0 in the equation while assuming d/dt = 0)
    source_variable_k = Source_term_variable_k(k0,mesh)

    for i in range(nt):

        Txx = sum(coeff_Txx[k] * np.roll(u0,-stencil[k]) for k in range(len(coeff_Txx)))
        Tx = sum(coeff_Tx[k] * np.roll(u0,-stencil[k]) for k in range(len(coeff_Tx)))
        kx = sum(coeff_k[k] * np.roll(k0,-stencil[k]) for k in range(len(coeff_k)))

        Txx = correction_periodicBC(Txx,u0,coeff_Txx)
        Tx = correction_periodicBC(Tx,u0,coeff_Tx)
        kx = correction_periodicBC(kx,k0,coeff_k)

        Txx = k0*Txx/dx**2
        kx_Tx = kx*Tx/dx**2
        RHS = kx_Tx + Txx + source_variable_k

        if i == 0:
            residual_0.append(np.linalg.norm(RHS))
  
        # compute the time discretization based on Heun's method (improved Euler method --> second order)
        u_bar = u0 + dt * RHS

        Txx_bar = sum(coeff_Txx[k] * np.roll(u_bar,-stencil[k]) for k in range(len(coeff_Txx)))
        Tx_bar = sum(coeff_Tx[k] * np.roll(u_bar,-stencil[k]) for k in range(len(coeff_Tx)))

        Txx_bar = correction_periodicBC(Txx_bar,u_bar,coeff_Txx)
        Tx_bar = correction_periodicBC(Tx_bar,u_bar,coeff_Tx)

        Txx_bar = k0*Txx_bar/dx**2
        kx_Tx_bar = kx*Tx_bar/dx**2
        RHS_bar = Txx_bar + kx_Tx_bar + source_variable_k

        # compute the new values for the next time step
        u_new = u0 + 0.5 * dt * (RHS + RHS_bar)

        # storing values
        u_res.append(u_new)
        
       # checking for convergence
        residual = np.linalg.norm(RHS) 
        residual = residual/residual_0[0]
        residual_vec.append(residual)

        if(residual < 1e-7):
            print(f"   ---> \033[34mMesh {len(mesh)}\033[0m : The solution has converged after {i} iterations")
            break
        elif(residual > 2):
            print(f"   ---> \033[34mMesh {len(mesh)}\033[0m : \033[31m/!\ The solution has diverged after {i} iterations\033[0m")
            break
        elif(i == nt-1):
            print(f"   ---> \033[34mMesh {len(mesh)}\033[0m : Residuals did not reach 1e-7 convergence tolerance after {i} iterations")

        # updating the initail conditions
        u0 = u_new


    return u_res,residual_vec








