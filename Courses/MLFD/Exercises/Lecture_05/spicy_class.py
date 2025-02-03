# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 14:23:42 2024

@author: mendez, ratz
"""

import numpy as np # used in all computations

# these functions are used for the clutering and collocation
from sklearn.neighbors import NearestNeighbors
# Function for the k means clusering
from sklearn.cluster import MiniBatchKMeans

# Note: there is a warning from kmeans when running on windows.
# This should fix it
import warnings
warnings.filterwarnings('ignore')

# Matplotlib for the plotting functions:
import matplotlib.pyplot as plt 

# function useful for computing smallsest and largest eig:
from scipy.sparse.linalg import eigsh
# we use scipy linalg for cholesky decomposition, solving linear systems etc
from scipy import linalg

# We use this function to handle polygonal refinements areas 
from shapely import geometry

    
class spicy:
    """
    SPICY (Super-resolution and Pressure from Image veloCimetrY) is a software 
    developed at the von Karman Institute to perform data assimilation by means 
    of Radial Basis Functions (RBF). The framework works both for structured and 
    unstructered data. Currently, the main application is to perform a regression
    of image velocimetry data and then solve the pressure equation. However, the framework 
    can be readily extended to regression of other fields (e.g. temperature fields).
    
    The original article by Sperotto et al. (2022) can be found at:
    https://arxiv.org/abs/2112.12752
    
    The software paper by Sperotto et al. (2024) can be found at:
    https://joss.theoj.org/papers/10.21105/joss.05749#
    
    YouTube channel with hands-on tutorials can be found at:
    https://www.youtube.com/@spicyVKI    
        
    """    
    # 1. Initialize the class with the data
    def __init__(self, data, grid_point):
        """
        Initialization of an instance of the spicy class. 
                                          
        :type data: list of 1D numpy.ndarray
        :param data:
            This list contains the target data; two arrays [u, v] for a 2D vector field.
                    
        :type grid point: list of 1D numpy.ndarray
        :param grid_point: 
            Is a list of arrays containing the grid points: [X_G , Y_G] in 2D                        
                                   
        General attributes:
            X_G, Y_G: coordinates of the point in which the data is available
            u, v: velocity field to learn, given in the data points
        
        If constraints are assigned:
            X_D, Y_D: coordinates of the points with Dirichlet (D) conditions
            c_D: values of the D conditions
            
            X_N, Y_N: coordinates of the points with Neumann (N) conditions
            n_x, n_y: normal versors where N conditions are introduced
            c_N_X, c_N_Y: values of the N conditions
            
            X_Div, Y_Div: coordinates of the points with Div conditions
                
        If clustering is done:
            r_mM: vector collecting minimum (m) and maximum (M) radious of the RBFs 
            eps_l: scalar controlling the value of an RBF at the closest RBF neighbor               
            X_C, Y_C : coordinates of the cluster centers/collocations 
            c_k: shape parameters of the RBFs 
            d_k: diameters of the rbfs  
       
        If problem is assembled:
            A: matrix A in the linear system
            B: matrix B in the linear system
            b_1: vector b_1 in the linear systems
            b_2: vector b_2 in the linear system       
        
        If computation is done: 
            weights: weights of the RBF regression 
            lambda: Lagrange multipliers of the RBF regression 
            
        """

        # Assign the values to the class
        self.type = '2D'
        self.model = 'laminar'
        # Define the internal grid and velocity data
        self.X_G = grid_point[0]
        self.Y_G = grid_point[1]
        self.u = data[0]
        self.v = data[1]
        self.n_p = len(self.X_G)
        
        return
    
    
    # 2. Clustering (this does not depend on the model, but only on the dimension).
    def clustering(self, n_K, Areas, r_mM=[0.01,0.3], eps_l=0.7):
        """
        This function defines the collocation of a set of RBFs using the multi-
        level clustering first introduced in the article. Note that we modified the slightly original formulation
        to ease the programming; see video tutorials for more.
        The function must be run before the constraint definition.
                 
        :type n_K: list 
        :param n_K:
            This contains the n_k vector in eq (33) in the paper; this is the 
            list of expected particles per RBF at each level. For example, if n_K=[4,10], 
            it means that the clustering will try to have a first level with RBFs whose size
            seeks to embrace 4 points, while the second level seeks to embrace
            10 points, etc. The length of this vector automatically defines the
            number of levels.
            
        :type Areas: list
        :param Areas:
            List of the refinement regions for each clustering level. If no 
            refinement is needed, then this should be a list of empty
            lists (default option).
        
        :type r_mM: list of two float values
        :param r_mM: default=[0.01, 0.3].
            This contains the minimum and the maximum RBF's radiuses. This is
            defined as the distance from the collocation point at which the RBF
            value is 0.5.
                
        :type float: float
        :param eps_l: default=0.7.
            This is the value that a RBF will have at its closest neighbour. It 
            is used to define the shape factor from the clustering results.
                   
        """
        
        # we assign the clustering parameters to self
        # they are needed in the constraints to set the shape parameters for the
        # RBFs which are located at constraint points
        
        self.r_mM = r_mM
        self.eps_l = eps_l
 
        # Number of levels
        n_l = len(n_K)  
        
        # Loop over the number of levels
        for l in range(n_l):
            # We look for the points that belongs to the given area:
            if Areas[l]:
                # This means a polygon object is given, so take only points
                # inside this:
                poly = Areas[l]    
                List = []  # prepare empty list
                for j in range(len(self.X_G)):  # fill list of points in poly
                    List.append(poly.contains(geometry.Point(self.X_G[j], self.Y_G[j])))
                # Take only these points as data matrix 
                X_G_c = self.X_G[List]
                Y_G_c = self.Y_G[List]
                Data_matrix = np.column_stack((X_G_c, Y_G_c))
                List = [] # delete the list for safety
            else: # if Areas is empty then all points should be included
                Data_matrix = np.column_stack((self.X_G, self.Y_G))
            
            # Define number of clusters
            Clust = int(np.ceil(np.shape(Data_matrix)[0]/ n_K[l])) 
            
            # Initialize the cluster function
            model = MiniBatchKMeans(n_clusters=Clust, random_state=0)  
            # Run the clustering and return the indices (optional)
            y_P = model.fit_predict(Data_matrix)
            # Obtaining the centers of the points
            Centers = model.cluster_centers_
            
            # Get the nearest neighbour of each center
            nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(Centers)
            distances, indices = nbrs.kneighbors(Centers)
            sigma1 = distances[:,1]
            
            # Remove all of the clusters which either have a distance of 
            # zero to the nearest neighbor (that would be the same RBF)
            # and the clusters with only one point in them
            count = np.bincount(y_P, minlength = Clust) 
            sigma1[sigma1 == 0] = np.amax(sigma1[sigma1 != 0]) 
            sigma1[count == 1] = np.amax(sigma1) 
            
            # Pre-assign the collocation points
            X_C1 = Centers[:,0]
            Y_C1 = Centers[:,1]
            list_Index=np.array([l]*len(X_C1)) # to use also hstack
            
            # Assign the results to a vector of collocation points
            if l == 0: # If this is the first layer, just assign:
                X_C = X_C1 
                Y_C = Y_C1 
                sigma = sigma1 
                l_list = list_Index
            else: # Stack onto the existing ones
                X_C = np.hstack((X_C, X_C1))
                Y_C = np.hstack((Y_C, Y_C1))
                sigma = np.hstack((sigma, sigma1))
                l_list=np.hstack((l_list, list_Index))
            print('Clustering level ' + str(l) + ' completed')
        
        # Assign to the class
        self.X_C = X_C
        self.Y_C = Y_C
        # For plotting purposes, we keep track of the scale at which
        # the RBF have been places
        self.Clust_list = l_list
          
        # We conclude with the computation of the shape factors
        # Set the max and min values of c_k
        c_min = 1/(r_mM[1])*np.sqrt(np.log(2))
        c_max = 1/(r_mM[0])*np.sqrt(np.log(2))
        # compute the c_k 
        c_k = np.sqrt(-np.log(eps_l))/sigma
        # crop to the minimum and maximum value
        c_k[c_k < c_min] = c_min
        c_k[c_k > c_max] = c_max
        # for plotting purposes, we store also the diameters
        d_k = 2/c_k*np.sqrt(np.log(2))
        
        self.c_k = c_k
        self.d_k = d_k
        
        print(str(len(X_C)) + ' RBFs placed through clustering')
        
        return
        
    # 3. Constraints.
    def vector_constraints(self, DIR=[], NEU=[], DIV=[], CURL=[], extra_RBF = True):
        """        
        # This functions sets the boundary conditions for a laminar problem. The
        function must be run after the clustering was carried out.
        
        :type DIR: list of 1D numpy.ndarray
        :param DIR:
            This contains the info for the Dirichlet conditions.
            In 2D, this has [X_D, Y_D, c_D_X, c_D_Y].
              
            Here X_D, Y_D are the coordinates of the poins where the value c_D_X,
            c_D_Y is set.
        
        :type NEU: list of 1D numpy.ndarray
        :param NEU:
            This contains the info for the Neuman conditions.
            If the model is 2D, then this has [X_N, Y_N, n_x, n_y, c_N_X, c_N_Y].
            
            Here X_N, Y_N are the coordinates of the poins where the value c_N_X,
            c_N_Y is set for the directional derivative along the
            normal direction n_x, n_y
            
        :type DIV: list of 1D numpy.ndarray
        :param DIV:
            This contains the info for the divergence-free conditions.
            In 2D, this has [X_Div, Y_Div].
            
            Here X_Div, Y_Div are the coordinates of the poins where the
            divergence-free condition is imposed.
            
        :type CURL: list of 1D numpy.ndarray
        :param CURL:
            This contains the info for the curl-free conditions.
            In 2D, this has [X_Curl, Y_Curl].
            
            Here X_Curl, Y_Curl are the coordinates of the poins where the
            curl-free condition is imposed.
        
        :type extra_RBF: bool
        :param extra_RBF: default=True
            This is a flag to put extra collocation points where a constraint is
            set. It can improve the solution of the linear system as constraints
            remove degrees of freedom

        """
        
        # Check for Dirichlet conditions
        if not DIR: 
            # We still assign empty arrays so that the assembly of the system is easier
            self.n_D = 0
            self.X_D = np.array([])
            self.Y_D = np.array([])
            self.c_D_X = np.array([])
            self.c_D_Y = np.array([])
        else:
            self.n_D = []
            self.X_D = []
            self.Y_D = []
            self.n_x_D = []
            self.n_y_D = []
            self.c_D = []
            for i in range(len(DIR)):
                self.n_D.append(len(DIR[0][0]))
                self.X_D.append(DIR[i][0])
                self.Y_D.append(DIR[i][1])
                self.n_x_D.append(DIR[i][2])
                self.n_y_D.append(DIR[i][3])
                self.c_D.append(DIR[i][4])
           
        # Check for Neumann conditions
        
        if not NEU:
            # We still assign empty arrays so that the assembly of the system is easier
            self.n_N = 0
            self.X_N = np.array([])
            self.Y_N = np.array([])
            self.c_N_X = np.array([])
            self.c_N_Y = np.array([])
            self.n_x_N = np.array([])
            self.n_y_N = np.array([])
        else:
            self.n_N = len(NEU[0])
            self.X_N = NEU[0]
            self.Y_N = NEU[1]
            self.n_x_N = NEU[2]
            self.n_y_N = NEU[3]
            self.c_N_X = NEU[4]
            self.c_N_Y = NEU[5]
                
        # Check for Divergence conditions
        if not DIV:
            # We still assign empty arrays so that the assembly of the system is easier
            self.n_Div = 0
            self.X_Div = []
            self.Y_Div = []
        else:
            self.n_Div = len(DIV[0])
            self.X_Div = DIV[0]
            self.Y_Div = DIV[1]
                
        # Check for Curl conditions
        if not CURL:
            # We still assign empty arrays so that the assembly of the system is easier
            self.n_Curl = 0
            self.X_Curl = []
            self.Y_Curl = []
        else:
            self.n_Curl = len(CURL[0])
            self.X_Curl = CURL[0]
            self.Y_Curl = CURL[1]
        
        # Finally, we add the extra RBFs in the constraint points if desired
        if extra_RBF == True:
            # Assemble all the points where we have constraints
            X_constr = np.concatenate((self.X_Div, self.X_Curl))
            Y_constr = np.concatenate((self.Y_Div, self.Y_Curl))
            for i in range(len(self.X_D)):
                X_constr = np.concatenate((X_constr, self.X_D[i]))
                Y_constr = np.concatenate((Y_constr, self.Y_D[i]))
            for i in range(len(self.X_N)):
                X_constr = np.concatenate((X_constr, self.X_N[i]))
                Y_constr = np.concatenate((Y_constr, self.Y_N[i]))
            # Get the unique values 
            unique_values = np.unique(np.column_stack((X_constr, Y_constr)), axis = 0)
            X_unique = unique_values[:,0]
            Y_unique = unique_values[:,1]
            # Get the additional RBF shape parameters
            c_k, d_k = add_constraint_collocations(X_unique, Y_unique,
                self.X_C, self.Y_C, self.r_mM, self.eps_l)
            # The new collocation points are stored on the l-th + 1 level
            self.Clust_list = np.concatenate((self.Clust_list, np.ones(c_k.shape)*self.Clust_list[-1]+1))
            # Concatenate them with the existing collocation points
            self.c_k = np.concatenate((self.c_k, c_k))
            self.d_k = np.concatenate((self.d_k, d_k))
            self.X_C = np.concatenate((self.X_C, X_unique))
            self.Y_C = np.concatenate((self.Y_C, Y_unique)) 
        
        # Summary output for the user
        print(str(np.array(self.n_D).sum()) + ' D conditions assigned') 
        print(str(np.array(self.n_N).sum()) + ' N conditions assigned')
        print(str(self.n_Div) + ' Div conditions assigned')
        print(str(self.n_Curl) + ' Curl conditions assigned')
        
        self.n_b = self.c_k.shape[0]
        
        return
        
    # 3.3 Plot the RBFs, this is just a visualization tool
    def plot_RBFs(self,l=0):
        """
        Utility function to check the spreading of the RBFs after the clustering.
        This function generates several plots. It produces no new variable in SPICY. 
        
        :type l: int
        :param l: 
            This defines the cluster level of RBF that will be visualized.

        """

        try:
            # We define the data that will be included
            X_Plot = self.X_C[np.argwhere(self.Clust_list==l)]
            Y_Plot = self.Y_C[np.argwhere(self.Clust_list==l)]
            d_K_Plot = self.d_k[np.argwhere(self.Clust_list==l)]

            fig, axs = plt.subplots(1, 2, figsize = (10, 5), dpi = 100)
            # First plot is the RBF distribution
            axs[0].set_title("RBF Collocation for l="+str(l))

            # Also show the data points
            axs[0].scatter(self.X_G, self.Y_G, c=np.sqrt(self.u**2 + self.v**2), s=10)

            for i in range(0,len(X_Plot),1):
                circle1 = plt.Circle((X_Plot[i], Y_Plot[i]), d_K_Plot[i]/2,
                                      fill=True, facecolor='g', edgecolor='k', alpha=0.1)
                axs[0].add_artist(circle1)

            # Also show the constraints if they are set
            for i in range(len(self.X_D)):
                axs[0].plot(self.X_D[i], self.Y_D[i],'ro')

            for i in range(len(self.X_N)):
                axs[0].plot(self.X_N[i], self.Y_N[i],'bs')

            axs[0].plot(self.X_Div, self.Y_Div, 'bd')

            axs[0].plot(self.X_Curl, self.Y_Curl, 'bv')

            # Second plot is the distribution of diameters:
            axs[1].stem(d_K_Plot)
            axs[1].set_xlabel('Basis index')
            axs[1].set_ylabel('Diameter')
            axs[1].set_title("Distribution of diameters for L="+str(l))
            fig.tight_layout()
            plt.show()

        except:
            raise ValueError('Problems in plotting. Set constraints and cluster!')

        return


    # 4. Assembly A, B, b_1, b_2
    def Assembly_Regression(self, alpha_div=None):
        """
        This function assembly the matrices A, B, C, D from the paper (see video tutorial 1).
           
        :type alpha_div: float
        :param alpha_div:
            This enables a divergence free penalty in the entire flow field.
            The higher this parameter, the more SPICY penalizes errors in the divergence-free 
            condition. This is particularly important to obtain good derivatives 
            for the pressure computation.
            
         """   
                  
        # Define the rescaling factor which is done based on the maximum
        # absolute velocity that is available in u and v
        data = np.concatenate((self.u, self.v))
        self.rescale = data[np.argmax(np.abs(data))]
        
        
        ### Curl-free constraints ###
        # Compute Phi_x on X_Curl
        Matrix_Phi_X_Curl_der_x = Phi_RBF_x(self.X_Curl, self.Y_Curl, self.X_C, self.Y_C, self.c_k)
        # Compute Phi_x on X_Curl
        Matrix_Phi_X_Curl_der_y = Phi_RBF_y(self.X_Curl, self.Y_Curl, self.X_C, self.Y_C, self.c_k)
        # Stack into block structure
        Matrix_D_Curl = np.hstack((-Matrix_Phi_X_Curl_der_y, Matrix_Phi_X_Curl_der_x))
        
        
        ### Divergence-free constraints ###
        # Compute Phi_x on X_Div
        Matrix_Phi_X_Div_der_x = Phi_RBF_x(self.X_Div, self.Y_Div, self.X_C, self.Y_C, self.c_k)
        # compute the derivatives in y
        Matrix_Phi_X_Div_der_y = Phi_RBF_y(self.X_Div, self.Y_Div, self.X_C, self.Y_C, self.c_k)
        # Stack into the block structure of equation (15)
        Matrix_D_Div = np.hstack((Matrix_Phi_X_Div_der_x, Matrix_Phi_X_Div_der_y)) 
        
        
        ### Dirichlet constraints ###
        Matrix_D = np.empty((0, 2*self.n_b))
        b_2_D = np.array([])
        # Loop over all the boundaries
        for i in range(len(self.X_D)):
            # Get Phi on X_D[i]
            Matrix_Phi = Phi_RBF(self.X_D[i], self.Y_D[i], self.X_C, self.Y_C, self.c_k)
            # Project Phi onto the x and y direction
            Matrix_Phi_Proj = np.hstack((
                Matrix_Phi*self.n_x_D[i][:, np.newaxis],
                Matrix_Phi*self.n_y_D[i][:, np.newaxis]
                ))
            # Concatenate the matrix block and the constraints
            Matrix_D = np.concatenate((Matrix_D, Matrix_Phi_Proj))
            b_2_D = np.concatenate((b_2_D, self.c_D[i]))
        # # Concatenate the constraints
        # b_2_D = np.concatenate(self.c_D)
        
        
        ### Neumann constraints ###
        # Compute Phi_x on X_N[i]
        Matrix_Phi_X_N_der_x = Phi_RBF_x(self.X_N, self.Y_N, self.X_C, self.Y_C, self.c_k)
        # Compute Phi_y on X_N[i]
        Matrix_Phi_X_N_der_y = Phi_RBF_y(self.X_N, self.Y_N, self.X_C, self.Y_C, self.c_k)
        # Project onto the x and y direction
        Matrix_Phi_N = Matrix_Phi_X_N_der_x*self.n_x_N[:, np.newaxis] +\
                       Matrix_Phi_X_N_der_y*self.n_y_N[:, np.newaxis]
        Matrix_D_N = np.block([
            [Matrix_Phi_N,np.zeros((self.n_N, self.n_b))],
            [np.zeros((self.n_N, self.n_b)), Matrix_Phi_N]
            ])
                
        # Assemble B and b_2, we also rescale b_2
        self.B = np.vstack((Matrix_D_Div, Matrix_D_Curl, Matrix_D, Matrix_D_N)).T
        self.b_2 = np.concatenate((np.zeros(self.n_Div),
                                   np.zeros(self.n_Curl),
                                   b_2_D,
                                   self.c_N_X, self.c_N_Y)) / self.rescale
        
        # Compute Phi on X_G
        Matrix_Phi_X = Phi_RBF(self.X_G, self.Y_G, self.X_C, self.Y_C, self.c_k)
        # Stack Phi.T@Phi into the block structure of equation (10)
        PhiT_dot_Phi = Matrix_Phi_X.T.dot(Matrix_Phi_X)
        self.A = 2*np.block([
            [PhiT_dot_Phi, np.zeros((self.n_b, self.n_b))],
            [np.zeros((self.n_b, self.n_b)), PhiT_dot_Phi]
            ])
        # Compute and rescale b_1
        self.b_1 = 2*np.concatenate((Matrix_Phi_X.T.dot(self.u),
                                     Matrix_Phi_X.T.dot(self.v))) / self.rescale
        
        # We check if alpha_div is None or 0
        if alpha_div is not None and alpha_div != 0: 
            # Compute Phi_x on X_G
            Matrix_Phi_X_der_x = Phi_RBF_x(self.X_G, self.Y_G, self.X_C, self.Y_C, self.c_k)
                
            # Compute Phi_y on X_G
            Matrix_Phi_X_der_y = Phi_RBF_y(self.X_G, self.Y_G, self.X_C, self.Y_C, self.c_k)
            
            # Compute the individual matrix products between x, y and z
            # For the diagonal
            PhiXT_dot_PhiX = Matrix_Phi_X_der_x.T.dot(Matrix_Phi_X_der_x)
            PhiYT_dot_PhiY = Matrix_Phi_X_der_y.T.dot(Matrix_Phi_X_der_y) 
            # For the off-diagonal elements
            PhiXT_dot_PhiY = Matrix_Phi_X_der_x.T.dot(Matrix_Phi_X_der_y)
            
            # And we add them into the A-matrix
            # Diagonal
            self.A[self.n_b*0:self.n_b*1,self.n_b*0:self.n_b*1] += 2*alpha_div*PhiXT_dot_PhiX
            self.A[self.n_b*1:self.n_b*2,self.n_b*1:self.n_b*2] += 2*alpha_div*PhiYT_dot_PhiY
            # Upper off-diagonal elements
            self.A[self.n_b*0:self.n_b*1,self.n_b*1:self.n_b*2] += 2*alpha_div*PhiXT_dot_PhiY
            # Lower off-diagonal elements
            self.A[self.n_b*1:self.n_b*2,self.n_b*0:self.n_b*1] += 2*alpha_div*PhiXT_dot_PhiY.T 
            
        return
        
    
    # 5 Solver using the Shur complement
    def Solve(self, K_cond=1e12):
        """
        This function solves the constrained quadratic problem A, B, b_1, b_2.
    
        The input parameters are the class itself and the desired condition 
        number of A which is fixed based on its largest and smallest eigenvalue
        
        The function assigns the weights 'w' and the Lagrange multipliers
        Lambda to the class. The weights are computed for the min/max scaled problem,
        i.e. the right hand-side of the linear system is normalized. The assigned
        weights are rescaled by self.rescale to get the real, physical quantities
        
        :type K_cond: float
        :param K_cond: Default 1e12.
          This is the regularization parameter. It fixes the condition number (see Video 1)
          The estimation is based such that the regularize matrix has the condition
          number k_cond. For this, we compute the max and the min eigenvalue.
        """
            
        # Step 1: Regularize the matrix A
        lambda_A = eigsh(self.A, 1, return_eigenvectors=False) # Largest eigenvalue
        alpha = (lambda_A-K_cond*2.2e-16) / K_cond
        self.A = self.A + alpha*np.eye(np.shape(self.A)[0])
        print('Matrix A regularized')

        # Step 2: Cholesky Decomposition of A
        L_A, low = linalg.cho_factor(self.A, overwrite_a=True, check_finite=False, lower=True)

        # Step 3: Solve for N
        N = linalg.cho_solve((L_A,low), self.B, check_finite=False)

        # Step 4: prepare M
        M = N.T@self.B

        # Step 5 + 6: Regularize M if needed, then compute chol factor
        try:
            # try without regularization
            L_M, low = linalg.cho_factor(M, overwrite_a=True, check_finite=False, lower=True)
            print('Chol factor of M WITHOUT regularization')
        except:
            # if it does not work, regularize M the same way as for A
            lambda_M = eigsh(M, 1, return_eigenvectors=False) # Largest eigenvalue
            alpha = (lambda_M - K_cond*2.2e-16) / K_cond
            M = M + alpha*np.eye(np.shape(M)[0])
            L_M, low = linalg.cho_factor(M, overwrite_a = True, check_finite = False, lower = True)
            print('Chol factor of M WITH regularization')

        # Step 7: Solve the system for lambda
        b2star = N.T.dot(self.b_1) - self.b_2
        self.lambdas = linalg.cho_solve((L_M, low), b2star, check_finite = False)
        print('Lambdas computed')

        # Step 8: Solve for w.
        b1_star = self.b_1 - self.B.dot(self.lambdas)
        self.weights = linalg.cho_solve((L_A, low), b1_star, check_finite=False) * self.rescale
        print('w computed')

        return 


    # 6. Evaluate solution on arbitrary grid
    
    # Here is a function to compute the solution on an arbitrary set of points
    # X_G, Y_G. We take w, lam from the solution, X_C, Y_C, c_k from the clustering.
    def Get_Sol(self, grid):
        """
        This function evaluates the solution of the linear system on an arbitrary
        set of points on the grid.
        
        :type K_cond: list
        :param grid:
            Contains the points at which the source term is evaluated
            If the model is 2D, then this has [X_P, Y_P].
            
        :return: U_sol, V_sol.  
            If laminar and 2D, the solution is U_sol, V_sol
        
        """   

        # Assign the grid
        X_P = grid[0]
        Y_P = grid[1]

        # Evaluate Phi on the grid X_P
        Phi = Phi_RBF(X_P, Y_P, self.X_C, self.Y_C, self.c_k)

        # Compute the individual components
        U_sol = Phi.dot(self.weights[0*self.n_b:1*self.n_b])
        V_sol = Phi.dot(self.weights[1*self.n_b:2*self.n_b])

        return U_sol, V_sol
         
    def Get_first_Derivatives(self, grid):
        """ 
        This function evaluates the first derivative of the solution of the
        linear system on an arbitrary set of points on the grid.
        
        :type grid: list
        :param grid:
            Contains the points at which the source term is evaluated
            If the model is 2D, then this has [X_P, Y_P].

        :return: dUdX, dUdY, dVdX, dVdY.
            If laminar and 2D, the solution is dUdx, dUdY, dVdX, dVdY
            
        """ 

        # Assign the grid
        X_P = grid[0]
        Y_P = grid[1]

        # We do it in 2 blocks: first all derivatives in x
        # Evaluate Phi on the grid X_P, Y_P
        Phi_deriv = Phi_RBF_x(X_P, Y_P, self.X_C, self.Y_C, self.c_k)

        # Compute dudx and dvdx on the new grid
        dUdX = Phi_deriv.dot(self.weights[0*self.n_b:1*self.n_b])
        dVdX = Phi_deriv.dot(self.weights[1*self.n_b:2*self.n_b])

        # Then we do it again for the derivatives of y.
        # Note however, that we re-use the same variables Phi_deriv
        # to limit the memory usage. This is pretty much copy-paste.
        Phi_deriv = Phi_RBF_y(X_P, Y_P, self.X_C, self.Y_C, self.c_k)

        # Compute dudy and dvdy on the new grid
        dUdY = Phi_deriv.dot(self.weights[0*self.n_b:1*self.n_b])
        dVdY = Phi_deriv.dot(self.weights[1*self.n_b:2*self.n_b])

        return dUdX, dUdY, dVdX, dVdY
            
# =============================================================================
#  Utilities functions
#  These functions are not needed/called by the user. They are simply helper 
#  functions required to assemble and solve the linear systems. In the current
#  release of SPICY, these are:
#  - RBF functions and their derivatives in 2D   
#  - Adding collocation points in the constraints 
# =============================================================================


def Phi_RBF(X_G, Y_G, X_C, Y_C, c_k):
    """
    Get the basis matrix at the points (X_G,Y_G) from RBFs at the collocation points
    at (X_C,Y_C), having shape factors c_k. The output is a matrix of side (n_p) x (n_c).
    The basis can be 'c4' or 'gauss'.
    """
    # This is the contribution of the RBF part
    n_b = len(X_C); n_p = len(X_G)
    Phi_RBF = np.zeros((n_p, n_b))

    # Iterate over all basis elements
    for r in range(n_b):
        # Compute the Gaussian
        gaussian = np.exp(-c_k[r]**2 * ((X_G-X_C[r])**2 + (Y_G-Y_C[r])**2))
        # Assemble into matrix
        Phi_RBF[:,r] = gaussian
    
    # Return the matrix
    return Phi_RBF


def Phi_RBF_x(X_G, Y_G, X_C, Y_C, c_k):
    """
    Get the derivative along x of the basis matrix at the points (X_G,Y_G) from 
    RBFs at the collocation points at (X_C,Y_C), having shape factors c_k. The
    output is a matrix of side (n_p) x (n_c). The basis can be 'c4' or 'gauss'.
    """
    # number of bases (n_b) and points (n_p)
    n_b = len(X_C); n_p = len(X_G)
    # Initialize the matrix
    Phi_RBF_x = np.zeros((n_p, n_b))

    # Iterate over all basis elements
    for r in range(n_b):
        # Compute the Gaussian
        gaussian = np.exp(-c_k[r]**2 * ((X_G-X_C[r])**2 + (Y_G-Y_C[r])**2))
        # Multiply with inner term and assemble into matrix
        Phi_RBF_x[:,r] = - 2*c_k[r]**2 * (X_G-X_C[r])*gaussian

    # Return the matrix
    return Phi_RBF_x


def Phi_RBF_y(X_G, Y_G, X_C, Y_C, c_k):
    """
    Get the derivative along y of the basis matrix at the points (X_G,Y_G) from 
    RBFs at the collocation points at (X_C,Y_C), having shape factors c_k. The
    output is a matrix of side (n_p) x (n_c). The basis can be 'c4' or 'gauss'.
    """
    # number of bases (n_b) and points (n_p)
    n_b = len(X_C); n_p = len(X_G)
    # Initialize the matrix
    Phi_RBF_y = np.zeros((n_p, n_b))

    # Iterate over all basis elements
    for r in range(n_b):
        # Compute the Gaussian
        gaussian = np.exp(-c_k[r]**2 * ((X_G-X_C[r])**2 + (Y_G-Y_C[r])**2))
        # Multiply with inner term and assemble into matrix
        Phi_RBF_y[:,r] = - 2*c_k[r]**2 * (Y_G-Y_C[r])*gaussian
            
    # Return the matrix
    return Phi_RBF_y

def add_constraint_collocations(X_constr, Y_constr, X_C, Y_C, r_mM, eps_l):
    """
    This function adds collocation points where constraints are set in 2D.
    
    ----------------------------------------------------------------------------------------------------------------
    Parameters
    ----------
    :param X_constr: np.ndarray
        X coordinates of the constraints
    :param Y_constr: np.ndarray
        Y coordinates of the constraints
    :param X_C: np.ndarray
        X coordinates of the collocation points
    :param Y_C: np.ndarray
        Y coordinates of the collocation points
    :param r_mM: list
        Minimum and maximum radius of the RBFs
    :param eps_l: float
        Value of the RBF at its closest neighbor
    """   
    # Get the number of constraints
    n_constr = X_constr.shape[0]
    # Initialize an empty array for the shape parameters
    c_ks = np.zeros(n_constr)
    
    # Set the max and min values of c_k  
    c_min = 1 / (r_mM[1]) * np.sqrt(np.log(2))
    c_max = 1 / (r_mM[0]) * np.sqrt(np.log(2))
    # Loop over all constraints
    for k in range(n_constr):
        # Get the distance to all collocation points
        dist_to_colloc = np.sqrt((X_C - X_constr[k])**2 + (Y_C - Y_constr[k])**2)
        # Get the distance to all constraints, except for itself
        dist_to_constr = np.sqrt((np.delete(X_constr, k) - X_constr[k])**2+\
                                 (np.delete(Y_constr, k) - Y_constr[k])**2)
        # Set the max and min values of c_k 
        c_k = np.sqrt(-np.log(eps_l)) / np.concatenate((dist_to_colloc, dist_to_constr))
        # crop to the minimum and maximum value
        c_k[c_k < c_min] = c_min
        c_k[c_k > c_max] = c_max
        # get the maximum value in the case of the Gaussian
        c_ks[k] = np.max(c_k)
    # for plotting purposes, we store also the diameters             
    d_k = 2/c_ks*np.sqrt(np.log(2))      
    
    return c_ks, d_k