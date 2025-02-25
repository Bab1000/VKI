import numpy as np
from scipy.signal import stft
from scipy import signal
from scipy.signal import filtfilt
import matplotlib.pyplot as plt
import glob
from colorama import Fore, Style, init
import os
from modulo_vki.utils.others import Plot_Field_TEXT_Cylinder
from tqdm import tqdm
import imageio as imageio
 
def LibraryManager():
   # List of required libraries
    libraries = [
        ("numpy", "np"),
        ("scipy.signal", "signal"),
        ("matplotlib.pyplot", "plt"),
        ("glob", None),
        ("colorama", None),
        ("os", None),
        ("modulo_vki.utils.others", "Plot_Field_TEXT_Cylinder"),
        ("tqdm", "tqdm"),
        ("imageio", "imageio")
    ]
 
    # Try to import each library and handle missing ones
    missing_libs = []
    for lib, alias in libraries:
        try:
            if alias:
                exec(f"import {lib} as {alias}")
            else:
                exec(f"import {lib}")
        except ImportError:
            missing_libs.append(lib)
            print(Fore.RED + f"ERROR: The library '{lib}' is missing! Please install it.")
            exit(1)  # Exit the script with an error
 
    # If all libraries are available, continue execution
    print("")
    print(Fore.GREEN + "All required libraries are installed. Proceeding...")
    print("")
 
    # Initialize colorama
    init(autoreset=True)
 
def ResFileProcessing(path,coord):
 
    # Processing the Res files to get the velocities at the required location
    print(Fore.WHITE + "--> Processing the 2D cylinder results files...")
 
    # Data path
    data_path = path + "/Tutorial_5_2D_Cylinder_Memory_Saving/data"
 
 
    # Pattern of the Res file name
    name_pattern = "/Res*.dat"
 
    # Creating a list with all the Res files
    Res_file_list = sorted(glob.glob(data_path+name_pattern))
 
    # Name of the mesh file
    mesh_file = data_path + "/MESH.dat"
 
    # Looking for the index of the coordinates to analyse
    with open(mesh_file,'r') as file:
        for index, line in enumerate(file):
           
            if index != 0:
           
                # Get the coordinates of each lines in the file
                mesh_coord = line.strip().split()
 
                mesh_coord_0 = float(mesh_coord[0])
                mesh_coord_1 = float(mesh_coord[1])
 
                if mesh_coord_0 == coord[0] and mesh_coord_1 == coord[1]:
 
                    coord_index = index
 
    # Gathering all the velocities
    U = []
    V = []
    time_list = []
    time = 0
 
    for file in Res_file_list:
 
        with open(file,'r') as file:
 
            # Adding the velocity in the list
            velocities = file.readlines()
            coord_vel = velocities[coord_index].strip().split()
            U.append(float(coord_vel[0]))
            V.append(float(coord_vel[1]))
 
    return U,V
 
def ButterworthHighpass(cutoff_freq, fs, order):
    b, a = signal.butter(order, cutoff_freq, 'hp', analog = False, fs=fs)
    return b, a
 
def Exercise1(u,v,fs,n_samples):
 
    # Computing the magnitude of the velocities
    magnitude = []
    for i in range(len(u)):
        magnitude.append(np.sqrt(u[i]**2 + v[i]**2))
   
 
    # Apply the high-pass filter to remove low frequencies
    cutoff_freq = 100  # Cutoff frequency
    order = 1         # order of the butterworth
    b, a = ButterworthHighpass(cutoff_freq, fs, order)
    magnitude_filtered = filtfilt(b, a, magnitude)  # Zero-phase filtering
 
    print(Fore.WHITE + "--> Computing the fft of the signal ...")
    # Computing the fft of the magnitude
    delta_t = 1/fs
 
    fft=np.fft.fft(magnitude_filtered)/np.sqrt(n_samples)
    freq=np.fft.fftfreq(n_samples,delta_t)
    freq_shift = np.fft.fftshift(freq)
    fft_shift = np.fft.fftshift(fft)
 
    return freq_shift,fft_shift,magnitude,magnitude_filtered
 
def PlotExercise1(freq,fft,path,magnitude,n_samples,fs):
 
    print(Fore.WHITE + "--> Plotting the results of Exercise 1 ...")
   
    # Creating the results path for Exercise 1
    os.makedirs(path + "/Results_exercise_1", exist_ok=True)
 
    # Creating the name of the pdf file to save the figure
    save_path_fft = os.path.join(path, "Results_exercise_1", "fft.png")
   
    # Ploting the results
    plt.figure()
    plt.plot(freq, np.abs(fft))
    plt.ylim([0, 8])
    plt.xlim([0, 600])
    plt.xlabel("f [Hz]")  # Label de l'axe x
    plt.ylabel("Normalized Amplitude")  # Label de l'axe y
    plt.grid(True)  # Ajout d'une grille pour une meilleure lisibilité
    # Saving the results in the pdf file
    plt.savefig(save_path_fft, format='png', dpi=300, bbox_inches='tight')
    plt.close()

    # Total time of the signal
    t_tot = n_samples * 1 / fs
    t_plot = np.linspace(0, t_tot, n_samples, endpoint=False)

    save_path_mag = os.path.join(path, "Results_exercise_1", "mag.png")
    # Ploting the results
    plt.figure(figsize=(5,3))
    plt.plot(t_plot, magnitude)
    plt.xlabel("time [s]")  # Label de l'axe x
    plt.ylabel("Velocity magnitude [m/s]")  # Label de l'axe y
    plt.grid(True)  # Ajout d'une grille pour une meilleure lisibilité
    # Saving the results in the pdf file
    plt.savefig(save_path_mag, format='png', dpi=300, bbox_inches='tight')
    plt.close()
   
    print(Fore.GREEN + f"The results have been successfully saved at: {save_path_fft} | {save_path_mag}")
 
def Exercise2(magnitude,fs):
    #Time frequency analysis
    print(Fore.WHITE + "--> Computing the sfft of the signal for time frequency analysis ...")
 
    # STFT to analyse frequency content over time
    freq, t, Amplitude = stft(magnitude, fs, window='hann',nperseg=200)
 
    return freq, t, Amplitude
 
def PlotExercise2(freq,t,Amplitude,path):
 
    print(Fore.WHITE + "--> Plotting the results of Exercise 2 ...")
   
    # Creating the results path for Exercise 2
    os.makedirs(path + "/Results_exercise_2", exist_ok=True)
 
    # Creating the name of the pdf file to save the figure
    save_path_sfft = os.path.join(path, "Results_exercise_2", "sfft.jpeg")
 
    # Plot the magnitude spectrum
    plt.figure(figsize=(5,3))
    plt.pcolormesh(t, freq, np.abs(Amplitude), shading='auto')
    plt.ylabel("Frequency [Hz]")
    plt.xlabel("Time [s]")
    plt.colorbar(label="Magnitude")
    # Saving the results in the pdf file
    plt.savefig(save_path_sfft, format='jpeg', dpi=600, bbox_inches='tight', transparent=True)
    plt.close()
 
    print(Fore.GREEN + f"The results have been successfully saved at: {save_path_sfft}")
 
def Exercise3(magnitude,fs,n_samples):
 
    # Multi-Resolution Analysis of the signal
    print(Fore.WHITE + "--> Computing the Multi-Resolution Analysis of the signal ...")
 
    # Frequency splitting vector for MRA analysis
    freq_splitting_vector = np.array([15, 270, 330, 420, 480])
 
    # Calling MRA function from Signal processing class (Mendez fct)
    U_MRA, H_MRA = MRA_SISO(magnitude, freq_splitting_vector, fs, 200)
 
    print(Fore.WHITE + "--> MRA analysis successfully executed !")
 
    return U_MRA, H_MRA, freq_splitting_vector
 
def MRA_SISO(u,F_V,f_s,N):
  """
  This function computes the MRA of a signal u
  using Hamming Windows  
  :param u: Input Signal
  :param F_V: Frequency Splitting Vectors (see notes)
  :param f_s: Sampling Frequency
  :param N: Order of the Filter  
  :return: U_MRA, n_M scale partitions of the signal
           H_MRA  n_M scale Amplitude Responses Functions
  """
  # Get number of scales and number of points
  n_M=len(F_V)+1; n_t=len(u)
  # Initialize the output
  U_MRA=np.zeros((n_t,n_M))
  # Initialize the Transfer Function Matrix
  H_MRA=np.zeros((512,n_M))  
  # Loop from highest frequency to lowest:
  F_V_o=-np.sort(-F_V)
  # Loop over the scales
  for m in range(0,n_M):
   if m==0:
    #Create first low pass kernel    
    h=signal.firwin(N, F_V_o[m], pass_zero=True, fs=f_s)
    # Get frequency response (for checking)
    w,H_L=signal.freqz(h)
    # This is the first large scale
    u_L=signal.fftconvolve(u,h,'same'); u_H=u-u_L;
    U_MRA[:,m]=u_H
    H_MRA[:,m]=1-abs(H_L)
   elif m>0 and m<n_M-1:
    #Create mth low pass kernel    
    h=signal.firwin(N, F_V_o[m], pass_zero=True, fs=f_s)
    w,H_L_new=signal.freqz(h)
    #Low-pass filter m+1
    u_L_new=signal.fftconvolve(u,h,'same')
    # Band Pass filtered and store
    u_H=u_L-u_L_new; U_MRA[:,m]=u_H
    u_L=u_L_new
    # Get Frequency Transf Function and store
    H_MRA[:,m]=abs(H_L)-abs(H_L_new); H_L=H_L_new
   else:
    U_MRA[:,m]=u_L_new
    # Get Frequency Transf Function and store
    H_MRA[:,m]=abs(H_L);
 
  return U_MRA,H_MRA
 
def PlotExercise3(U_MRA,freq_splitting_vector,n_samples,fs, path):
 
    print(Fore.WHITE + "--> Plotting the results of Exercise 3 ...")
   
    # Creating the results path for Exercise 3
    os.makedirs(path + "/Results_exercise_3", exist_ok=True)
 
    # Total time of the signal
    t_tot = n_samples * 1 / fs
    t_plot = np.linspace(0, t_tot, n_samples, endpoint=False)
 
    for i in range(len(freq_splitting_vector) + 1):
 
        # Creating the name of the pdf file to save the figure
        save_path_MRA = os.path.join(path, "Results_exercise_3", f"MRA_range_{i}.pdf")
       
        plt.figure()  
        plt.plot(t_plot, U_MRA[:, i])
        plt.xlabel("Time (s)")
        plt.ylabel("Amplitude")
        plt.title(f"Multi-Resolution Analysis range {i}")
        plt.grid(True)  
        plt.savefig(save_path_MRA, format='pdf', dpi=600, bbox_inches='tight', transparent=True)
        plt.close()
 
    # Generation of the gifs
    GifGenerator(path,freq_splitting_vector,fs,n_samples)
 
def GifGenerator(path,freq_splitting_vector,fs,n_samples):
 
    print(Fore.WHITE + "--> Generation of the animation of the different modes ...")
 
    # Load the data from Snapshot_Matrices.npz file
    data = np.load('Snapshot_Matrices.npz')
 
    # Access the stored arrays
    D_U = data['D_U']
    D_V = data['D_V']
    Xg = data['Xg']
    Yg = data['Yg']
 
    # Number of velocity points
    n_x,n_y=np.shape(Xg)
 
    # Number of points in one image
    n_p = np.shape(D_U)[0]
 
    # Computing the magnitude Mag=np.sqrt(D_U**2+D_V**2)
    Mag=np.sqrt(D_U**2+D_V**2)
 
    # Prepare the range which needs to be plotted
    scale1=np.zeros((n_samples,n_p))
    scale2=np.zeros((n_samples,n_p))
    scale3=np.zeros((n_samples,n_p))
 
    # MRA analysis of each point of each image
    print(Fore.WHITE + "--> MRA analysis of the full images ...")
    for i in np.arange(0,n_p,1):
        # Gives point by point in all the images into the function to analyse the frequency content of each point in the image
        U_MRA, _= MRA_SISO(Mag[i,:],freq_splitting_vector,fs,200)
        scale1[:,i]=U_MRA[:,-1]
        scale2[:,i]=U_MRA[:,-3]
        scale3[:,i]=U_MRA[:,-5]
 
    print(Fore.WHITE + "--> MRA analysis successfully executed !")
 
    # -----------------------------------------------------------------------------------
 
    # Generation of images scale 1
 
    # Creating images storing directory
    os.makedirs(path + "/Results_exercise_3/Images/Scale1", exist_ok=True)
 
    # Save path of the figure
    save_path_image = os.path.join(path, "Results_exercise_3", "Images","Scale1")
 
    Images_list=[]
    for i in tqdm(np.arange(0, n_samples, 300), desc="--> Generating images for animation of scale 1", unit="img", ncols=120):
 
        plt.figure()  
        # Handle NaN values by replacing them with interpolated values or zero
        mode_data = scale1[i, :].reshape(n_y, n_x).T  
        mode_data = np.nan_to_num(mode_data, nan=np.nanmean(mode_data))  # Replace NaN with the mean value
        vmin, vmax = mode_data.min(), mode_data.max()
        plt.contourf(Xg, Yg, mode_data, levels=np.linspace(0, 14.5, 500), cmap="jet")
        plt.colorbar(label="Partition velocity")
        filename = save_path_image + f"/{i}_1.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight', format='png')
        plt.close()
        Images_list.append(save_path_image + "/" + str(i)+'_1'+'.png')
 
    with imageio.get_writer(path + "/Results_exercise_3/Gif_scale1.gif", mode='I', duration=0.7) as writer:
        for image in Images_list:
            img = imageio.imread(image)  # Leggi ogni immagine
            writer.append_data(img)
 
    print(Fore.GREEN + "--> Gif succesfully generated for scale 1 !")
 
    # -----------------------------------------------------------------------------------
 
    # Generation of images scale 2
 
    # Creating images storing directory
    os.makedirs(path + "/Results_exercise_3/Images/Scale2", exist_ok=True)
 
    # Save path of the figure
    save_path_image = os.path.join(path, "Results_exercise_3", "Images","Scale2")
 
    Images_list=[]
    for i in tqdm(np.arange(0, n_samples, 300), desc="--> Generating images for animation of scale 2", unit="img", ncols=120):
 
        plt.figure()  
        # Handle NaN values by replacing them with interpolated values or zero
        mode_data = scale2[i, :].reshape(n_y, n_x).T  
        mode_data = np.nan_to_num(mode_data, nan=np.nanmean(mode_data))  # Replace NaN with the mean value
        vmin, vmax = mode_data.min(), mode_data.max()
        plt.contourf(Xg, Yg, mode_data, levels=np.linspace(-3, 2.5, 500), cmap="jet")
        plt.colorbar(label="Partition velocity")
        filename = save_path_image + f"/{i}_2.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight', format='png')
        plt.close()
        Images_list.append(save_path_image + "/" + str(i)+'_2'+'.png')
 
    with imageio.get_writer(path + "/Results_exercise_3/Gif_scale2.gif", mode='I', duration=0.7) as writer:
        for image in Images_list:
            img = imageio.imread(image)  # Leggi ogni immagine
            writer.append_data(img)
 
    print(Fore.GREEN + "--> Gif succesfully generated for scale 2 !")
 
    # --------------------------------------------------------------------------------------
 
    # Generation of images scale 3
 
    # Creating images storing directory
    os.makedirs(path + "/Results_exercise_3/Images/Scale3", exist_ok=True)
 
    # Save path of the figure
    save_path_image = os.path.join(path, "Results_exercise_3", "Images","Scale3")
 
    Images_list=[]
    for i in tqdm(np.arange(0, n_samples, 300), desc="--> Generating images for animation of scale 3", unit="img", ncols=120):
 
        plt.figure()  
        # Handle NaN values by replacing them with interpolated values or zero
        mode_data = scale3[i, :].reshape(n_y, n_x).T  
        mode_data = np.nan_to_num(mode_data, nan=np.nanmean(mode_data))  # Replace NaN with the mean value
        vmin, vmax = mode_data.min(), mode_data.max()
        plt.contourf(Xg, Yg, mode_data, levels=np.linspace(-4, 3.7, 500), cmap="jet")
        plt.colorbar(label="Partition velocity")
        filename = save_path_image + f"/{i}_3.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight', format='png')
        plt.close()
        Images_list.append(save_path_image + "/" + str(i)+'_3'+'.png')
 
    with imageio.get_writer(path + "/Results_exercise_3/Gif_scale3.gif", mode='I', duration=0.7) as writer:
        for image in Images_list:
            img = imageio.imread(image)  # Leggi ogni immagine
            writer.append_data(img)
 
    print(Fore.GREEN + "--> Gif succesfully generated for scale 3 !")
 