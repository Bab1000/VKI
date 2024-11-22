import pandas as pd
import numpy as np

# import data
matT = pd.read_csv("/home/jpe/VKI/ELAB/2023-2024/matrix_T.csv").to_numpy()
matH = pd.read_csv("/home/jpe/VKI/ELAB/2023-2024/matrix_h.csv").to_numpy()

def nat_conv_interp(yIndex:int, T:float) -> float:
    """Computes h(y,T) by linear interpolation, from experimental curves
    for natural convection h(y) and T(y)

    Parameters
    ----------
    yIndex : int
        index of the y position
    T : float
        temperature at this position

    Returns
    -------
    float
        natural convection heat transfer coefficient
    """

    # values corresponding to this y index
    T_y = matT[yIndex,:]
    h_y = matH[yIndex,:]

    # find closest T index
    diff = T - T_y
    i = len(diff[diff >= 0]) - 1
    if i < 0:
        return h_y[0]
    
    if i == len(T_y)-1:
        return h_y[i]

    # ratio for interpolation
    ratio = (T - T_y[i])/(T_y[i+1] - T_y[i])

    # interpolated value for h
    h_nat = h_y[i] + ratio*(h_y[i+1] - h_y[i])

    return h_nat