"""
Thermocouple Data Acquisition with NI DAQ and Real-Time Plotting

Description:
This script acquires temperature data from two thermocouples connected to a National Instruments (NI) cDAQ device
with an NI 9211 module, using Python and the NI-DAQmx library.
The data is read continuously from both channels at a defined sampling rate,
saved to a text file with timestamps, and plotted in real-time.

Components:
- **Data Acquisition**: Configured for two channels (ai0 and ai1) with K-type thermocouples.
- **File Saving**: Saves acquired data with timestamps to "thermocouple_data.txt".
- **Real-Time Plotting**: Displays the temperature readings from both channels on a matplotlib plot that updates in real-time.

Setup Requirements:
1. NI cDAQ chassis with an NI 9211 module.
2. Python environment with NI-DAQmx, numpy, matplotlib, and other dependencies.

Acquisition Parameters:
- `Fs`: Sampling rate in Hz (must be compatible with NI 9211's specifications).
- `N_seconds`: Acquisition window duration (seconds).
- `buffer_length`: Buffer size for circular buffer.

Data Flow:
1. The `acquire_data` function reads data from both thermocouple channels and saves it to buffers and a text file.
2. The `update_plot` function updates the real-time plot with the latest buffered data.
3. A background thread handles data acquisition, while the main thread updates the plot.
"""
import nidaqmx
from nidaqmx.constants import AcquisitionType, ThermocoupleType, TemperatureUnits
import threading
import time
import numpy as np
import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from collections import deque

# Acquisition parameters
Fs = 4  # Sample rate (in Hz) appropriate for thermocouples
N_seconds = 2  # Number of seconds per acquisition window
nbr_sample = int(Fs * N_seconds)
buffer_length = Fs * 4  # Set buffer size for circular buffer

# Thermocouple configuration
ni_board = "cDAQ1Mod1/ai0:1"  # Two channels: ai0 and ai1
thermocouple_type = ThermocoupleType.K  # Adjust if using a different type

# Data buffers
data_buffer_ch0 = deque(maxlen=buffer_length)
data_buffer_ch1 = deque(maxlen=buffer_length)

# Set up plot
plt.ion()
fig, ax = plt.subplots(figsize=(10, 6))
line_ch0, = ax.plot([], [], color="b", label="Channel 0")
line_ch1, = ax.plot([], [], color="r", label="Channel 1")
ax.set_xlabel('Time (s)')
ax.set_ylabel('Temperature (Â°C)')
ax.set_title('Thermocouple Temperature Data')

# File path for saving data
file_path = "thermocouple_data.txt"


# Data acquisition function
def acquire_data():
    with nidaqmx.Task() as task, open(file_path, "w") as f:
        # Configure thermocouple channels
        task.ai_channels.add_ai_thrmcpl_chan(
            ni_board,
            thermocouple_type=thermocouple_type,
            units=TemperatureUnits.DEG_C
        )
        task.timing.cfg_samp_clk_timing(rate=Fs, sample_mode=AcquisitionType.CONTINUOUS)

        # Write header to file
        f.write("Timestamp (s)\tChannel 0 (Â°C)\tChannel 1 (Â°C)\n")

        while True:
            # Read data from NI device (2 channels)
            data = task.read(number_of_samples_per_channel=nbr_sample)
            timestamps = np.linspace(time.time(), time.time() + N_seconds, nbr_sample)

            # Separate data by channel and add to respective buffers
            data_ch0 = np.array(data[0])
            data_ch1 = np.array(data[1])
            data_buffer_ch0.extend(data_ch0)
            data_buffer_ch1.extend(data_ch1)

            # Write data to file with timestamps
            for ts, temp0, temp1 in zip(timestamps, data_ch0, data_ch1):
                f.write(f"{ts:.3f}\t{temp0}\t{temp1}\n")

            f.flush()  # Ensure data is written to disk
            time.sleep(1)


# Update plot function
def update_plot(frame):
    # Convert deque to array and define time axis
    data_ch0 = np.array(data_buffer_ch0)
    data_ch1 = np.array(data_buffer_ch1)
    t = np.linspace(-N_seconds * len(data_ch0) / Fs, 0, len(data_ch0))

    # Update plot
    ax.clear()
    ax.plot(t, data_ch0, color="b", label="Channel 0")
    ax.plot(t, data_ch1, color="r", label="Channel 1")
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Temperature (Â°C)')
    ax.set_title('Thermocouple Temperature Data')
    ax.legend()
    plt.tight_layout()


# Start data acquisition in a separate thread
thread_acquisition = threading.Thread(target=acquire_data)
thread_acquisition.daemon = True
thread_acquisition.start()

# Animate plot
ani = FuncAnimation(fig, update_plot, interval=1000)
plt.show()