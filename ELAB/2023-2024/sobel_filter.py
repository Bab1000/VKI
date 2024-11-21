import numpy as np
import os
import csv
import cv2
import fig_management
import matplotlib.pyplot as plt

# Specify the path to the folder containing your CSV file
data_folder = 'Natural Convection/'

# Specify the filename of your CSV file
csv_filename = 'Natural_conv_26V_3A.csv'

# Construct the full path to the CSV file
csv_filepath = os.path.join(data_folder, csv_filename)

# Load the data from the CSV file using the csv module
with open(csv_filepath, 'r') as file:
    # Assuming your CSV file has no header and contains numerical values separated by commas
    reader = csv.reader(file, delimiter=';')
    matrix = [[float(num) for num in row] for row in reader]

# Convert the matrix to a NumPy array
image = np.array(matrix)

# Perform edge detection (replace with your actual edge detection code)
# Example using Sobel operators
sobel_horizontal = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]])
sobel_vertical = np.array([[-1, -2, -1], [0, 0, 0], [1, 2, 1]])
gradient_x = cv2.filter2D(image, -1, sobel_horizontal)
gradient_y = cv2.filter2D(image, -1, sobel_vertical)
gradient_magnitude = np.sqrt(gradient_x**2 + gradient_y**2)
threshold = 1100.0
edges = (gradient_magnitude > threshold).astype(np.uint8)

# Find contours in the edge-detected image
contours, _ = cv2.findContours(edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

# Assuming the largest contour corresponds to the circular region
largest_contour = max(contours, key=cv2.contourArea)

# Create a mask for the largest contour
mask = np.zeros_like(image, dtype=np.uint8)
cv2.drawContours(mask, [largest_contour], -1, 255, thickness=cv2.FILLED)

# Apply the mask to the original image
values_inside_circle = image * (mask > 0)

"""
# Display the original matrix
plt.subplot(121)
plt.imshow(image, cmap='gray')
plt.title('Original Matrix')

# Display the edge detection result
plt.subplot(122)
plt.imshow(edges, cmap='gray')
plt.title('Edge Detection Result')

# Display the mask for the largest contour
plt.subplot(133)
plt.imshow(mask, cmap='gray')
plt.title('Largest Contour Mask')

plt.show()
"""
vmin1, vmax1 = 3000,8000

# Display the values inside the circular region
plt.figure()
plt.plot()
im1 = plt.imshow(values_inside_circle, cmap='jet', vmin=vmin1, vmax=vmax1)
plt.colorbar(im1, label='IU', fraction=0.046, pad=0.04)
plt.xlabel("Horizontal position [pixels]")
plt.ylabel("Vertical position [pixels]")
plt.savefig("heatmap.pdf")

plt.figure()
plt.plot()
im2 = plt.imshow(image, cmap='jet', vmin=vmin1, vmax=vmax1)
plt.colorbar(im2, label='IU', fraction=0.046, pad=0.04)
plt.xlabel("Horizontal position [pixels]")
plt.ylabel("Vertical position [pixels]")
plt.savefig("heatmap_edge.pdf")

plt.figure()
plt.plot()
im2 = plt.imshow(values_inside_circle[10:41,150:201], cmap='jet', vmin=vmin1, vmax=vmax1, extent=[150,200,40,10])
plt.colorbar(im2, label='IU', fraction=0.046, pad=0.04)
plt.xlabel("Horizontal position [pixels]")
plt.ylabel("Vertical position [pixels]")
plt.savefig("heatmap_edge_zoom.pdf")
