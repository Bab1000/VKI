import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from smt.surrogate_models import KRG
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression

# === 1. Load and preprocess the data ===
df = pd.read_csv("/home/jpe/VKI/Courses/PCMAF/Results/Same_dx/Res_tuning_Pstat_Stagline_dx_conv_only.csv")

# Extract input features (Pstat, Tinlet) and target variable (log10 of dx_conv)
X_raw = df[["Pstat", "Tinlet"]].values
y_raw = np.log10(df["dx_conv"].values)

# Optional: Mask out extreme or irrelevant values (disabled here)
# mask = (y_raw > -1.5) & (y_raw < -0.3)
# X_raw, y_raw = X_raw[mask], y_raw[mask]

# === 2. Train/test split ===
X_train, X_test, y_train, y_test = train_test_split(X_raw, y_raw, test_size=0.2, random_state=42)

# === 3. Normalize input features ===
scaler_X = StandardScaler()
X_train_scaled = scaler_X.fit_transform(X_train)
X_test_scaled = scaler_X.transform(X_test)

# === 4. Train the SMT Kriging model ===
# Matern 5/2 correlation, constant trend, small noise level for numerical stability
kriging_model = KRG(
    print_global=True,
    corr="matern52",
    poly="constant",
    theta0=[1.0]      # initial guess for theta 
)

kriging_model.set_training_values(X_train_scaled, y_train)
kriging_model.train()

# === 5. Predict on test set ===
y_pred = kriging_model.predict_values(X_test_scaled).flatten()

# === 6. Evaluation metrics ===
mae = mean_absolute_error(y_test, y_pred)
rmse = np.sqrt(mean_squared_error(y_test, y_pred))
r2 = r2_score(y_test, y_pred)

print("\n SMT-Kriging results (log10(dx_conv)):")
print(f"→ MAE   = {mae:.4f}")
print(f"→ RMSE  = {rmse:.4f}")
print(f"→ R²    = {r2:.4f}")

# === 7. Baseline: linear regression for comparison ===
lr_model = LinearRegression()
lr_model.fit(X_train_scaled, y_train)
y_lr_pred = lr_model.predict(X_test_scaled)
r2_lr = r2_score(y_test, y_lr_pred)

print("\n Linear Regression baseline:")
print(f"→ R² (linear) = {r2_lr:.4f}")

# === 8. Scatter plot: True vs Predicted ===
plt.figure(figsize=(6, 6))
plt.scatter(y_test, y_pred, alpha=0.8, label="Kriging", color='blue')
plt.scatter(y_test, y_lr_pred, alpha=0.4, label="Linear Regression", color='orange')
plt.plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], 'r--', lw=2)
plt.xlabel("True log10(dx_conv)")
plt.ylabel("Predicted")
plt.title("Prediction Comparison")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("Performance.png", dpi=300)
plt.show()



# === 6b. Normalized RMSE ===
range_y = y_test.max() - y_test.min()
mean_abs_y = np.mean(np.abs(y_test))

nrmse_range = (rmse / range_y) * 100
nrmse_mean = (rmse / mean_abs_y) * 100

print(f"→ nRMSE (range-based) = {nrmse_range:.2f}%")
print(f"→ nRMSE (mean-based)  = {nrmse_mean:.2f}%")

# === 10. Kriging surface vs real data (3D plot) ===
from mpl_toolkits.mplot3d import Axes3D

# Create grid of inputs in normalized space
pstat_range = np.linspace(X_train[:, 0].min(), X_train[:, 0].max(), 50)
tinlet_range = np.linspace(X_train[:, 1].min(), X_train[:, 1].max(), 50)
PSTAT, TINLET = np.meshgrid(pstat_range, tinlet_range)

# Flatten and normalize the input grid
X_grid = np.column_stack([PSTAT.ravel(), TINLET.ravel()])
X_grid_scaled = scaler_X.transform(X_grid)

# Predict values on the grid
Y_pred_grid = kriging_model.predict_values(X_grid_scaled).reshape(PSTAT.shape)

# Create 3D plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Scatter true training points
ax.scatter(X_train[:, 0], X_train[:, 1], y_train, color='black', label='Training data')

# Surface from Kriging
ax.plot_surface(PSTAT, TINLET, Y_pred_grid, cmap='viridis', alpha=0.6, edgecolor='none', label='Kriging Surface')

ax.set_xlabel("Pstat [mbar]")
ax.set_ylabel("Tinlet [K]")
ax.set_zlabel("log10(dx_conv)")
ax.set_title("Kriging Prediction Surface vs Training Data")
plt.tight_layout()
plt.savefig("kriging_surface_vs_data.png", dpi=300)
plt.show()


# === 11. Linear Regression surface vs real data (3D plot) ===

# Predict values from the linear regression model on the same grid
Y_lr_grid = lr_model.predict(X_grid_scaled).reshape(PSTAT.shape)

# Create 3D plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Scatter true training points
ax.scatter(X_train[:, 0], X_train[:, 1], y_train, color='black', label='Training data')

# Surface from Linear Regression
ax.plot_surface(PSTAT, TINLET, Y_lr_grid, cmap='plasma', alpha=0.6, edgecolor='none', label='Linear Surface')

ax.set_xlabel("Pstat [mbar]")
ax.set_ylabel("Tinlet [K]")
ax.set_zlabel("log10(dx_conv)")
ax.set_title("Linear Regression Surface vs Training Data")
plt.tight_layout()
plt.savefig("linear_regression_surface_vs_data.png", dpi=300)
plt.show()





