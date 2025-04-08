import pymc as pm
import numpy as np
import matplotlib.pyplot as plt
import arviz as az

# Generate synthetic data
np.random.seed(42)
x = np.linspace(0, 10, 50)
true_a = 2.5
true_b = 1.0
sigma = 1.0
y = true_a * x + true_b + np.random.normal(0, sigma, size=len(x))

# Bayesian model
with pm.Model() as model:
    # Priors for a and b
    a = pm.Normal("a", mu=0, sigma=10)
    b = pm.Normal("b", mu=0, sigma=10)
    # Noise std (can be known or estimated)
    sigma_obs = pm.HalfNormal("sigma", sigma=1)
    
    # Forward model
    mu = a * x + b
    
    # Likelihood
    y_obs = pm.Normal("y_obs", mu=mu, sigma=sigma_obs, observed=y)
    
    # Inference
    trace = pm.sample(2000, tune=1000, return_inferencedata=True)

# Plot posterior distributions
az.plot_trace(trace)
plt.show()

# Summary
print(az.summary(trace, hdi_prob=0.95))
