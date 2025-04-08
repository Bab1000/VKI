import pymc as pm
import numpy as np
import matplotlib.pyplot as plt
import arviz as az

class BayesianLinearModel:
    def __init__(self, x, y, name="LinearModel"):
        self.x = x
        self.y = y
        self.name = name
        self.model = None
        self.trace = None
        self.posterior_predictive = None

    def build_model(self):
        with pm.Model() as model:
            # Priors
            a = pm.Normal("a", mu=0, sigma=10)
            b = pm.Normal("b", mu=0, sigma=10)
            sigma = pm.HalfNormal("sigma", sigma=1)

            # Forward model
            mu = a * self.x + b

            # Likelihood
            y_obs = pm.Normal("y_obs", mu=mu, sigma=sigma, observed=self.y)

            self.model = model

    def fit(self, draws=2000, tune=1000):
        if self.model is None:
            self.build_model()
        with self.model:
            self.trace = pm.sample(draws=draws, tune=tune, return_inferencedata=True)

    def predict(self, x_new, n_samples=1000):
        if self.trace is None:
            raise ValueError("Model not fitted. Call `.fit()` first.")

        # Extraire les échantillons du trace
        a_samples = self.trace.posterior["a"].stack(draws=("chain", "draw")).values
        b_samples = self.trace.posterior["b"].stack(draws=("chain", "draw")).values
        sigma_samples = self.trace.posterior["sigma"].stack(draws=("chain", "draw")).values

        idx = np.random.choice(len(a_samples), size=n_samples, replace=False)
        a_s = a_samples[idx]
        b_s = b_samples[idx]
        sigma_s = sigma_samples[idx]

        preds = np.array([
            a * x_new + b + np.random.normal(0, sigma)
            for a, b, sigma in zip(a_s, b_s, sigma_s)
        ])

        return preds

    def plot_trace(self):
        if self.trace is not None:
            az.plot_trace(self.trace)
            plt.show()
        else:
            print("⚠️ Trace is empty. Run .fit() first.")

    def summary(self, hdi_prob=0.95):
        if self.trace is not None:
            return az.summary(self.trace, hdi_prob=hdi_prob)
        else:
            print("⚠️ No trace available. Run .fit() first.")
            return None

    def plot_predictions(self, x_grid=None):
        if x_grid is None:
            x_grid = np.linspace(self.x.min(), self.x.max(), 100)

        preds = self.predict(x_grid)
        mean_preds = preds.mean(axis=0)
        hdi_lower = np.percentile(preds, 2.5, axis=0)
        hdi_upper = np.percentile(preds, 97.5, axis=0)

        plt.figure(figsize=(10, 6))
        plt.scatter(self.x, self.y, label="Observed data")
        plt.plot(x_grid, mean_preds, color='C1', label="Mean prediction")
        plt.fill_between(x_grid, hdi_lower, hdi_upper, alpha=0.3, label="95% HDI")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Posterior Predictive Plot")
        plt.legend()
        plt.grid(True)
        plt.show()

# Exemple d'utilisation
if __name__ == "__main__":
    # Données synthétiques
    np.random.seed(42)
    x = np.linspace(0, 10, 50)
    true_a = 2.5
    true_b = 1.0
    sigma = 1.0
    y = true_a * x + true_b + np.random.normal(0, sigma, size=len(x))

    # Modèle bayésien
    model = BayesianLinearModel(x, y)
    model.fit()
    model.plot_trace()
    print(model.summary())

    # Prédictions
    model.plot_predictions()
