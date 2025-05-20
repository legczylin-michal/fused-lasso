import numpy as np
import matplotlib.pyplot as plt

from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Lasso

np.random.seed(666)

class LinRegLasso:
    __slots__ = ["__beta_tilde", "__gamma", "__max_iter", "__tol", "__mean_X", "__std_X", "__mean_y"]

    @property
    def coef_(self):
        return self.__beta_tilde

    @property
    def intercept_(self):
        return self.__mean_y
    
    @staticmethod
    def soft_threshold(a, b):
        return np.sign(a) * np.maximum((np.abs(a) - b), 0)
    
    def __init__(self, C=1.0, max_iter=1000, tol=0.0001):
        self.__beta_tilde = np.array([])
        
        self.__gamma = float("+inf") if C == 0 else 1 / C
        self.__max_iter = max_iter
        self.__tol = tol

        self.__mean_X = 0
        self.__std_X = 1
        self.__mean_y = 0
        
        return

    def fit(self, X_standarized, y_centered):
        # n, p = X.shape
        
        # self.__mean_X = X.mean(axis=0)
        # self.__std_X = X.std(axis=0)

        # X_standarized = (X - self.__mean_X) / (np.sqrt(n - 1) * self.__std_X)

        # self.__mean_y = y.mean()
        # y_centered = y - self.__mean_y

        self.__beta_tilde = LinRegLasso.soft_threshold(X_standarized.T @ y_centered, self.__gamma)
        
        converged = False
        for _ in range(self.__max_iter):
            old_beta_tilde = self.__beta_tilde.copy()
            
            self.__beta_tilde = LinRegLasso.soft_threshold(self.__beta_tilde + X_standarized.T @ (y_centered - X_standarized @ self.__beta_tilde), self.__gamma)

            if np.linalg.norm(old_beta_tilde - self.__beta_tilde, 2) < self.__tol:
                converged = True
                
                break

        if not converged:
            print("algorithm did not converge")
        
        return self

    def predict(self, X_standarized):
        # n, p = X.shape
        
        # X_standarized = (X - self.__mean_X) / (np.sqrt(n - 1) * self.__std_X)
        
        return X_standarized @ self.__beta_tilde #+ self.__mean_y
    
class LinRegFusedLasso:
    __slots__ = ["__beta_tilde", "__lambda_1", "__lambda_2", "__max_iter", "__tol", "__mean_X", "__std_X", "__mean_y"]

    @property
    def coef_(self):
        return self.__beta_tilde

    @property
    def intercept_(self):
        return self.__mean_y
    
    @staticmethod
    def help_solve(d, a, b, c, lambda_a, lambda_b, lambda_c):
        if d + lambda_a + lambda_b + lambda_c <= a:
            return d + lambda_a + lambda_b + lambda_c
        elif a < d - lambda_a + lambda_b + lambda_c <= b:
            return d - lambda_a + lambda_b + lambda_c
        elif b < d - lambda_a - lambda_b + lambda_c <= c:
            return d - lambda_a - lambda_b + lambda_c
        elif c < d - lambda_a - lambda_b - lambda_c:
            return d - lambda_a - lambda_b - lambda_c

        return None

    def __init__(self, lambda_1=1.0, lambda_2=1.0, max_iter=1000, tol=0.0001):
        self.__beta_tilde = None
        
        self.__lambda_1 = lambda_1
        self.__lambda_2 = lambda_2
        self.__max_iter = max_iter
        self.__tol = tol

        self.__mean_X = 0
        self.__std_X = 1
        self.__mean_y = 0
        
        return

    def __loss(self, beta, X, y):
        return 0.5 * ((y - X @ beta) ** 2).sum() + self.__lambda_1 * abs(beta).sum() + self.__lambda_2 * abs(beta[1:] - beta[:-1]).sum()

    def fit(self, X_standarized, y_centered):
        # n, p = X.shape
        
        # self.__mean_X = X.mean(axis=0)
        # self.__std_X = X.std(axis=0)

        # X_standarized = (X - self.__mean_X) / (np.sqrt(n - 1) * self.__std_X)

        # self.__mean_y = y.mean()
        # y_centered = y - self.__mean_y

        n, p = X_standarized.shape

        self.__beta_tilde = np.zeros(p)

        converged = False
        for _ in range(self.__max_iter):
            old_beta_tilde = self.__beta_tilde.copy()

            d = X_standarized.T @ (y_centered - X_standarized @ self.__beta_tilde) + self.__beta_tilde

            for j in range(p):                    
                breakpoints = np.array([0, 0 if j == 0 else self.__beta_tilde[j - 1], 0 if j == p - 1 else self.__beta_tilde[j + 1]])
                penalties = np.array([self.__lambda_1, 0 if j == 0 else self.__lambda_2, 0 if j == p - 1 else self.__lambda_2])
                s = np.argsort(breakpoints)
                breakpoints = breakpoints[s]
                penalties = penalties[s]

                a, b, c = breakpoints[0], breakpoints[1], breakpoints[2]
                lambda_a, lambda_b, lambda_c = penalties[0], penalties[1], penalties[2]

                sol = LinRegFusedLasso.help_solve(d[j], a, b, c, lambda_a, lambda_b, lambda_c)

                min_v = float("+inf")
                if sol is None:
                    for _breakpoint in breakpoints:
                        self.__beta_tilde[j] = _breakpoint
    
                        v = self.__loss(self.__beta_tilde, X_standarized, y_centered)
                        if v < min_v:
                            min_v = v
                            sol = _breakpoint
                    
                self.__beta_tilde[j] = sol

            if np.linalg.norm(old_beta_tilde - self.__beta_tilde, 2) < self.__tol:
                converged = True
                
                break

        if not converged:
            print("algorithm did not converge")
        
        return self

    def predict(self, X_standarized):
        # n, p = X.shape
        
        # X_standarized = (X - self.__mean_X) / (np.sqrt(n - 1) * self.__std_X)
        
        return X_standarized @ self.__beta_tilde #+ self.__mean_y
    
