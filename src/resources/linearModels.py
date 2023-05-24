import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy.stats import chi2



class LinearRegression:
    """Authored by professor Emanuel Gonçalves, this call allows fo the computation of the linear regression coefs for a set of features. 
    Do check if the response variable is Normally distributed, since this is nor a GLM
    """
    
    

    def __init__(
        self,
        Y,
        X,
        M,
        M2=None,
        normalize=False,
        fit_intercept=True,
        copy_X=True,
        n_jobs=4,
        verbose=1,
    ):
        self.samples = set.intersection(
            set(Y.index),
            set(X.index),
            set(M.index),
            set(Y.index) if M2 is None else set(M2.index),
        )

        self.X = X.loc[self.samples]
        self.X = self.X.loc[:, self.X.count() > (M.shape[1] + (1 if M2 is None else 2))]
        self.X_ma = np.ma.masked_invalid(self.X.values)

        self.Y = Y.loc[self.samples]
        self.Y = self.Y.loc[:, self.Y.std() > 0]

        self.M = M.loc[self.samples]

        self.M2 = M2.loc[self.samples, self.X.columns] if M2 is not None else M2

        self.normalize = normalize
        self.fit_intercept = fit_intercept
        self.copy_X = copy_X
        self.n_jobs = n_jobs

        self.verbose = verbose
        # self.log = logging.getLogger("Crispy")

    def model_regressor(self):
        regressor = LinearRegression(
            fit_intercept=self.fit_intercept,
            normalize=self.normalize,
            copy_X=self.copy_X,
            n_jobs=self.n_jobs,
        )
        return regressor

    @staticmethod
    def loglike(y_true, y_pred):
        nobs = len(y_true)
        nobs2 = nobs / 2.0

        ssr = np.power(y_true - y_pred, 2).sum()

        llf = -nobs2 * np.log(2 * np.pi) - nobs2 * np.log(ssr / nobs) - nobs2

        return llf

    @staticmethod
    def multipletests_per(
        associations, method="fdr_bh", field="pval", fdr_field="fdr", index_cols=None
    ):
        index_cols = ["y_id"] if index_cols is None else index_cols

        d_unique = {tuple(i) for i in associations[index_cols].values}

        df = associations.set_index(index_cols)

        df = pd.concat(
            [
                df.loc[i]
                .assign(fdr=multipletests(df.loc[i, field], method=method)[1])
                .rename(columns={"fdr": fdr_field})
                for i in d_unique
            ]
        ).reset_index()

        return df

    def fit_matrix(self):
        lms = []

        for x_idx, x_var in enumerate(self.X):
            if self.verbose > 0:
                self.log.info(f"LM={x_var} ({x_idx})")

            # Mask NaNs
            x_ma = np.ma.mask_rowcols(self.X_ma[:, [x_idx]], axis=0)

            # Build matrices
            x = self.X.iloc[~x_ma.mask.any(axis=1), [x_idx]]
            y = self.Y.iloc[~x_ma.mask.any(axis=1), :]

            # Covariate matrix (remove invariable features and add noise)
            m = self.M.iloc[~x_ma.mask.any(axis=1), :]
            if self.M2 is not None:
                m2 = self.M2.iloc[~x_ma.mask.any(axis=1), [x_idx]]
                m = pd.concat([m2, m], axis=1)
            m = m.loc[:, m.std() > 0]
            m += np.random.normal(0, 1e-4, m.shape)

            # Fit covariate model
            lm_small = self.model_regressor().fit(m, y)
            lm_small_ll = self.loglike(y, lm_small.predict(m))

            # Fit full model: covariates + feature
            lm_full_x = np.concatenate([m, x], axis=1)
            lm_full = self.model_regressor().fit(lm_full_x, y)
            lm_full_ll = self.loglike(y, lm_full.predict(lm_full_x))

            # Log-ratio test
            lr = 2 * (lm_full_ll - lm_small_ll)
            lr_pval = chi2(1).sf(lr)

            # Assemble + append results
            res = pd.DataFrame(
                dict(
                    y_id=y.columns,
                    x_id=x_var,
                    n=y.attrs["nan_mask"].loc[y.columns, x.index].sum(1) if "nan_mask" in y.attrs else len(x),
                    beta=lm_full.coef_[:, -1],
                    lr=lr.values,
                    covs=m.shape[1],
                    pval=lr_pval,
                    fdr=multipletests(lr_pval, method="fdr_bh")[1],
                )
            )

            lms.append(res)

        lms = pd.concat(lms, ignore_index=True).sort_values("pval")

        return lms

    @staticmethod
    def lm_residuals(y, x, fit_intercept=True, add_intercept=False):
        # Prepare input matrices
        ys = y.dropna()

        xs = x.loc[ys.index].dropna()
        xs = xs.loc[:, xs.std() > 0]

        ys = ys.loc[xs.index]

        if ys.shape[0] <= xs.shape[1]:
            return None

        # Linear regression models
        lm = LinearRegression(fit_intercept=fit_intercept).fit(xs, ys)

        # Calculate residuals
        residuals = ys - lm.predict(xs) - lm.intercept_

        # Add intercept
        if add_intercept:
            residuals += lm.intercept_

        return residuals
    

