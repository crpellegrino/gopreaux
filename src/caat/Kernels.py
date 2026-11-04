import warnings

import numpy as np
from sklearn.gaussian_process.kernels import Kernel as SklearnKernel
from sklearn.gaussian_process.kernels import KernelOperator

warnings.filterwarnings("ignore")


class Kernel(SklearnKernel):
    """
    A custom Kernel class used to store kernel hyperparameters
    and use them in the sklearn GaussianProcessRegression class.
    This is intended to be used to produce linear combinations of
    the basic Sklearn kernels.

    This class inherits the sklearn `Kernel` object and overwrites
    some of its methods to maintain compatibility with the base class
    and the `GaussianProcessRegressor` class while enabling users to
    input custom kernels.
    """
    def __init__(self, kernel: SklearnKernel):
        """
        Initialize a Kernel object with an sklearn `Kernel`.

        Args:
            kernel (SklearnKernel): The sklearn `Kernel` object.
        """
        self.kernel = kernel

    def __call__(self, X, Y=None, eval_gradient=False):
        return self.kernel.__call__(X, Y=Y, eval_gradient=eval_gradient)

    def is_stationary(self):
        return self.kernel.is_stationary()

    def diag(self, X):
        return self.kernel.diag(X)

    def _recursively_get_component(self, kernel=None):
        """
        Recursively get all components of self.kernel
        
        Args:
            kernel (SklearnKernel | Kernel | None, optional): The 
                sklearn `Kernel` object.
        """
        if kernel is None:
            kernel = self.kernel

        if not isinstance(kernel, KernelOperator):
            return np.asarray([kernel])

        return np.hstack(
            [
                self._recursively_get_component(kernel.k1),
                self._recursively_get_component(kernel.k2),
            ]
        )

    @property
    def components(self):
        return self._recursively_get_component()

    def recursively_set_params(
        self,
        values: list,
        bounds: list | str,
        kernel=None,
    ):
        """
        Recursively set parameters for all component kernels.
        Gets the valid param keywords for the kernel and uses
        them to set the values of the correct parameters.
        If Kernel is a product or sum of kernels, assumes the
        input values and bounds are in the same order as the kernels.

        NOTE: `values` are the hyperparameter values themselves (e.g., length
        scales), not the log-transformed values returned by the sklearn
        `theta` attribute. To set hyperparameters from `theta`, use
        `fix_kernel_from_theta`.

        Args:
            values (list): The kernel parameters
            bounds (list | str): The kernel parameter bounds
            kernel (_type_, optional): The kernel to set parameters
                of. Defaults to None.
        """
        if kernel is None:
            kernel = self.kernel

        if not isinstance(kernel, KernelOperator):
            valid_param_keys = list(kernel.get_params().keys())
            valid_params = {valid_param_keys[0]: values}
            valid_params[valid_param_keys[1]] = bounds
            kernel.set_params(**valid_params)

        else:
            k1_dim = kernel.k1.n_dims
            self.recursively_set_params(
                values[:k1_dim][0] if k1_dim == 1 else values[:k1_dim],
                bounds if isinstance(bounds, str) else bounds[:k1_dim],
                kernel=kernel.k1,
            )
            self.recursively_set_params(
                values[k1_dim:][0] if kernel.k2.n_dims == 1 else values[k1_dim:],
                bounds if isinstance(bounds, str) else bounds[k1_dim:],
                kernel=kernel.k2,
            )


def fix_kernel_from_theta(kernel: Kernel | SklearnKernel, theta: np.ndarray):
    """
    Clone a kernel with its free hyperparameters set from `theta`, and fix them.

    `theta` is the log-transformed vector of free (non-fixed) hyperparameters that
    sklearn stores in the `theta` attribute of a kernel, e.g. the `kernel_.theta` of
    a fitted `GaussianProcessRegressor`. It is applied through sklearn's own
    `clone_with_theta`, so it works for any single or compound (sum / product)
    kernel, as long as `theta` has the same layout as `kernel.theta`.

    Args:
        kernel (Kernel | SklearnKernel): The kernel to clone. It is not modified.
        theta (np.ndarray): The log-transformed values of the free hyperparameters.

    Returns:
        Kernel | SklearnKernel: A new kernel of the same type as `kernel`, with the
            hyperparameters set from `theta` and the bounds of every free
            hyperparameter set to "fixed".
    """
    inner_kernel = kernel.kernel if isinstance(kernel, Kernel) else kernel
    fixed_kernel = inner_kernel.clone_with_theta(np.asarray(theta))
    for hyperparameter in fixed_kernel.hyperparameters:
        if not hyperparameter.fixed:
            fixed_kernel.set_params(**{f"{hyperparameter.name}_bounds": "fixed"})

    return Kernel(fixed_kernel) if isinstance(kernel, Kernel) else fixed_kernel
