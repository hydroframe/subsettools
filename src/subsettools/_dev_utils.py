"""This module contains utilities for subsettools developers."""

import warnings
import functools


def replace_kwargs(kwarg_map):
    """Deprecate and replace function kwargs without changing the function."""

    def decorator(func):
        @functools.wraps(func)
        def wrapped(*args, **kwargs):
            new_kwargs = {}
            for k, v in kwargs.items():
                if k in kwarg_map:
                    warnings.warn(
                        f"keyword argument '{k}' is no longer valid. "
                        f"Use '{kwarg_map[k]}' instead.",
                        DeprecationWarning,
                        stacklevel=2,
                    )
                new_kwargs[kwarg_map.get(k, k)] = v
            return func(*args, **new_kwargs)

        return wrapped

    return decorator
