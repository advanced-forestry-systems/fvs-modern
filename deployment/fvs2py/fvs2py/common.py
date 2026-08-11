from __future__ import annotations

import ctypes as ct
import functools
import logging
import os
from collections.abc import Callable
from typing import TypeVar, cast

try:  # ParamSpec moved into typing in 3.10; fall back for 3.9 runtimes (e.g. OSC default)
    from typing import ParamSpec
except ImportError:
    from typing_extensions import ParamSpec

P = ParamSpec("P")
T = TypeVar("T")
ClassT = TypeVar("ClassT", bound=type)


def load_dll(dll_path: str | os.PathLike) -> ct.CDLL:
    """Loads a Dynamic Link Library.

    Supports linux shared objects, not Windows DLLs.

    Args:
        dll_path (str | os.PathLike): the path to the DLL
            to be loaded.

    Returns:
        the loaded DLL as a ctypes CDLL instance
    """
    dll_path_str = os.fspath(dll_path)
    logging.info(f"Loading library from {dll_path_str}")
    dll = ct.CDLL(dll_path_str, mode=1)  # RTLD_LAZY: defer undefined symbol resolution
    logging.info(f"Library loaded successfully from {dll_path_str}")
    return dll


def unload_dll(dll: ct.CDLL) -> None:
    """Unloads a Dynamic Link Library.

    Args:
        dll: the loaded DLL.
    """
    close_func = dll.dlclose
    close_func.argtypes = (ct.c_void_p,)
    close_func.restype = ct.c_int
    result = close_func(dll._handle)
    if result != 0:
        msg = f"Failed to unload DLL: {result}"
        raise RuntimeError(msg)
    logging.debug("Library unloaded successfully.")


def function_requires_fvs_library() -> Callable[
    [Callable[P, T]], Callable[P, T]
]:
    """Decorator ensuring a ._lib attribute is present and not None.

    This decorator is intended to be applied to instance methods that rely upon
    accessing attributes or routines frome the FVS library. If the FVS library has been
    unloaded without this protection, a segmentation fault would occur.
    """

    def decorator(func: Callable[P, T]) -> Callable[P, T]:
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs) -> T:
            if not hasattr(self, "_lib") or self._lib is None:
                msg = f"FVS library not loaded, unable to call {func.__name__}."
                raise RuntimeError(msg)
            return func(self, *args, **kwargs)

        return cast(Callable[P, T], wrapper)

    return decorator


def fvs_property(func: Callable[..., T]) -> property:
    """A property decorator that checks if the FVS library is loaded before accessing.

    This combines @property with a check for self._lib. This decorator will help ensure
    that the FVS library is loaded before accessing the property. If the FVS library
    has been unloaded without this protection, a segmentation fault would occur.

    Args:
        func: The property getter function to wrap.

    Returns:
        A property descriptor with FVS library checking.
    """

    @functools.wraps(func)
    def wrapper(self) -> T:
        if not hasattr(self, "_lib") or self._lib is None:
            msg = f"FVS library not loaded, unable to access {func.__name__} property."
            raise RuntimeError(msg)
        return func(self)

    return property(wrapper)


def class_requires_fvs_library(cls: ClassT) -> ClassT:
    """Decorator enforcing `function_requires_fvs_library` for all public methods.

    This decorator applies the `function_requires_fvs_library` decorator to all public
    methods of the class.
    """
    for name, method in cls.__dict__.items():
        if callable(method) and not name.startswith("_"):
            decorated_method = function_requires_fvs_library()(method)
            setattr(cls, name, decorated_method)
    return cls
