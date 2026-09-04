import importlib
_TS_UTILS_ATTRS = {
    "get_bottom_layer_indices": ".ts_utils",
    "freeze_bottom_layers": ".ts_utils",
    "set_seed": ".ts_utils",
    "load_state": ".ts_utils",
    "build_temperature_tensor": ".ts_utils",
    "get_d3_xc_params": ".ts_utils",
    "build_d3_params": ".ts_utils",
    "build_maced3_model": ".ts_utils",
    "make_weight_reporter": ".ts_utils",
    "make_reporter": ".ts_utils",
    "run_langevin": ".ts_utils",
    "run_sspd": ".ts_utils",
    "run_minimization": ".ts_utils",
    "unwrap_positions": ".ts_convert",
    "load_h5md": ".ts_convert",
    "convert_h5_xyz": ".ts_convert",
}

__all__ = list(_TS_UTILS_ATTRS)


def __getattr__(name):
    try:
        submodule_name = _TS_UTILS_ATTRS[name]
    except KeyError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
    module = importlib.import_module(submodule_name, __name__)
    value = getattr(module, name)
    globals()[name] = value  # cache on the package so later lookups skip __getattr__
    return value


def __dir__():
    return sorted(list(globals()) + __all__)
