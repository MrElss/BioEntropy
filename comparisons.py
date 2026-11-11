from contextlib import contextmanager
from copy import deepcopy
import hashlib
import types
import numpy as np
import comparison_statistics as main_statistics
import visit_statistics as raw


def _code_key(value):
    if isinstance(value, types.CodeType):
        return (
            value.co_code,
            value.co_names,
            value.co_varnames,
            value.co_freevars,
            value.co_cellvars,
            value.co_argcount,
            value.co_kwonlyargcount,
            tuple(_code_key(x) for x in value.co_consts),
        )
    if isinstance(value, (tuple, list)):
        return tuple(_code_key(x) for x in value)
    return value


def _array_key(value):
    a = np.ascontiguousarray(value)
    return (a.shape, a.dtype.str, hashlib.blake2b(a.view(np.uint8), digest_size=20).digest())


def _value_key(value):
    if isinstance(value, np.ndarray):
        return ("array", _array_key(value))
    if isinstance(value, (tuple, list)):
        return ("list", tuple(_value_key(x) for x in value))
    if isinstance(value, dict):
        return ("dict", tuple(sorted((k, _value_key(v)) for k, v in value.items())))
    return ("value", value)


@contextmanager
def shared_evaluations(enabled=True):
    cache = {}
    stats = dict(evaluations=0, reused=0)
    original = []
    if not enabled:
        yield stats
        return

    def wrap(function, category):
        semantic = _code_key(function.__code__)

        def call(*args, **kwargs):
            key = (category, semantic, tuple(_value_key(x) for x in args), _value_key(kwargs))
            if key not in cache:
                cache[key] = deepcopy(function(*args, **kwargs))
                stats["evaluations"] += 1
            else:
                stats["reused"] += 1
            return deepcopy(cache[key])

        return call

    try:
        for module in (main_statistics, raw):
            for name in ("entropy_labelings", "energy_labelings"):
                function = getattr(module, name)
                original.append((module, name, function))
                setattr(module, name, wrap(function, name))
        yield stats
    finally:
        for module, name, function in original:
            setattr(module, name, function)


def all_comparisons(result):
    return result.get("comparisons_all", result.get("comparisons", []))


def _stamp(rows, source, mode):
    for row in rows:
        row["Software_version"] = "V1.0.2"
        row["Inference_source"] = source
        if mode == "animal":
            scope = "动物图01–05；原始值比较"
        elif source == "分周原值":
            scope = "临床图10–11；当周原始值比较"
        elif row["Contrast"] == "change":
            scope = "临床图02、04、05；相对0W变化之差"
        else:
            scope = "临床图01、03；当周原始值比较"
        row["Figure_scope"] = scope
    return rows


def compare_policies(result, permutations=9999, progress=print):
    dual = result["mode"] == "clinical" and result.get("selected_method", "Shannon") != "none"
    with shared_evaluations(dual) as stats:
        progress(result["name"] + "：计算图01–05的既定检验方案")
        first = main_statistics.compare(result, permutations, progress)
        first_settings = deepcopy(result["inference_settings"])
        first_numeric = deepcopy(result.get("permutation_numeric_checks", []))
        second, second_settings, second_numeric = ([], {}, [])
        if dual:
            progress(result["name"] + "：计算图10–11的当周原值检验方案")
            second = raw.compare(result, permutations, progress)
            second_settings = deepcopy(result["inference_settings"])
            second_numeric = deepcopy(result.get("permutation_numeric_checks", []))
    first = _stamp(first, "常规", result["mode"])
    second = _stamp(second, "分周原值", result["mode"])
    result["comparisons"] = first
    result["comparisons_visit"] = second
    result["comparisons_all"] = first + second
    result["permutation_numeric_checks"] = first_numeric
    result["permutation_numeric_checks_visit"] = second_numeric
    result["inference_settings_main"] = first_settings
    result["inference_settings_visit"] = second_settings
    result["inference_settings"] = dict(
        version="V1.0.2",
        permutations=permutations,
        policies={"常规": first_settings, **({"分周原值": second_settings} if dual else {})},
        matching_rule="Metric + Week + Contrast + Inference_source",
        figure_policy="animal 01–05,07,08: 常规; clinical 01–05: 常规; clinical 10–11: 分周原值",
        shared_kernel_evaluations=stats,
        status="完成" if permutations else "用户关闭",
    )
    return result["comparisons_all"]
