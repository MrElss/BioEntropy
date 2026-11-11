METHODS = {
    "Shannon熵（默认）": "Shannon",
    "Rényi二阶熵": "Renyi2",
    "能量距离（正常参照）": "Energy",
    "不计算（仅NMD）": "none",
}
FIELDS = {"Shannon": "Entropy_Shannon", "Renyi2": "Entropy_Renyi2", "Energy": "Energy_distance"}
LABELS = {
    "NMD": "NMD（正常参照马氏距离）",
    "Entropy_Shannon": "Shannon熵",
    "Entropy_Renyi2": "Rényi二阶熵",
    "Energy_distance": "能量距离",
}


def validate_method(method):
    if method not in (*FIELDS, "none"):
        raise ValueError("请选择Shannon、Renyi2、Energy或none。")
    return method


def metrics(result):
    method = result.get("selected_method")
    if method is not None:
        validate_method(method)
        return ["NMD"] + ([FIELDS[method]] if method != "none" else [])
    columns = set(result["summary"][0]) if result.get("summary") else set()
    return ["NMD"] + [field for field in FIELDS.values() if field in columns]


def method_name(metric):
    return "Energy" if metric == "Energy_distance" else metric.removeprefix("Entropy_")


def score_field(metric):
    return {
        "NMD": "NMD",
        "Entropy_Shannon": "Shannon_information_score",
        "Energy_distance": "Energy_contribution_score",
    }.get(metric)


def unit(metric):
    return "（nats）" if metric.startswith("Entropy_") else ""
