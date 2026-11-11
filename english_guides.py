from help_texts import INPUT_GUIDE_EN, method_guide as localized_method_guide


def input_guide():
    return INPUT_GUIDE_EN


def method_guide(method):
    return localized_method_guide(method, "en")
