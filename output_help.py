from help_texts import method_guide
from localization import current_language


def guide_for(method):
    return method_guide(method, current_language())
