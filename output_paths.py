from pathlib import Path
import os
import re


def desktop_directory():
    if os.name == "nt":
        import ctypes

        buffer = ctypes.create_unicode_buffer(32768)
        if ctypes.windll.shell32.SHGetFolderPathW(None, 16, None, 0, buffer) == 0 and buffer.value:
            return Path(buffer.value)
    return Path.home() / "Desktop"


def experiment_folder_name(experiments):
    names = list(dict.fromkeys(e.name for e in experiments if e.enabled))
    name = "+".join(names) or "BioEntropy_结果"
    name = re.sub('[<>:"/\\\\|?*\\x00-\\x1f]', "_", name).strip(" .")[:100].rstrip(" .")
    if not name:
        return "BioEntropy_结果"
    if re.fullmatch("(CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])(?:\\..*)?", name, re.I):
        name = "实验_" + name
    return name


def numbered_path(base, number):
    base = Path(base)
    return base if number == 1 else base.with_name(f"{base.name}-{number}")


def available_output(base):
    number = 1
    while True:
        path = numbered_path(base, number)
        if not os.path.lexists(path):
            return path
        number += 1


def reserve_output(base):
    base = Path(base).expanduser().resolve()
    base.parent.mkdir(parents=True, exist_ok=True)
    number = 1
    while True:
        path = numbered_path(base, number)
        try:
            path.mkdir()
            return path
        except FileExistsError:
            number += 1
