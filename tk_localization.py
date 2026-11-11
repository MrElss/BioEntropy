from functools import wraps
from pathlib import Path
import tkinter as tk
from tkinter import ttk, messagebox as native_messagebox, filedialog as native_filedialog
from localization import translate, presentation, current_language


def identities(owner):
    values = []
    for experiment in getattr(owner, "experiments", []):
        values.extend([experiment.name, *experiment.groups])
        values.extend(
            Path(path).name
            for path in [*experiment.files, *experiment.normal_files, *experiment.range_files]
        )
    return values


class ChoiceVar(tk.StringVar):
    def __init__(self, owner, value="", **kwargs):
        self.owner = owner
        self.aliases = {}
        self.canonical_choices = []
        self.labels = {}
        self.display_choices = []
        super().__init__(master=owner.root, value=value, **kwargs)
        if not hasattr(owner, "_choice_vars"):
            owner._choice_vars = {}
        owner._choice_vars[str(self)] = self

    def get(self):
        value = super().get()
        return self.aliases.get(value, value)

    def set(self, value):
        label = self.labels.get(value, self.label(value))
        self.aliases[label] = value
        super().set(label)

    def label(self, value):
        return translate(str(value), self.owner.lang, identities(self.owner))

    def refresh(self, values):
        selected = self.get()
        canonical = (
            [self.aliases.get(value, value) for value in values]
            if values == self.display_choices
            else values
        )
        self.canonical_choices = canonical
        labels = [self.label(value) for value in canonical]
        duplicates = {label for label in labels if labels.count(label) > 1}
        if duplicates:
            from experiment_config import NONE, EXTERNAL

            labels = [
                label + (" [setting]" if value in {NONE, EXTERNAL} else " [group]")
                if label in duplicates
                else label
                for value, label in zip(canonical, labels)
            ]
        self.labels = dict(zip(canonical, labels))
        self.display_choices = labels
        self.aliases.update(zip(labels, canonical))
        self.set(selected)
        return labels


class DisplayVar(tk.StringVar):
    def __init__(self, owner, value="", **kwargs):
        self.owner = owner
        self.original = value
        super().__init__(master=owner.root, value=value, **kwargs)
        if not hasattr(owner, "_display_vars"):
            owner._display_vars = []
        owner._display_vars.append(self)

    def set(self, value):
        self.original = value
        self.refresh()

    def refresh(self):
        super().set(translate(self.original, self.owner.lang, identities(self.owner)))


def update_value(widget, key, value, language, names):
    cache = getattr(widget, "_translation_cache", {})
    previous = cache.get(key)
    original = previous[0] if previous and value == previous[1] else value
    rendered = translate(original, language, names)
    cache[key] = (original, rendered)
    widget._translation_cache = cache
    return rendered


def refresh_tree(owner):
    try:
        if not owner.root.winfo_exists():
            return
    except tk.TclError:
        return
    language, names = owner.lang, identities(owner)
    for variable in getattr(owner, "_display_vars", []):
        variable.refresh()

    def visit(widget):
        try:
            if isinstance(widget, (tk.Tk, tk.Toplevel)):
                value = widget.title()
                widget.title(update_value(widget, "window_title", value, language, names))
            if isinstance(widget, ttk.Combobox):
                variable = getattr(owner, "_choice_vars", {}).get(str(widget.cget("textvariable")))
                if variable:
                    widget.configure(values=variable.refresh(list(widget.cget("values"))))
                    choices = set(variable.canonical_choices)
                    if choices and choices <= {"Shannon熵", "Rényi二阶熵", "能量距离", "仅NMD"}:
                        widget.configure(width=30 if language == "en" else 18)
            if "text" in widget.keys() and not isinstance(widget, (tk.Entry, ttk.Entry)):
                value = widget.cget("text")
                widget.configure(text=update_value(widget, "text", value, language, names))
            if isinstance(widget, ttk.Treeview):
                for column in widget.cget("columns"):
                    value = widget.heading(column, "text")
                    widget.heading(
                        column,
                        text=update_value(widget, "heading_" + column, value, language, names),
                    )
                for row in widget.get_children():
                    values = list(widget.item(row, "values"))
                    widget.item(
                        row,
                        values=[
                            update_value(widget, f"item_{row}_{i}", value, language, names)
                            for i, value in enumerate(values)
                        ],
                    )
            if isinstance(widget, tk.Listbox) and widget is not getattr(owner, "orderlist", None):
                values = list(widget.get(0, "end"))
                rendered = [
                    update_value(widget, f"list_{i}", value, language, names)
                    for i, value in enumerate(values)
                ]
                if rendered != values:
                    widget.delete(0, "end")
                    for value in rendered:
                        widget.insert("end", value)
            if isinstance(widget, tk.Text) and widget is not getattr(owner, "logs", None):
                value = widget.get("1.0", "end-1c")
                rendered = update_value(widget, "content", value, language, names)
                if rendered != value:
                    state = widget.cget("state")
                    widget.configure(state="normal")
                    widget.delete("1.0", "end")
                    widget.insert("1.0", rendered)
                    widget.configure(state=state)
            for child in widget.winfo_children():
                visit(child)
        except tk.TclError:
            return

    visit(owner.root)
    if hasattr(owner, "methodbox"):
        owner.methodbox.configure(width=38 if language == "en" else 26)
    if hasattr(owner, "inferencebox"):
        owner.inferencebox.configure(width=38 if language == "en" else 31)
    if hasattr(owner, "files"):
        widths = (220, 165, 140, 105, 75, 350) if language == "en" else (230, 100, 120, 90, 70, 310)
        for column, width in zip(owner.files.cget("columns"), widths):
            owner.files.column(column, width=width)
    if hasattr(owner, "exptree"):
        widths = (55, 150, 160, 170, 190, 250) if language == "en" else (55, 150, 90, 170, 170, 230)
        for column, width in zip(owner.exptree.cget("columns"), widths):
            owner.exptree.column(column, width=width)
    if hasattr(owner, "lang_button"):
        owner.lang_button.configure(text="中文" if language == "en" else "EN")
    if hasattr(owner, "logs"):
        owner.logs.configure(state="normal")
        owner.logs.delete("1.0", "end")
        owner.logs.insert(
            "1.0",
            "".join(
                translate(str(value), language, names) + "\n"
                for value in getattr(owner, "_log_messages", [])
            ),
        )
        owner.logs.see("end")
        owner.logs.configure(state="disabled")


def localized_action(function):
    @wraps(function)
    def run(self, *args, **kwargs):
        if getattr(self, "_localizing", False):
            return function(self, *args, **kwargs)
        self._localizing = True
        refresh = (
            function.__name__ != "poll"
            or not getattr(self, "events", None)
            or not self.events.empty()
        )
        try:
            with presentation(getattr(self, "lang", "zh"), identities(self)):
                return function(self, *args, **kwargs)
        finally:
            self._localizing = False
            if refresh:
                refresh_tree(self)

    return run


def localized_callback(owner, function):
    @wraps(function)
    def run(*args, **kwargs):
        with presentation(owner.lang, identities(owner)):
            try:
                return function(*args, **kwargs)
            finally:
                refresh_tree(owner)

    return run


def localize_class(cls):
    for base in cls.__mro__:
        if base is object:
            continue
        for name, value in list(base.__dict__.items()):
            if name == "__init__" or name.startswith("__") or not callable(value):
                continue
            if name not in cls.__dict__ or base is cls:
                setattr(cls, name, localized_action(value))
    return cls


class Dialogs:
    def __getattr__(self, name):
        original = getattr(native_messagebox, name)

        def show(title=None, message=None, **options):
            if current_language() != "en":
                return original(title=title, message=message, **options)
            title, message = translate(title or ""), translate(message or "")
            parent = options.get("parent")
            window = tk.Toplevel(parent)
            window.title(title)
            window.transient(parent)
            window.resizable(False, False)
            frame = ttk.Frame(window, padding=22)
            frame.pack(fill="both", expand=True)
            ttk.Label(frame, text=message, wraplength=680, justify="left").pack(anchor="w")
            buttons = ttk.Frame(frame)
            buttons.pack(fill="x", pady=(20, 0))
            answer = {
                "value": False if name in {"askyesno", "askokcancel", "askretrycancel"} else "ok"
            }

            def finish(value):
                answer["value"] = value
                window.destroy()

            if name in {"askyesno", "askokcancel", "askretrycancel"}:
                positive, negative = {
                    "askyesno": ("Yes", "No"),
                    "askokcancel": ("OK", "Cancel"),
                    "askretrycancel": ("Retry", "Cancel"),
                }[name]
                ttk.Button(buttons, text=negative, command=lambda: finish(False)).pack(
                    side="right", padx=(8, 0)
                )
                ok = ttk.Button(buttons, text=positive, command=lambda: finish(True))
                ok.pack(side="right")
                window.bind("<Return>", lambda event: finish(True))
            else:
                ok = ttk.Button(buttons, text="OK", command=lambda: finish("ok"))
                ok.pack(side="right")
                window.bind("<Return>", lambda event: finish("ok"))
            window.bind("<Escape>", lambda event: finish(False if name.startswith("ask") else "ok"))
            window.protocol(
                "WM_DELETE_WINDOW", lambda: finish(False if name.startswith("ask") else "ok")
            )
            window.update_idletasks()
            if parent:
                x = parent.winfo_rootx() + max(
                    0, (parent.winfo_width() - window.winfo_reqwidth()) // 2
                )
                y = parent.winfo_rooty() + max(
                    0, (parent.winfo_height() - window.winfo_reqheight()) // 2
                )
                window.geometry(f"+{x}+{y}")
            window.grab_set()
            ok.focus_set()
            window.wait_window()
            return answer["value"]

        return show


class FileDialogs:
    def __getattr__(self, name):
        original = getattr(native_filedialog, name)

        def show(*args, **options):
            if "title" in options:
                options["title"] = translate(options["title"])
            if "filetypes" in options:
                options["filetypes"] = [
                    (translate(label), pattern) for label, pattern in options["filetypes"]
                ]
            return original(*args, **options)

        return show


messagebox = Dialogs()
filedialog = FileDialogs()
