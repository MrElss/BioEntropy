from pathlib import Path
import json, os
import tkinter as tk
from tkinter import ttk, scrolledtext
from tk_localization import ChoiceVar, localized_callback, messagebox, filedialog
from localization import translate, sheet_title
from metric_selection import METHODS

NAVY = "#18334e"
BG = "#edf1f7"
TEAL = "#00a087"
BLUE = "#3d5488"
INK = "#20344b"


class InterfaceHelpers:
    def _style(self):
        s = ttk.Style(self.root)
        s.theme_use("clam")
        s.configure(".", font=("Microsoft YaHei UI", 10))
        s.configure("TFrame", background=BG)
        s.configure("TLabel", background=BG, foreground=INK)
        s.configure("Card.TFrame", background="white")
        s.configure("Card.TLabel", background="white", foreground=INK)
        s.configure(
            "Card.TLabelframe",
            background="white",
            bordercolor="#d5dde6",
            relief="solid",
            borderwidth=1,
        )
        s.configure(
            "Card.TLabelframe.Label",
            background=BG,
            foreground=NAVY,
            font=("Microsoft YaHei UI", 11, "bold"),
        )
        s.configure("TButton", padding=(11, 7), background="#e9eef7", foreground=INK, borderwidth=0)
        s.map("TButton", background=[("active", "#d8e3f2")])
        s.configure(
            "Accent.TButton",
            background=TEAL,
            foreground="white",
            font=("Microsoft YaHei UI", 11, "bold"),
            padding=(18, 10),
        )
        s.map("Accent.TButton", background=[("active", "#008b78"), ("disabled", "#94c6bd")])
        s.configure(
            "Treeview", rowheight=25, background="white", fieldbackground="white", foreground=INK
        )
        s.configure(
            "Treeview.Heading",
            background=BLUE,
            foreground="white",
            font=("Microsoft YaHei UI", 10, "bold"),
        )
        s.map("Treeview.Heading", background=[("active", NAVY)])
        s.configure("TProgressbar", background=TEAL, troughcolor="#dce3ee")

    def text(self, w, zh, en):
        self.trans.append((w, zh, en))
        w.configure(text=en if self.lang == "en" else zh)
        return w

    def button(self, parent, zh, en, command, **kw):
        return self.text(ttk.Button(parent, command=command, **kw), zh, en)

    def section(self, parent, zh, en):
        f = self.text(ttk.LabelFrame(parent, style="Card.TLabelframe", padding=10), zh, en)
        f.pack(fill="x", padx=16, pady=(0, 12))
        return f

    def toggle_language(self):
        self.lang = "en" if self.lang == "zh" else "zh"
        for w, zh, en in self.trans:
            w.configure(text=en if self.lang == "en" else zh)
        self.lang_button.configure(text="中文" if self.lang == "en" else "EN")
        titles = (
            ["File", "Type", "Experiment", "Week / note", "Records", "Note"]
            if self.lang == "en"
            else ["文件名", "类型", "实验号", "周数/说明", "记录数", "说明"]
        )
        for key, title in zip(self.files["columns"], titles):
            self.files.heading(key, text=title)

    def log(self, text):
        self.logs.configure(state="normal")
        if not hasattr(self, "_log_messages"):
            self._log_messages = []
        self._log_messages.append(str(text))
        self.logs.insert("end", translate(str(text)) + "\n")
        self.logs.see("end")
        self.logs.configure(state="disabled")

    def show_text(self, title, text):
        w = tk.Toplevel(self.root)
        w.title(title)
        w.geometry("850x550")
        box = scrolledtext.ScrolledText(
            w, wrap="word", font=("Microsoft YaHei UI", 11), padx=18, pady=14
        )
        box.pack(fill="both", expand=True)
        box.insert("1.0", text)
        box.configure(state="disabled")

    def g_dialog(self):
        if self.running:
            messagebox.showinfo("正在计算", "请等待当前计算完成。", parent=self.root)
            return
        w = tk.Toplevel(self.root)
        w.title("Hedges’ g值计算")
        w.geometry("780x260")
        w.configure(bg=BG)
        source = tk.StringVar()
        dest = tk.StringVar(value=self.output.get())
        frame = ttk.Frame(w, padding=20)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(1, weight=1)
        for i, (label, var) in enumerate([("G值数据目录", source), ("输出目录", dest)]):
            ttk.Label(frame, text=label).grid(row=i, column=0, pady=8, padx=(0, 10))
            ttk.Entry(frame, textvariable=var).grid(row=i, column=1, sticky="ew")

            def choose(v=var):
                p = filedialog.askdirectory(parent=w)
                if p:
                    v.set(p)

            ttk.Button(frame, text="选择…", command=localized_callback(self, choose)).grid(
                row=i, column=2, padx=(8, 0)
            )
        ttk.Label(frame, text="单指标、跨模型指标及NMD的G值计算与森林图。", wraplength=700).grid(
            row=2, column=0, columnspan=3, sticky="w", pady=12
        )

        def launch_g():
            if not source.get().strip():
                messagebox.showerror("G值输入", "请选择原G值数据目录。", parent=w)
                return
            a, b = (source.get(), dest.get())
            w.destroy()
            self.start_job("gvalues", [a], b)

        ttk.Button(
            frame,
            text="开始G值计算",
            style="Accent.TButton",
            command=localized_callback(self, launch_g),
        ).grid(row=3, column=1, sticky="e")

    def preview(self):
        path = Path(self.last) if self.last else None
        if path is None or not any(
            (path / name).is_file() for name in ("结果表.xlsx", "图表索引.json", "查看结果.html")
        ):
            p = filedialog.askdirectory(parent=self.root, title="选择已生成的结果目录")
            if not p:
                return
            path = Path(p)
        if not (path / "图表索引.json").is_file() and not (path / "结果表.xlsx").is_file():
            if (path / "查看结果.html").is_file():
                os.startfile(str(path / "查看结果.html"))
            else:
                messagebox.showinfo("结果预览", "此目录没有结果图。", parent=self.root)
            return
        if (path / "结果表.xlsx").is_file():
            from openpyxl import load_workbook

            book = load_workbook(path / "结果表.xlsx", read_only=True)
            available = "图表索引" in book.sheetnames or "Figure index" in book.sheetnames
            book.close()
            if not available:
                os.startfile(str(path / "查看结果.html"))
                return
        show_preview(self.root, path, METHODS[self.method.get()])


def show_preview(parent, path, default="Shannon"):
    from PIL import Image, ImageTk

    owner = getattr(parent, "_bioentropy_owner", None)
    if (path / "结果表.xlsx").is_file():
        import pandas as pd

        with pd.ExcelFile(path / "结果表.xlsx") as book:
            sheet = "图表索引" if "图表索引" in book.sheet_names else "Figure index"
            entries = pd.read_excel(
                book, sheet_name=sheet, dtype={"Experiment": str, "Method": str}
            ).to_dict("records")
    else:
        entries = json.loads((path / "图表索引.json").read_text("utf-8"))
    w = tk.Toplevel(parent)
    w.title("BioEntropy · 结果图")
    w.geometry("1120x820")
    available = list(dict.fromkeys(x["Method"] for x in entries if x["Method"] != "NMD")) or [
        "none"
    ]
    names = {"Shannon": "Shannon熵", "Renyi2": "Rényi二阶熵", "Energy": "能量距离", "none": "仅NMD"}
    lookup_methods = {names[key]: key for key in available}
    chosen = default if default in available else available[0]
    controls = ttk.Frame(w, padding=12)
    controls.pack(fill="x")
    exp = tk.StringVar(value=entries[0]["Experiment"])
    method = ChoiceVar(owner, value=names[chosen]) if owner else tk.StringVar(value=names[chosen])
    selected = ChoiceVar(owner) if owner else tk.StringVar()
    ttk.Label(controls, text="实验").pack(side="left")
    eb = ttk.Combobox(
        controls,
        textvariable=exp,
        values=list(dict.fromkeys(x["Experiment"] for x in entries)),
        state="readonly",
        width=24,
    )
    eb.pack(side="left", padx=(6, 18))
    ttk.Label(controls, text="本次方法").pack(side="left")
    mb = ttk.Combobox(
        controls, textvariable=method, values=list(lookup_methods), state="readonly", width=18
    )
    mb.pack(side="left", padx=(6, 18))
    row = ttk.Frame(w, padding=(12, 0, 12, 8))
    row.pack(fill="x")
    picker = ttk.Combobox(row, textvariable=selected, state="readonly", width=75)
    picker.pack(side="left")
    image_label = ttk.Label(w, anchor="center")
    image_label.pack(fill="both", expand=True, padx=10, pady=10)
    lookup = {}
    cache = {"image": None, "timer": None}

    def render(*args):
        item = lookup.get(selected.get())
        if not item:
            return
        if item["File"].endswith(".html"):
            image_label.configure(image="", text="三维交互图：点击右侧“打开当前图”可拖动旋转。")
            return
        image_label.configure(text="")
        with Image.open(path / item["File"]) as im:
            im = im.copy()
            im.thumbnail(
                (max(200, w.winfo_width() - 30), max(200, w.winfo_height() - 125)),
                Image.Resampling.LANCZOS,
            )
        cache["image"] = ImageTk.PhotoImage(im)
        image_label.configure(image=cache["image"])

    def refresh(*args):
        key = lookup_methods[method.get()]
        previous = selected.get()[:3]
        lookup.clear()
        for item in sorted(
            (x for x in entries if x["Experiment"] == exp.get() and x["Method"] in ("NMD", key)),
            key=lambda x: x["Number"],
        ):
            lookup[f"图{item['Number']:02d}  {item['Title']}"] = item
        picker.configure(values=list(lookup))
        selected.set(next((x for x in lookup if x.startswith(previous)), next(iter(lookup))))
        render()

    def resized(e):
        if e.widget != w:
            return
        if cache["timer"]:
            w.after_cancel(cache["timer"])
        cache["timer"] = w.after(120, render)

    if owner:
        refresh = localized_callback(owner, refresh)
        render = localized_callback(owner, render)
    eb.bind("<<ComboboxSelected>>", refresh)
    mb.bind("<<ComboboxSelected>>", refresh)
    picker.bind("<<ComboboxSelected>>", render)
    w.bind("<Configure>", resized)
    ttk.Button(
        row,
        text="打开当前图",
        command=lambda: os.startfile(str(path / lookup[selected.get()]["File"])),
    ).pack(side="right", padx=8)
    ttk.Button(row, text="打开结果文件夹", command=lambda: os.startfile(str(path))).pack(
        side="right"
    )
    w.update_idletasks()
    refresh()
