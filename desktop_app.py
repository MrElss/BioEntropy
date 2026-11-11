from pathlib import Path
import datetime, json, os, queue, threading, traceback, uuid, copy, math
import tkinter as tk
from tkinter import ttk, scrolledtext
from tk_localization import (
    ChoiceVar,
    DisplayVar,
    localize_class,
    localized_callback,
    messagebox,
    filedialog,
)
from localization import presentation
from ui_components import InterfaceHelpers as Appearance, METHODS, NAVY, BG, TEAL
from analysis_engine import VERSION, BASE, write_json
from experiment_config import (
    discover,
    candidates,
    Experiment,
    EXTERNAL,
    NONE,
    save_project,
    load_project,
)
from analysis_workflow import preflight, run_project
from output_paths import desktop_directory, experiment_folder_name, available_output, reserve_output

SETTINGS = Path.home() / ".bioentropy_settings.json"
INFERENCE_OPTIONS = {
    "标准检验（9,999次置换）": 9999,
    "快速预览（999次置换）": 999,
    "关闭显著性检验": 0,
}


@localize_class
class Desktop(Appearance):
    def __init__(self, execute, smoke=None, language=None):
        self.execute = execute
        self.root = tk.Tk()
        self.root._bioentropy_owner = self
        self.root.title("BioEntropy 桌面版 · " + VERSION)
        self.root.iconbitmap(default=str(BASE / "BioEntropy.ico"))
        self.root.geometry(
            f"{min(1260, int(self.root.winfo_screenwidth() * 0.9))}x{min(960, int(self.root.winfo_screenheight() * 0.9))}"
        )
        self.root.minsize(950, 660)
        self.root.configure(bg=BG)
        self.events = queue.Queue()
        self.running = False
        self.last = None
        self.inputs = []
        self.experiments = []
        self.active = None
        self.lang = language or "zh"
        self._requested_language = language
        self.trans = []
        self.loading = False
        self.automatic_output = True
        self.setting_output = False
        self.output_base = desktop_directory() / "BioEntropy_结果"
        self.output_name = "BioEntropy_结果"
        self.method = ChoiceVar(self, value=next(iter(METHODS)))
        self.output = tk.StringVar(value=str(available_output(self.output_base)))
        self.inference = ChoiceVar(self, value=next(iter(INFERENCE_OPTIONS)))
        self.output.trace_add("write", self.output_edited)
        self.status = DisplayVar(self, value="请选择数据，程序自动识别实验。")
        self.path = DisplayVar(self)
        self.reference_note = DisplayVar(self)
        self._style()
        self._build()
        self.root.protocol("WM_DELETE_WINDOW", self.close)
        if smoke:
            self.root.withdraw()
            self.root.update_idletasks()
            self.root.update()
            write_json(
                smoke,
                dict(
                    version=VERSION,
                    language=self.lang,
                    gui_created=True,
                    analysis_type_selector=False,
                    automatic_detection=True,
                    role_selectors=list(self.role_vars),
                    group_order_drag=True,
                    group_order_buttons=True,
                    project_save_load=True,
                    normal_range_template=True,
                    gvalues_entry=True,
                    entropy_choices=list(METHODS),
                    inference_choices=list(INFERENCE_OPTIONS),
                    shannon_information_scores=True,
                    weekly_joint_panels=True,
                    figure_policy="animal and clinical 01-05: 常规; animal 07-08: 常规; clinical 10-11: 分周原值",
                    g_engine_preserved=True,
                ),
            )
            self.root.destroy()
            return
        self.restore()
        self.poll()
        self.root.mainloop()

    def _build(self):
        header = tk.Frame(self.root, bg=NAVY)
        header.pack(fill="x")
        logo = tk.Canvas(header, width=120, height=65, bg=NAVY, highlightthickness=0)
        logo.pack(side="left", padx=18, pady=8)
        for i in range(42):
            angle = i * 2.39996
            r = 5 + math.sqrt(i) * 6
            x = 60 + math.cos(angle) * r
            y = 32 + math.sin(angle) * r * 0.65
            logo.create_oval(x - 2, y - 2, x + 2, y + 2, fill="#63cbbb", outline="")
        title = tk.Frame(header, bg=NAVY)
        title.pack(side="left", pady=12)
        tk.Label(
            title, text="BioEntropy", font=("Microsoft YaHei UI", 24, "bold"), bg=NAVY, fg="white"
        ).pack(anchor="w")
        self.text(
            tk.Label(title, bg=NAVY, fg="#c3d2e3"),
            "多指标熵 · 正常参照马氏距离分析",
            "Multivariate entropy · Normal-reference Mahalanobis distance",
        ).pack(anchor="w")
        tk.Label(header, text=VERSION, bg=NAVY, fg="white", font=("Microsoft YaHei UI", 12)).pack(
            side="right", padx=22
        )
        self.lang_button = ttk.Button(header, text="EN", command=self.toggle_language)
        self.lang_button.pack(side="right")
        tk.Frame(self.root, bg=TEAL, height=3).pack(fill="x")
        outer = ttk.Frame(self.root)
        outer.pack(fill="both", expand=True)
        self.canvas = tk.Canvas(outer, bg=BG, highlightthickness=0)
        self.canvas.pack(side="left", fill="both", expand=True)
        bar = ttk.Scrollbar(outer, command=self.canvas.yview)
        bar.pack(side="right", fill="y")
        self.canvas.configure(yscrollcommand=bar.set)
        body = ttk.Frame(self.canvas)
        item = self.canvas.create_window((0, 0), window=body, anchor="nw")
        body.bind(
            "<Configure>", lambda ev: self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        )
        self.canvas.bind("<Configure>", lambda ev: self.canvas.itemconfigure(item, width=ev.width))

        def wheel(ev):
            if isinstance(ev.widget, (tk.Text, ttk.Treeview, tk.Listbox, ttk.Combobox)):
                return
            self.canvas.yview_scroll(-int(ev.delta / 120), "units")

        self.root.bind_all("<MouseWheel>", wheel)
        ttk.Frame(body, height=12).pack()
        files = self.section(body, "1）文件、模板与输出目录", "1) Files, templates and output")
        self.controls = []
        actions = [
            ("添加Excel文件", "Add Excel files", self.select_files),
            ("添加数据文件夹", "Add data folder", self.select_folder),
            ("清空", "Clear", self.clear_files),
            ("检查输入", "Check inputs", self.check_inputs),
            ("保存分析设置", "Save project", self.save_project_ui),
            ("载入分析设置", "Load project", self.load_project_ui),
        ]
        row = ttk.Frame(files, style="Card.TFrame")
        row.pack(fill="x")
        for zh, en, fn in actions:
            b = self.button(row, zh, en, fn)
            b.pack(side="left", padx=(0, 5))
            self.controls.append(b)
        row = ttk.Frame(files, style="Card.TFrame")
        row.pack(fill="x", pady=7)
        for zh, en, fn in [
            ("填写说明", "Input guide", self.show_guide),
            ("导出模板…", "Export templates…", self.templates_dialog),
            ("计算与核查规则", "Methods and review rules", self.show_methods),
            ("Hedges’ g值计算", "Hedges’ g calculation", self.g_dialog),
        ]:
            b = self.button(row, zh, en, fn)
            b.pack(side="left", padx=(0, 5))
            self.controls.append(b)
        form = ttk.Frame(files, style="Card.TFrame")
        form.pack(fill="x")
        form.columnconfigure(1, weight=1)
        self.text(
            ttk.Label(form, style="Card.TLabel"), "熵/距离方法：", "Entropy / distance:"
        ).grid(row=0, column=0, pady=6, sticky="w")
        self.methodbox = ttk.Combobox(
            form, textvariable=self.method, values=list(METHODS), state="readonly", width=26
        )
        self.methodbox.grid(row=0, column=1, sticky="w")
        self.text(
            ttk.Label(form, style="Card.TLabel", foreground="#66758a"),
            "自动识别实验类型，无需手动选择。",
            "Experiment type is detected automatically.",
        ).grid(row=0, column=2, columnspan=2, sticky="e")
        self.text(ttk.Label(form, style="Card.TLabel"), "结果目录：", "Output folder:").grid(
            row=1, column=0, sticky="w"
        )
        self.output_entry = ttk.Entry(form, textvariable=self.output)
        self.output_entry.grid(row=1, column=1, sticky="ew")
        self.output_button = self.button(form, "选择…", "Choose…", self.choose_output)
        self.output_button.grid(row=1, column=2, padx=5)
        self.controls.append(self.output_button)
        self.button(form, "打开结果", "Open results", self.open_output).grid(row=1, column=3)
        ttk.Label(form, textvariable=self.path, style="Card.TLabel", wraplength=1100).grid(
            row=2, column=0, columnspan=4, sticky="w", pady=(6, 0)
        )
        self.text(
            ttk.Label(form, style="Card.TLabel", foreground="#66758a"),
            "默认：桌面／实验名；重名自动加 -2、-3。开始计算时创建目录。",
            "Default: Desktop / experiment name; existing names get -2, -3. Created when analysis starts.",
        ).grid(row=3, column=0, columnspan=4, sticky="w", pady=(3, 0))
        self.text(ttk.Label(form, style="Card.TLabel"), "显著性：", "Inference:").grid(
            row=4, column=0, sticky="w", pady=6
        )
        self.inferencebox = ttk.Combobox(
            form,
            textvariable=self.inference,
            values=list(INFERENCE_OPTIONS),
            state="readonly",
            width=31,
        )
        self.inferencebox.grid(row=4, column=1, sticky="w", pady=6)
        self.text(
            ttk.Label(form, style="Card.TLabel", foreground="#66758a"),
            "探索性组间检验；标准计算耗时较长。",
            "Exploratory comparison; standard testing takes longer.",
        ).grid(row=4, column=2, columnspan=2, sticky="e")
        selected = self.section(body, "2）已识别文件", "2) Recognized files")
        self.files = ttk.Treeview(
            selected,
            columns=("file", "type", "experiment", "week", "records", "note"),
            show="headings",
            height=4,
        )
        for key, label, width in [
            ("file", "文件名", 230),
            ("type", "类型", 100),
            ("experiment", "实验", 120),
            ("week", "周数/说明", 90),
            ("records", "记录数", 70),
            ("note", "说明", 310),
        ]:
            self.files.heading(key, text=label)
            self.files.column(key, width=width, minwidth=60)
        self.files.pack(side="left", fill="both", expand=True)
        filebar = ttk.Scrollbar(selected, command=self.files.yview)
        filebar.pack(side="right", fill="y")
        self.files.configure(yscrollcommand=filebar.set)
        groups = self.section(
            body, "3）实验分组、正常参照与显示顺序", "3) Roles, reference and display order"
        )
        treeframe = ttk.Frame(groups, style="Card.TFrame")
        treeframe.pack(fill="x")
        self.exptree = ttk.Treeview(
            treeframe,
            columns=("run", "experiment", "type", "treatment", "comparator", "reference"),
            show="headings",
            height=4,
        )
        for key, label, width in [
            ("run", "计算", 55),
            ("experiment", "实验", 150),
            ("type", "类型", 90),
            ("treatment", "治疗组", 170),
            ("comparator", "比较组/模型组", 170),
            ("reference", "正常参照", 230),
        ]:
            self.exptree.heading(key, text=label)
            self.exptree.column(key, width=width, minwidth=50)
        self.exptree.pack(side="left", fill="x", expand=True)
        expbar = ttk.Scrollbar(treeframe, command=self.exptree.yview)
        expbar.pack(side="right", fill="y")
        self.exptree.configure(yscrollcommand=expbar.set)
        self.exptree.bind("<<TreeviewSelect>>", self.select_experiment)
        panel = ttk.Frame(groups, style="Card.TFrame")
        panel.pack(fill="x", pady=(10, 0))
        left = ttk.Frame(panel, style="Card.TFrame")
        left.pack(side="left", fill="both", expand=True)
        self.enabled = tk.BooleanVar(value=True)
        self.enablebox = ttk.Checkbutton(
            left, text="计算所选实验", variable=self.enabled, command=self.commit
        )
        self.enablebox.grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 5))
        self.role_vars = {}
        self.role_boxes = {}
        for i, (key, zh, en) in enumerate(
            [
                ("treatment", "治疗组：", "Treatment:"),
                ("comparator", "比较组/模型组：", "Comparator / model:"),
                ("reference", "正常参照：", "Normal reference:"),
            ],
            1,
        ):
            self.text(ttk.Label(left, style="Card.TLabel"), zh, en).grid(
                row=i, column=0, sticky="w", pady=3
            )
            var = ChoiceVar(self, value=NONE)
            box = ttk.Combobox(left, textvariable=var, state="disabled", width=34)
            box.grid(row=i, column=1, sticky="w")
            box.bind("<<ComboboxSelected>>", self.commit)
            self.role_vars[key] = var
            self.role_boxes[key] = box
        buttons = ttk.Frame(left, style="Card.TFrame")
        buttons.grid(row=4, column=0, columnspan=2, sticky="w", pady=6)
        for zh, en, fn in [
            ("Normal样本…", "Normal samples…", lambda: self.attach("normal_files")),
            ("正常范围…", "Normal ranges…", lambda: self.attach("range_files")),
            ("参照文件清单", "Reference files", self.reference_files),
        ]:
            b = self.button(buttons, zh, en, fn)
            b.pack(side="left", padx=(0, 5))
            self.controls.append(b)
        ttk.Label(
            left,
            textvariable=self.reference_note,
            style="Card.TLabel",
            wraplength=610,
            foreground="#66758a",
        ).grid(row=5, column=0, columnspan=2, sticky="w")
        right = ttk.Frame(panel, style="Card.TFrame")
        right.pack(side="right", fill="both", padx=(15, 0))
        self.text(
            ttk.Label(right, style="Card.TLabel"),
            "显示顺序（拖动或使用上下按钮）",
            "Display order (drag or move)",
        ).pack(anchor="w")
        orderrow = ttk.Frame(right, style="Card.TFrame")
        orderrow.pack(fill="both", expand=True)
        self.orderlist = tk.Listbox(
            orderrow,
            height=6,
            width=29,
            exportselection=False,
            font=("Microsoft YaHei UI", 10),
            activestyle="none",
        )
        self.orderlist.pack(side="left", fill="both", expand=True)
        self.drag_index = None
        self.orderlist.bind("<Button-1>", self.drag_start)
        self.orderlist.bind("<B1-Motion>", self.drag_move)
        self.orderlist.bind("<ButtonRelease-1>", self.commit)
        bs = ttk.Frame(orderrow, style="Card.TFrame")
        bs.pack(side="right", padx=6)
        for zh, en, fn in [
            ("上移", "Up", lambda: self.move(-1)),
            ("下移", "Down", lambda: self.move(1)),
            ("恢复顺序", "Reset order", self.reset_order),
        ]:
            b = self.button(bs, zh, en, fn)
            b.pack(fill="x", pady=2)
            self.controls.append(b)
        action = self.text(
            ttk.LabelFrame(self.root, style="Card.TLabelframe", padding=10),
            "开始分析",
            "Run analysis",
        )
        action.pack(side="bottom", fill="x", padx=16, pady=(0, 10))
        self.run_button = self.button(
            action,
            "开始分析并生成结果",
            "Run analysis and export",
            self.start,
            style="Accent.TButton",
        )
        self.run_button.pack(side="left")
        self.button(action, "查看结果图", "View figures", self.preview).pack(side="left", padx=10)
        self.progress = ttk.Progressbar(action, mode="indeterminate")
        self.progress.pack(side="left", fill="x", expand=True, padx=10)
        log = self.section(body, "运行日志", "Run log")
        self.logs = scrolledtext.ScrolledText(
            log,
            height=5,
            wrap="word",
            font=("Consolas", 10),
            bg="#0f1722",
            fg="#d6e2f0",
            state="disabled",
        )
        self.logs.pack(fill="x")
        ttk.Label(body, textvariable=self.status, wraplength=1100).pack(
            anchor="w", padx=16, pady=(0, 14)
        )
        self.log("程序已启动。可同时添加临床分周表和动物分组表；选择实验后可修改组别和显示顺序。")

    def current(self):
        return next((e for e in self.experiments if e.name == self.active), None)

    def set_output(self, path):
        self.setting_output = True
        try:
            self.output.set(str(path))
        finally:
            self.setting_output = False

    def output_edited(self, *args):
        if self.setting_output:
            return
        self.automatic_output = False
        self.output_base = Path(self.output.get().strip()).expanduser()

    def refresh_output(self, reset=False):
        if reset:
            self.automatic_output = True
        if self.automatic_output:
            self.output_name = experiment_folder_name(self.experiments)
            self.output_base = desktop_directory() / self.output_name
            self.set_output(available_output(self.output_base))

    def choose_output(self):
        if self.running:
            return
        p = filedialog.askdirectory(
            parent=self.root,
            title="选择保存位置（将在其下新建实验结果文件夹）",
            initialdir=str(desktop_directory()),
        )
        if p:
            self.automatic_output = False
            self.output_base = Path(p) / experiment_folder_name(self.experiments)
            self.set_output(available_output(self.output_base))

    def open_output(self):
        p = self.last or self.output.get()
        if Path(p).is_dir():
            os.startfile(str(p))
        else:
            messagebox.showinfo(
                "结果目录", "开始计算时会自动创建结果目录：\n" + self.output.get(), parent=self.root
            )

    def commit(self, *args):
        if self.loading or self.running:
            return
        e = self.current()
        if not e:
            return
        for k, v in self.role_vars.items():
            setattr(e, k, v.get())
        e.enabled = self.enabled.get()
        e.order = list(self.orderlist.get(0, "end"))
        self.update_tree(e)
        if self.automatic_output and experiment_folder_name(self.experiments) != self.output_name:
            self.refresh_output()
        self.reference_note.set(
            (
                "Normal样本："
                f"{len(e.normal_files)}"
                "份；正常范围："
                f"{len(e.range_files)}"
                "份。正常范围用于核查；NMD采用Normal样本。\n治疗/比较组用于组间比较；正常参照用于NMD；显示顺序只改变图表排列。"
            )
        )

    def update_tree(self, e):
        vals = (
            "✓" if e.enabled else "—",
            e.name,
            "临床/多时间点" if e.kind == "clinical" else "动物/单时间点",
            e.treatment,
            e.comparator,
            e.reference,
        )
        if self.exptree.exists(e.name):
            self.exptree.item(e.name, values=vals)
        else:
            self.exptree.insert("", "end", iid=e.name, values=vals)

    def select_experiment(self, *args):
        if self.running:
            return
        chosen = self.exptree.selection()
        if not chosen:
            return
        self.commit()
        self.active = chosen[0]
        e = self.current()
        self.loading = True
        self.enabled.set(e.enabled)
        for k, v in self.role_vars.items():
            values = [NONE, *e.groups] if k != "reference" else [NONE, EXTERNAL, *e.groups]
            self.role_boxes[k].configure(values=values, state="readonly")
            v.set(getattr(e, k))
        self.orderlist.delete(0, "end")
        for g in e.order:
            self.orderlist.insert("end", g)
        self.loading = False
        self.commit()

    def render_catalog(self):
        self.files.delete(*self.files.get_children())
        self.exptree.delete(*self.exptree.get_children())
        self.active = None
        for e in self.experiments:
            for r in sorted(
                e.rows,
                key=lambda r: (
                    r["kind"] in ("normal", "range"),
                    r["week"] if isinstance(r["week"], int) else 9999,
                ),
            ):
                self.files.insert(
                    "",
                    "end",
                    values=(
                        Path(r["file"]).name,
                        {
                            "clinical": "临床/多时间点",
                            "animal": "动物/单时间点",
                            "normal": "Normal样本",
                            "range": "正常范围",
                        }[r["kind"]],
                        r["experiment"],
                        r["week"],
                        r["count"],
                        r["note"],
                    ),
                )
            self.update_tree(e)
        if self.experiments:
            self.exptree.selection_set(self.experiments[0].name)
            self.select_experiment()
        self.path.set(
            f"已识别 {len(self.experiments)} 个实验；{sum(len(e.files) for e in self.experiments)} 份实验数据。请选择实验核对组别。"
        )

    def add_inputs(self, paths):
        self.commit()
        previous = {e.name: e for e in self.experiments}
        selected = list(dict.fromkeys([*self.inputs, *[str(p) for p in candidates(paths)]]))
        try:
            found = discover(selected)
            for e in found:
                old = previous.get(e.name)
                if old and set(old.groups) == set(e.groups):
                    for key in ["treatment", "comparator", "reference", "order", "enabled"]:
                        setattr(e, key, copy.deepcopy(getattr(old, key)))
                    for key in ["normal_files", "range_files"]:
                        setattr(e, key, list(dict.fromkeys([*getattr(e, key), *getattr(old, key)])))
            self.inputs = selected
            self.experiments = found
            self.last = None
            self.render_catalog()
            self.refresh_output()
            self.status.set("识别完成。核对组别后即可开始。")
        except Exception as exc:
            messagebox.showerror("文件识别", str(exc), parent=self.root)

    def select_files(self):
        if self.running:
            return
        paths = filedialog.askopenfilenames(
            parent=self.root,
            title="选择实验Excel文件（可多选）",
            filetypes=[("Excel", "*.xlsx *.xls")],
        )
        if paths:
            self.add_inputs(paths)

    def select_folder(self):
        if self.running:
            return
        p = filedialog.askdirectory(
            parent=self.root, title="选择数据目录（识别本目录及下一层文件夹）"
        )
        if p:
            self.add_inputs([p])

    def clear_files(self):
        if self.running:
            return
        self.inputs = []
        self.experiments = []
        self.active = None
        self.render_catalog()
        self.orderlist.delete(0, "end")
        self.reference_note.set("")
        for box in self.role_boxes.values():
            box.configure(state="disabled")
        self.last = None
        self.refresh_output(reset=True)

    def attach(self, key):
        if self.running:
            return
        e = self.current()
        if not e:
            messagebox.showinfo("选择实验", "请先载入并选中实验。", parent=self.root)
            return
        paths = filedialog.askopenfilenames(
            parent=self.root,
            title=e.name + "：选择" + ("Normal样本" if key == "normal_files" else "正常范围"),
            filetypes=[("Excel", "*.xlsx *.xls")],
        )
        if paths:
            setattr(e, key, list(dict.fromkeys([*getattr(e, key), *paths])))
            if key == "normal_files":
                self.role_vars["reference"].set(EXTERNAL)
            self.commit()
            self.log(e.name + "：已指定参照文件，请运行前检查输入。")

    def reference_files(self):
        if self.running:
            return
        e = self.current()
        if not e:
            return
        w = tk.Toplevel(self.root)
        w.title(e.name + "：参照文件")
        w.geometry("840x330")
        box = tk.Listbox(w, width=110, height=10)
        box.pack(fill="both", expand=True, padx=12, pady=12)
        items = [(k, p) for k in ("normal_files", "range_files") for p in getattr(e, k)]
        for k, p in items:
            box.insert("end", ("样本 | " if k == "normal_files" else "范围 | ") + p)

        def remove():
            if not box.curselection():
                return
            i = box.curselection()[0]
            k, p = items.pop(i)
            getattr(e, k).remove(p)
            box.delete(i)
            if not e.normal_files and e.reference == EXTERNAL:
                self.role_vars["reference"].set(NONE)
            self.commit()

        ttk.Button(
            w, text="移除所选参照（不删除原文件）", command=localized_callback(self, remove)
        ).pack(pady=8)

    def drag_start(self, ev):
        self.drag_index = self.orderlist.nearest(ev.y)

    def drag_move(self, ev):
        if self.running or self.drag_index is None:
            return
        j = self.orderlist.nearest(ev.y)
        if j != self.drag_index:
            text = self.orderlist.get(self.drag_index)
            self.orderlist.delete(self.drag_index)
            self.orderlist.insert(j, text)
            self.orderlist.selection_clear(0, "end")
            self.orderlist.selection_set(j)
            self.drag_index = j
            self.commit()

    def move(self, delta):
        if self.running or not self.orderlist.curselection():
            return
        i = self.orderlist.curselection()[0]
        j = i + delta
        if 0 <= j < self.orderlist.size():
            text = self.orderlist.get(i)
            self.orderlist.delete(i)
            self.orderlist.insert(j, text)
            self.orderlist.selection_clear(0, "end")
            self.orderlist.selection_set(j)
            self.commit()

    def reset_order(self):
        e = self.current()
        if self.running or not e:
            return
        self.orderlist.delete(0, "end")
        for g in e.groups:
            self.orderlist.insert("end", g)
        self.commit()

    def check_inputs(self):
        if self.running:
            return
        self.commit()
        try:
            notes = preflight(self.experiments)
            self.show_text("输入核查完成", "\n".join(notes) + "\n\n核查通过。原始文件保持原样。")
            self.log("输入核查通过。")
        except Exception as exc:
            messagebox.showerror("输入核查", str(exc), parent=self.root)

    def save_project_ui(self):
        self.commit()
        p = filedialog.asksaveasfilename(
            parent=self.root,
            title="保存分析设置",
            initialfile="BioEntropy_分析设置.json",
            defaultextension=".json",
            filetypes=[("项目设置", "*.json")],
        )
        if p:
            save_project(
                p,
                self.experiments,
                METHODS[self.method.get()],
                str(self.output_base),
                INFERENCE_OPTIONS[self.inference.get()],
            )
            data = json.loads(Path(p).read_text("utf-8"))
            data["automatic_output"] = self.automatic_output
            data["language"] = self.lang
            write_json(p, data)
            self.log("已保存设置：" + p)

    def load_project_ui(self):
        p = filedialog.askopenfilename(
            parent=self.root, title="载入分析设置", filetypes=[("项目设置", "*.json")]
        )
        if not p:
            return
        try:
            self.experiments, data = load_project(p)
            self.inputs = [p for e in self.experiments for p in e.files]
            self.inference.set(
                next(
                    (k for k, v in INFERENCE_OPTIONS.items() if v == data.get("permutations")),
                    next(iter(INFERENCE_OPTIONS)),
                )
            )
            self.last = None
            self.automatic_output = data.get("automatic_output", True)
            if not self.automatic_output and data.get("output"):
                self.output_base = Path(data["output"])
                self.set_output(available_output(self.output_base))
            self.method.set(
                next(
                    (k for k, v in METHODS.items() if v == data.get("entropy_display")),
                    next(iter(METHODS)),
                )
            )
            self.render_catalog()
            self.refresh_output()
        except Exception as exc:
            messagebox.showerror("载入设置", str(exc), parent=self.root)

    def show_guide(self):
        if self.lang == "en":
            from english_guides import input_guide

            self.show_text("Input naming and format", input_guide())
            return
        self.show_text(
            "命名与填写说明",
            (
                "可同时选择临床与动物Excel。程序按内容和命名自动识别分析类型。\n\n单时间点：实验名.xlsx，前两列Group、Sample "
                "ID，之后为数值指标；每行一个样本。\n"
                "多时间点：实验名-0W.xlsx、实验名-6W.xlsx等，另需真实USUBJID列。该列用于身份核对和关联访视。105、106及"
                "105and106使用固定指标与尺度；其他实验使用自身0W固定尺度并记录配置。\n\n"
                "独立正常样本：实验名-Normal.xlsx；同目录只有一个实验时可叫Normal.xlsx。多个正常样本文件可在实验设置中手动添"
                "加并合并。正常范围：实验名-NormalRange.xlsx，Feature、Normal_min、Normal_max三列，名称"
                "及单位应完全一致。范围用于核查。\n\n"
                "选中实验后，可更改治疗组、比较组、正常参照；拖动或上下移动组名调整全部图表顺序。取消“计算所选实验”可跳过该实验。\n"
                "治疗/比较设置不改变原始估计；正常参照变更会改变NMD及选择的能量距离。无可用正常参照时NMD和能量距离写NA并说明原因。\n\n"
                "缺失及无效记录按现有纳入规则处理，原始记录保持原样。所有结果均写入新的结果目录。"
            ),
        )

    def show_methods(self):
        if self.lang == "en":
            from english_guides import method_guide

            self.show_text(
                "Selected method and validation rules", method_guide(METHODS[self.method.get()])
            )
            return
        from analysis_outputs import guide_for

        self.show_text("本次方法与核查规则", guide_for(METHODS[self.method.get()]))

    def templates_dialog(self):
        w = tk.Toplevel(self.root)
        w.title("导出模板与填写说明")
        w.geometry("500x350")
        kind = ChoiceVar(self, value="动物数据")
        exp = ChoiceVar(self, value="105")
        ttk.Label(w, text="选择模板类型", padding=12).pack(anchor="w")
        ttk.Combobox(
            w,
            textvariable=kind,
            values=["动物数据", "临床分周数据", "Normal样本", "Normal范围", "输入填写说明"],
            state="readonly",
            width=30,
        ).pack(padx=15, anchor="w")
        ttk.Label(w, text="临床/参照模板的指标样式", padding=12).pack(anchor="w")
        ttk.Combobox(
            w,
            textvariable=exp,
            values=["105", "106", "105and106", "通用"],
            state="readonly",
            width=30,
        ).pack(padx=15, anchor="w")

        def export():
            from table_exports import workbook
            import pandas as pd

            k = kind.get()
            e = exp.get()
            normal = k == "Normal样本"
            headers = ["Group", "Sample ID", "指标1（单位）", "指标2（单位）", "指标3（单位）"]
            if e != "通用" and k in ("临床分周数据", "Normal样本"):
                headers = json.loads((BASE / "template_headers.json").read_text("utf-8"))[e][
                    "normal" if normal else "study"
                ]
            elif k == "临床分周数据":
                headers = ["Group", "Sample ID", "USUBJID", *headers[2:]]
            filename = (
                "Normal.xlsx"
                if normal
                else f"{e}-NormalRange.xlsx"
                if k == "Normal范围"
                else f"{e}-0W.xlsx"
                if k == "临床分周数据"
                else "动物实验.xlsx"
            )
            if k == "输入填写说明":
                p = filedialog.asksaveasfilename(
                    parent=w, defaultextension=".txt", initialfile="BioEntropy_填写说明.txt"
                )
                if p:
                    from english_guides import input_guide

                    Path(p).write_text(
                        input_guide()
                        if self.lang == "en"
                        else (BASE / "使用说明.md").read_text("utf-8"),
                        "utf-8",
                    )
            else:
                p = filedialog.asksaveasfilename(
                    parent=w,
                    defaultextension=".xlsx",
                    initialfile=filename,
                    filetypes=[("Excel", "*.xlsx")],
                )
                if p:
                    with presentation("zh"):
                        workbook(
                            Path(p),
                            {
                                "数据": pd.DataFrame(
                                    columns=["Feature", "Normal_min", "Normal_max"]
                                    if k == "Normal范围"
                                    else headers
                                )
                            },
                            translated=False,
                        )
            if p:
                self.log("已导出：" + p)
                w.destroy()

        ttk.Button(
            w, text="导出", style="Accent.TButton", command=localized_callback(self, export)
        ).pack(pady=25)

    def start(self):
        if self.running:
            return
        self.commit()
        if not self.experiments:
            messagebox.showinfo("输入", "请先添加数据。", parent=self.root)
            return
        self.start_job(
            "auto",
            copy.deepcopy(self.experiments),
            self.output.get(),
            method=METHODS[self.method.get()],
        )

    def start_job(self, mode, inputs, destination, scope="auto", method="Shannon"):
        if self.running:
            return
        if not str(destination).strip():
            messagebox.showerror("结果目录", "请选择结果目录。", parent=self.root)
            return
        try:
            if mode == "auto":
                out = reserve_output(self.output_base)
                self.set_output(out)
            else:
                out = Path(destination) / (
                    "BioEntropy_Results_"
                    + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                    + "_"
                    + uuid.uuid4().hex[:5]
                )
        except OSError as exc:
            messagebox.showerror("结果目录", "无法创建结果目录：\n" + str(exc), parent=self.root)
            return
        permutations = INFERENCE_OPTIONS[self.inference.get()]
        language = self.lang
        self.lang_button.configure(state="disabled")
        self.inferencebox.configure(state="disabled")
        self.persist()
        self.running = True
        self.run_button.configure(state="disabled")
        self.methodbox.configure(state="disabled")
        self.enablebox.configure(state="disabled")
        self.orderlist.configure(state="disabled")
        self.output_entry.configure(state="disabled")
        for widget in [*self.controls, *self.role_boxes.values()]:
            widget.configure(state="disabled")
        self.progress.start(15)
        self.status.set("正在分析…")
        self.log("结果将保存到：" + str(out))

        def worker():
            try:
                log = lambda s: self.events.put(("log", s))
                result = (
                    run_project(inputs, out, method, log, permutations, language=language)
                    if mode == "auto"
                    else self.execute(
                        mode, inputs, out, scope, method, log, permutations, language=language
                    )
                )
                self.events.put(("done", str(result)))
            except Exception:
                detail = traceback.format_exc()
                out.mkdir(parents=True, exist_ok=True)
                (out / "错误信息.txt").write_text(detail, "utf-8")
                self.events.put(("error", detail))

        threading.Thread(target=worker, daemon=True).start()

    def poll(self):
        try:
            while True:
                kind, value = self.events.get_nowait()
                if kind == "log":
                    self.log(value)
                else:
                    self.lang_button.configure(state="normal")
                    self.inferencebox.configure(state="readonly")
                    self.running = False
                    self.progress.stop()
                    self.run_button.configure(state="normal")
                    self.methodbox.configure(state="readonly")
                    self.enablebox.configure(state="normal")
                    self.orderlist.configure(state="normal")
                    self.output_entry.configure(state="normal")
                    for widget in self.controls:
                        widget.configure(state="normal")
                    for box in self.role_boxes.values():
                        box.configure(state="readonly" if self.current() else "disabled")
                    if kind == "done":
                        self.last = value
                        self.status.set("计算完成：" + value)
                        self.log("完成。可点击“查看结果图”或“打开结果”。")
                        messagebox.showinfo("分析完成", "结果已保存：\n" + value, parent=self.root)
                    else:
                        self.status.set("分析未完成，请查看日志。")
                        self.log(value)
                        messagebox.showerror("分析未完成", value.splitlines()[-1], parent=self.root)
                    self.persist()
        except queue.Empty:
            pass
        self.root.after(150, self.poll)

    def persist(self):
        if not self.running:
            self.commit()
        write_json(
            SETTINGS,
            dict(
                geometry=self.root.geometry(),
                output=self.output.get(),
                output_base=str(self.output_base),
                automatic_output=self.automatic_output,
                method=self.method.get(),
                inference=self.inference.get(),
                language=self.lang,
                last=self.last,
                experiments=[e.record() for e in self.experiments],
            ),
        )

    def restore(self):
        try:
            data = json.loads(SETTINGS.read_text("utf-8"))
            self.automatic_output = data.get("automatic_output", True)
            if not self.automatic_output:
                self.output_base = Path(data.get("output_base") or data["output"])
                self.set_output(available_output(self.output_base))
            if data.get("method") in METHODS:
                self.method.set(data["method"])
            if data.get("inference") in INFERENCE_OPTIONS:
                self.inference.set(data["inference"])
            if data.get("geometry"):
                self.root.geometry(data["geometry"])
            desired_language = self._requested_language or data.get("language", "zh")
            if desired_language in ("zh", "en") and desired_language != self.lang:
                self.toggle_language()
            self.last = data.get("last")
            self.experiments = [
                Experiment(**x)
                for x in data.get("experiments", [])
                if all(Path(p).is_file() for p in x["files"])
            ]
            self.inputs = [p for e in self.experiments for p in e.files]
            self.render_catalog()
            self.refresh_output()
            if self.experiments:
                self.log("已恢复上次组别与显示顺序；开始分析前会重新核查文件。")
        except (OSError, ValueError, TypeError):
            pass

    def close(self):
        if self.running and (
            not messagebox.askyesno("正在计算", "关闭将中断计算，是否关闭？", parent=self.root)
        ):
            return
        self.persist()
        self.root.destroy()


def launch(execute, smoke=None, language=None):
    Desktop(execute, smoke, language)
