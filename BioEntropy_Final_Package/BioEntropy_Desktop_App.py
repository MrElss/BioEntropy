
from __future__ import annotations

import json
import os
import queue
import shutil
import subprocess
import sys
import threading
import traceback
from datetime import datetime
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk

from bioentropy_core_generalized import (
    WINSOR_Q_DEFAULT,
    analyze_all,
    build_excel_template,
    build_input_guide_html,
    build_normal_range_template,
    build_normal_sample_template,
    default_group_display_order,
    generate_output_bundle,
    get_experiments,
    get_groups,
    load_input_bundle,
    suggest_roles,
    suggest_reference_group,
    NORMAL_REF_NONE,
    UPLOADED_TRUE_NORMAL_LABEL,
)
from gvalue_recompute import analyze_gvalues

APP_NAME = "BioEntropy 桌面版"
# Keep in sync with SOFTWARE_VERSION in bioentropy_constants.
APP_VERSION = "V1.6"


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def unique_path(path: Path) -> Path:
    if not path.exists():
        return path
    stem = path.stem
    suffix = path.suffix
    parent = path.parent
    i = 2
    while True:
        candidate = parent / f"{stem}_{i}{suffix}"
        if not candidate.exists():
            return candidate
        i += 1


class BioEntropyDesktopApp:
    def __init__(self) -> None:
        self.root = tk.Tk()
        self.root.title(f"{APP_NAME} - {APP_VERSION}")
        self._set_window_icon()
        self.root.geometry("1260x930")
        # Content area scrolls, so a short window still reaches the Run button.
        self.root.minsize(1000, 560)

        self.file_paths: list[str] = []
        self.role_vars: dict[str, dict[str, tk.StringVar]] = {}
        self.group_order_listboxes: dict[str, tk.Listbox] = {}
        self.last_result_dir: Path | None = None
        self.last_zip_path: Path | None = None
        self.current_package_name: str | None = None
        self.is_running = False
        self.log_queue: queue.Queue[str] = queue.Queue()
        self.intro_expanded = False
        self.mixed_overview = None
        self.range_overview = None
        self.normal_overview = None
        self.range_map = {}
        self.true_normal_map = {}

        # Output folder is created on demand when an analysis runs, not at launch.
        self.output_dir = self.default_output_dir()

        self.status_var = tk.StringVar(value="请选择 Excel 文件。")
        self.output_var = tk.StringVar(value=str(self.output_dir))
        self.winsor_var = tk.StringVar(value=str(WINSOR_Q_DEFAULT))
        self.intro_button_var = tk.StringVar(value="展开软件简介")

        self.lang = "zh"
        self._i18n: list = []
        self._log_pristine = True
        self._apply_saved_settings()
        self._build_ui()
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)
        self._poll_log_queue()

    def _settings_path(self) -> Path:
        return Path.home() / ".bioentropy_desktop.json"

    def _apply_saved_settings(self) -> None:
        """Restore preferences from the last session (best effort)."""
        try:
            data = json.loads(self._settings_path().read_text(encoding="utf-8"))
        except Exception:
            return
        try:
            out = data.get("output_dir")
            if out and Path(out).parent.exists():
                # Restore the saved path only; do not create it at launch.
                self.output_dir = Path(out)
                self.output_var.set(str(self.output_dir))
        except Exception:
            pass
        try:
            wq = data.get("winsor_q")
            if wq is not None:
                float(wq)  # validate it parses
                self.winsor_var.set(str(wq))
        except Exception:
            pass
        try:
            if data.get("language") in ("zh", "en"):
                self.lang = data["language"]
        except Exception:
            pass
        try:
            geo = data.get("geometry")
            if isinstance(geo, str) and "x" in geo:
                w, h = (int(v) for v in geo.split("+")[0].split("x"))
                sw, sh = self.root.winfo_screenwidth(), self.root.winfo_screenheight()
                if 600 <= w <= sw and 400 <= h <= sh:
                    self.root.geometry(geo)
        except Exception:
            pass

    def _persist_settings(self) -> None:
        try:
            self._settings_path().write_text(
                json.dumps({
                    "output_dir": self.output_var.get(),
                    "winsor_q": self.winsor_var.get(),
                    "language": self.lang,
                    "geometry": self.root.geometry(),
                }, ensure_ascii=False, indent=2),
                encoding="utf-8",
            )
        except Exception:
            pass

    def _on_close(self) -> None:
        self._persist_settings()
        for attr in ("_motif_after", "_poll_after"):
            after_id = getattr(self, attr, None)
            if after_id is not None:
                try:
                    self.root.after_cancel(after_id)
                except Exception:
                    pass
                setattr(self, attr, None)
        try:
            self.root.destroy()
        except Exception:
            pass

    def default_output_dir(self) -> Path:
        # Desktop if it exists, else the working directory. No dedicated output
        # folder: each run writes its own timestamped BioEntropy_Results_* folder.
        desktop = Path.home() / "Desktop"
        if desktop.exists():
            return desktop
        return Path.cwd()

    # --------------------------------------------------------------- i18n
    def _t(self, key: str) -> str:
        pair = self.TR.get(key)
        if not pair:
            return key
        return pair[0] if self.lang == "zh" else pair[1]

    def _reg(self, widget, key):
        """Register a text-bearing widget for live language switching."""
        self._i18n.append((widget, key))
        try:
            widget.configure(text=self._t(key))
        except Exception:
            pass
        return widget

    def _toggle_language(self) -> None:
        self.lang = "en" if self.lang == "zh" else "zh"
        self._apply_language()

    def _apply_language(self) -> None:
        self.root.title(self._t("title"))
        for widget, key in self._i18n:
            try:
                widget.configure(text=self._t(key))
            except Exception:
                pass
        for col, key in getattr(self, "_tree_cols", []):
            try:
                self.file_tree.heading(col, text=self._t(key))
            except Exception:
                pass
        if hasattr(self, "lang_btn"):
            self.lang_btn.configure(text="EN" if self.lang == "zh" else "中文")
        self.intro_button_var.set(self._t("intro_hide" if self.intro_expanded else "intro_show"))
        self._set_intro_text()
        if self._log_pristine:
            self._render_welcome()

    def _set_intro_text(self) -> None:
        self.intro_text.configure(state="normal")
        self.intro_text.delete("1.0", tk.END)
        self.intro_text.insert("1.0", self._t("intro_body"))
        self.intro_text.configure(state="disabled")

    def _render_welcome(self) -> None:
        self.log_text.configure(state="normal")
        self.log_text.delete("1.0", tk.END)
        self.log_text.insert(tk.END, self._t("welcome"))
        self.log_text.configure(state="disabled")

    # ------------------------------------------------------ visual theme
    # Academic palette: a deep scholarly navy with a calm teal "entropy" accent
    # (matches the NPG-style figure palette used in the reports).
    NAVY = "#16304f"
    INK = "#1f2933"
    ACCENT = "#00A087"
    ACCENT_DK = "#0b7d6e"
    BLUE = "#3C5488"
    BG = "#eef1f6"
    CARD = "#ffffff"
    MUTED = "#6b7280"
    BORDER = "#d6dce6"
    SUBTLE = "#b9c4dc"

    def _setup_style(self) -> None:
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except Exception:
            pass
        self.root.configure(bg=self.BG)
        f = "Microsoft YaHei UI"
        style.configure(".", font=(f, 10))
        style.configure("TFrame", background=self.BG)
        style.configure("Card.TFrame", background=self.CARD)
        style.configure("Header.TFrame", background=self.NAVY)
        style.configure("TLabel", background=self.BG, foreground=self.INK)
        style.configure("Card.TLabel", background=self.CARD, foreground=self.INK)
        style.configure("Hint.TLabel", background=self.CARD, foreground=self.MUTED, font=(f, 9))
        style.configure("Status.TLabel", background=self.BG, foreground=self.BLUE, font=(f, 10, "bold"))
        style.configure("HeaderTitle.TLabel", background=self.NAVY, foreground="#ffffff", font=(f, 23, "bold"))
        style.configure("HeaderSub.TLabel", background=self.NAVY, foreground=self.SUBTLE, font=(f, 10))
        style.configure("Card.TLabelframe", background=self.CARD, bordercolor=self.BORDER,
                        relief="solid", borderwidth=1)
        style.configure("Card.TLabelframe.Label", background=self.CARD, foreground=self.NAVY,
                        font=(f, 11, "bold"))
        style.configure("TButton", font=(f, 9), padding=(11, 6), relief="flat",
                        background="#eaeff7", foreground=self.INK, bordercolor=self.BORDER)
        style.map("TButton",
                  background=[("active", "#dce5f2"), ("disabled", "#f0f2f6")],
                  foreground=[("disabled", "#a6adba")])
        style.configure("Accent.TButton", font=(f, 10, "bold"), padding=(18, 9),
                        background=self.ACCENT, foreground="#ffffff")
        style.map("Accent.TButton",
                  background=[("active", self.ACCENT_DK), ("disabled", "#9fccc3")])
        style.configure("Lang.TButton", font=(f, 9, "bold"), padding=(10, 4),
                        background=self.NAVY, foreground="#ffffff", bordercolor=self.SUBTLE)
        style.map("Lang.TButton", background=[("active", "#2c3f63")])
        style.configure("TEntry", fieldbackground="#ffffff", bordercolor=self.BORDER)
        style.configure("Treeview", background=self.CARD, fieldbackground=self.CARD,
                        foreground=self.INK, rowheight=26, bordercolor=self.BORDER)
        style.configure("Treeview.Heading", background=self.BLUE, foreground="#ffffff",
                        font=(f, 10, "bold"), relief="flat")
        style.map("Treeview.Heading", background=[("active", self.NAVY)])
        style.configure("TProgressbar", background=self.ACCENT, troughcolor="#dce3ee",
                        bordercolor="#dce3ee", lightcolor=self.ACCENT, darkcolor=self.ACCENT)

    def _init_motif(self, canvas: tk.Canvas) -> None:
        """A gently "breathing" Gaussian point-cloud — a live nod to system entropy:
        the dispersion of the cloud expands and contracts, evoking entropy change."""
        import random
        rnd = random.Random(7)
        self._motif_canvas = canvas
        # each particle: unit-Gaussian base offset + a random phase for jitter
        self._motif_pts = [(rnd.gauss(0, 1), rnd.gauss(0, 1), rnd.random()) for _ in range(48)]
        self._motif_phase = 0.0
        self._animate_motif()

    def _animate_motif(self) -> None:
        import math
        c = getattr(self, "_motif_canvas", None)
        if c is None or not c.winfo_exists():
            return
        w = int(c["width"]); h = int(c["height"]); cx, cy = w / 2.0, h / 2.0
        self._motif_phase += 0.045
        breathe = 1.0 + 0.20 * math.sin(self._motif_phase)
        sx, sy = w * 0.14 * breathe, h * 0.22 * breathe
        c.delete("all")
        for gx, gy, ph in self._motif_pts:
            jit = 1.0 + 0.10 * math.sin(self._motif_phase * 1.7 + ph * 6.283)
            x, y = cx + gx * sx * jit, cy + gy * sy * jit
            d = gx * gx + gy * gy
            if d < 0.5:
                shade, r = "#9af0e4", 2.4
            elif d < 1.6:
                shade, r = "#5fc4ba", 1.9
            else:
                shade, r = "#3f7f86", 1.5
            c.create_oval(x - r, y - r, x + r, y + r, fill=shade, outline="")
        self._motif_after = c.after(60, self._animate_motif)

    def _set_window_icon(self) -> None:
        """Replace Tk's default feather icon with an entropy-themed icon: a
        Gaussian point cloud (the same 'system-entropy dispersion' motif used in
        the header) on a rounded navy tile. Rendered once at startup; failures
        are swallowed so a missing backend can never block launch."""
        try:
            import tempfile
            import random
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.patches import FancyBboxPatch

            rnd = random.Random(7)
            fig = plt.figure(figsize=(2.56, 2.56), dpi=100)
            ax = fig.add_axes([0, 0, 1, 1])
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis("off")
            ax.add_patch(FancyBboxPatch(
                (0.07, 0.07), 0.86, 0.86,
                boxstyle="round,pad=0,rounding_size=0.20",
                fc=self.NAVY, ec="none"))
            for _ in range(46):
                gx, gy = rnd.gauss(0, 1), rnd.gauss(0, 1)
                x, y = 0.5 + gx * 0.125, 0.5 + gy * 0.125
                d = gx * gx + gy * gy
                if d < 0.5:
                    col, size = "#9af0e4", 62
                elif d < 1.6:
                    col, size = "#5fc4ba", 44
                else:
                    col, size = "#3f7f86", 30
                ax.scatter([x], [y], s=size, c=col, edgecolors="none", zorder=3)
            icon_path = Path(tempfile.gettempdir()) / "bioentropy_app_icon.png"
            fig.savefig(icon_path, transparent=True)
            plt.close(fig)
            self._app_icon = tk.PhotoImage(file=str(icon_path))
            self.root.iconphoto(True, self._app_icon)
        except Exception:
            pass

    # ----------------------------------------------------------- layout
    def _build_ui(self) -> None:
        self.TR = {
            "title": ("BioEntropy 桌面版", "BioEntropy Desktop"),
            "tagline": ("多指标系统熵 · 疾病偏离度分析框架",
                        "Multi-index system entropy · disease-deviation framework"),
            "version": ("版本", "Version"),
            "about": ("软件简介", "About"),
            "about_hint": ("简要了解本软件的用途、分析框架和适用场景。",
                           "A quick look at the purpose, analysis framework and use cases."),
            "intro_show": ("展开软件简介", "Show overview"),
            "intro_hide": ("收起软件简介", "Hide overview"),
            "intro_body": (
                "BioEntropy 用于多分组、多时间点生物指标数据的疾病偏离度—熵值分析。软件可批量读取规范 Excel 文件，"
                "自动完成共同指标识别、稳健预处理、Entropy 计算、异常值质量控制、时间趋势总结和结果导出。\n\n"
                "如果上传真实的独立 Normal/健康样本文件（例如“实验号-Normal.xlsx”），软件会把这些样本当作真正的健康参考云，"
                "使用稳健封顶 Normal-reference distance 计算 NMD，并保留原始协方差距离用于追溯。\n\n"
                "如果没有真实 Normal 样本但有外部正常区间，可上传“实验号-NormalRange.xlsx”，软件会额外输出 NRBS / NRPS。"
                "结果包括 Excel 工作簿、9 个主图/主文件、更多辅助图、自动报告和结果解释报告，统一输出为一个结果文件夹。",
                "BioEntropy analyses disease deviation and system entropy across multi-group, "
                "multi-timepoint biomarker data. It batch-reads standardised Excel files and "
                "performs common-feature detection, robust preprocessing, entropy estimation, "
                "outlier QC, time-trend summaries and result export.\n\n"
                "If you upload real independent Normal/healthy-sample files (e.g. "
                "\"<exp>-Normal.xlsx\"), they are treated as a true healthy reference cloud and "
                "used to compute a robustly-capped Normal-reference distance (NMD), keeping the "
                "raw covariance distance for traceability.\n\n"
                "Without real Normal samples but with an external normal range, upload "
                "\"<exp>-NormalRange.xlsx\" to also get NRBS / NRPS. Output includes an Excel "
                "workbook, 9 main figures/files, extra auxiliary figures, automatic reports and an "
                "explanation report, delivered as a single results folder."),
            "sec_files": ("1）文件、模板与输出目录", "1) Files, templates & output"),
            "grp_actions": ("数据与结果", "Data & results"),
            "grp_templates": ("模板与说明", "Templates & guides"),
            "btn_select": ("选择 Excel 文件", "Select Excel files"),
            "btn_clear": ("清空文件", "Clear files"),
            "btn_outdir": ("选择输出目录", "Choose output folder"),
            "btn_openout": ("打开结果文件夹", "Open results folder"),
            "btn_guide": ("下载命名与填写说明", "Naming & format guide"),
            "btn_tpl_data": ("下载实验数据模板", "Data template"),
            "btn_tpl_normal": ("下载 Normal 样本模板", "Normal-sample template"),
            "btn_tpl_range": ("下载 Normal 范围模板", "Normal-range template"),
            "btn_outlier": ("查看异常值算法说明", "Outlier method notes"),
            "btn_gvalues": ("Hedges'g值计算", "Hedges' g calculation"),
            "lbl_outdir": ("输出目录：", "Output folder:"),
            "lbl_winsor": ("稳健截尾比例 winsor_q：", "Robust winsorization winsor_q:"),
            "hint_winsor": ("设置原则：仅在少数极端值明显影响结果稳定性时再小幅上调；数值越大，截尾越强。建议先用默认值并比较结果是否稳定。",
                            "Only raise this slightly when a few extreme values clearly affect "
                            "stability; larger = stronger trimming. Start with the default and "
                            "compare for stability."),
            "sec_selected": ("2）已选择文件", "2) Selected files"),
            "col_file": ("文件名", "File"),
            "col_type": ("类型", "Type"),
            "col_exp": ("实验号", "Experiment"),
            "col_week": ("周数/说明", "Week / note"),
            "col_records": ("记录数", "Records"),
            "col_note": ("说明", "Note"),
            "sec_roles": ("3）每个实验的治疗组 / 比较组 / 健康参考组",
                          "3) Treatment / comparator / healthy-reference per experiment"),
            "sec_run": ("4）开始分析", "4) Run analysis"),
            "btn_run": ("开始分析并生成结果", "Run analysis & generate results"),
            "status_select": ("请选择 Excel 文件。", "Please select Excel files."),
            "sec_log": ("运行日志", "Run log"),
            "welcome": (
                "程序已启动。\n"
                "请点击“选择 Excel 文件”，可同时选择实验时间点文件、独立 Normal 样本文件和外部 Normal 范围文件。\n"
                "实验时间点文件命名：实验号-周数W.xlsx。\n"
                "独立 Normal 样本文件命名：实验号-Normal.xlsx（可上传多个，软件会自动合并）。\n"
                "外部 Normal 范围文件命名：实验号-NormalRange.xlsx。\n"
                "上传真实 Normal 样本会计算 NMD；只有外部区间文件则计算 NRBS / NRPS。\n",
                "Application started.\n"
                "Click \"Select Excel files\" to choose timepoint files, independent "
                "Normal-sample files and external Normal-range files together.\n"
                "Timepoint file name: <exp>-<week>W.xlsx.\n"
                "Independent Normal-sample file name: <exp>-Normal.xlsx (multiple allowed, "
                "auto-merged).\n"
                "External Normal-range file name: <exp>-NormalRange.xlsx.\n"
                "Real Normal samples compute NMD; external ranges only compute NRBS / NRPS.\n"),
        }

        self._setup_style()
        self.root.columnconfigure(0, weight=1)
        # Row 0 = fixed header; row 1 = the scrollable content area (expands).
        self.root.rowconfigure(1, weight=1)
        f = "Microsoft YaHei UI"

        # ---- header band (navy) with title, tagline, animated entropy motif
        header = ttk.Frame(self.root, style="Header.TFrame", padding=(22, 14, 22, 0))
        header.grid(row=0, column=0, sticky="ew")
        header.columnconfigure(1, weight=1)

        motif = tk.Canvas(header, width=210, height=72, bg=self.NAVY, highlightthickness=0, bd=0)
        motif.grid(row=0, column=0, rowspan=2, padx=(0, 18))
        self._init_motif(motif)

        ttk.Label(header, text="BioEntropy", style="HeaderTitle.TLabel").grid(row=0, column=1, sticky="sw")
        self._reg(ttk.Label(header, style="HeaderSub.TLabel"), "tagline").grid(row=1, column=1, sticky="nw")

        right = ttk.Frame(header, style="Header.TFrame")
        right.grid(row=0, column=2, rowspan=2, sticky="e")
        ttk.Label(header, text=APP_VERSION, style="HeaderSub.TLabel").grid(row=0, column=3, rowspan=2, sticky="e", padx=(12, 0))
        self.lang_btn = ttk.Button(right, text="EN", style="Lang.TButton", width=6,
                                   cursor="hand2", command=self._toggle_language)
        self.lang_btn.pack(side="right")
        # teal accent rule under the header
        tk.Frame(header, bg=self.ACCENT, height=3).grid(row=2, column=0, columnspan=4, sticky="ew", pady=(12, 0))

        # ---- scrollable content area: everything below the fixed header lives
        # inside a canvas so a vertical scrollbar appears on the right and the
        # Run button stays reachable on short / low-resolution screens.
        scroll_outer = ttk.Frame(self.root)
        scroll_outer.grid(row=1, column=0, sticky="nsew")
        scroll_outer.columnconfigure(0, weight=1)
        scroll_outer.rowconfigure(0, weight=1)
        self._scroll_canvas = tk.Canvas(scroll_outer, bg=self.BG, highlightthickness=0, bd=0)
        self._scroll_canvas.grid(row=0, column=0, sticky="nsew")
        vbar = ttk.Scrollbar(scroll_outer, orient="vertical", command=self._scroll_canvas.yview)
        vbar.grid(row=0, column=1, sticky="ns")
        self._scroll_canvas.configure(yscrollcommand=vbar.set)
        body = ttk.Frame(self._scroll_canvas)
        body.columnconfigure(0, weight=1)
        self._scroll_body_id = self._scroll_canvas.create_window((0, 0), window=body, anchor="nw")
        body.bind(
            "<Configure>",
            lambda e: self._scroll_canvas.configure(scrollregion=self._scroll_canvas.bbox("all")),
        )
        self._scroll_canvas.bind(
            "<Configure>",
            lambda e: self._scroll_canvas.itemconfigure(self._scroll_body_id, width=e.width),
        )

        def _on_mousewheel(event):
            self._scroll_canvas.yview_scroll(int(-event.delta / 120), "units")

        self._scroll_canvas.bind_all("<MouseWheel>", _on_mousewheel)

        # ---- about / intro
        intro_wrap = self._reg(ttk.LabelFrame(body, style="Card.TLabelframe", padding=10), "about")
        intro_wrap.grid(row=0, column=0, sticky="ew", padx=16, pady=(12, 8))
        intro_wrap.columnconfigure(0, weight=1)
        intro_top = ttk.Frame(intro_wrap, style="Card.TFrame")
        intro_top.grid(row=0, column=0, sticky="ew")
        intro_top.columnconfigure(0, weight=1)
        self._reg(ttk.Label(intro_top, style="Hint.TLabel"), "about_hint").grid(row=0, column=0, sticky="w")
        ttk.Button(intro_top, textvariable=self.intro_button_var, cursor="hand2", command=self.toggle_intro).grid(row=0, column=1, sticky="e")
        self.intro_body = ttk.Frame(intro_wrap, style="Card.TFrame")
        self.intro_text = tk.Text(self.intro_body, height=9, wrap=tk.WORD, relief=tk.FLAT,
                                  bg=self.CARD, fg=self.INK, font=(f, 10), padx=4, pady=4)
        self.intro_text.pack(fill="both", expand=True)

        # ---- section 1: files, templates, output
        controls = self._reg(ttk.LabelFrame(body, style="Card.TLabelframe", padding=12), "sec_files")
        controls.grid(row=1, column=0, sticky="ew", padx=16, pady=(0, 8))
        controls.columnconfigure(0, weight=1)

        def button_row(specs):
            bar = ttk.Frame(controls, style="Card.TFrame")
            for key, cmd in specs:
                b = ttk.Button(bar, command=cmd, cursor="hand2")
                self._reg(b, key)
                b.pack(side="left", padx=(0, 8), pady=3)
            return bar

        button_row([
            ("btn_select", self.select_files), ("btn_clear", self.clear_files),
            ("btn_outdir", self.choose_output_dir), ("btn_openout", self.open_last_result_dir),
        ]).grid(row=0, column=0, sticky="w")
        templates_bar = ttk.Frame(controls, style="Card.TFrame")
        for key, cmd in [
            ("btn_guide", self.export_input_guide), ("btn_tpl_data", self.export_excel_template),
            ("btn_tpl_normal", self.export_normal_sample_template),
            ("btn_tpl_range", self.export_normal_range_template),
            ("btn_outlier", self.show_outlier_rules),
        ]:
            self._reg(ttk.Button(templates_bar, command=cmd, cursor="hand2"), key).pack(side="left", padx=(0, 8), pady=3)
        self.gvalues_btn = ttk.Button(templates_bar, command=self.start_gvalues_analysis, cursor="hand2")
        self._reg(self.gvalues_btn, "btn_gvalues")
        self.gvalues_btn.pack(side="left", padx=(16, 8), pady=3)
        templates_bar.grid(row=1, column=0, sticky="w")

        fields = ttk.Frame(controls, style="Card.TFrame")
        fields.grid(row=2, column=0, sticky="ew", pady=(10, 0))
        fields.columnconfigure(1, weight=1)
        self._reg(ttk.Label(fields, style="Card.TLabel"), "lbl_outdir").grid(row=0, column=0, sticky="w", pady=4)
        ttk.Entry(fields, textvariable=self.output_var).grid(row=0, column=1, columnspan=3, sticky="ew", pady=4)
        self._reg(ttk.Label(fields, style="Card.TLabel"), "lbl_winsor").grid(row=1, column=0, sticky="w", pady=4)
        ttk.Entry(fields, textvariable=self.winsor_var, width=10).grid(row=1, column=1, sticky="w", pady=4)
        self._reg(ttk.Label(fields, style="Hint.TLabel", wraplength=900), "hint_winsor").grid(row=2, column=1, columnspan=3, sticky="w")

        # ---- section 2: selected files
        files_frame = self._reg(ttk.LabelFrame(body, style="Card.TLabelframe", padding=10), "sec_selected")
        files_frame.grid(row=2, column=0, sticky="ew", padx=16, pady=(0, 8))
        files_frame.columnconfigure(0, weight=1)
        files_frame.rowconfigure(0, weight=1)
        self._tree_cols = [("c_file", "col_file"), ("c_type", "col_type"), ("c_exp", "col_exp"),
                           ("c_week", "col_week"), ("c_records", "col_records"), ("c_note", "col_note")]
        col_ids = [c for c, _ in self._tree_cols]
        self.file_tree = ttk.Treeview(files_frame, columns=col_ids, show="headings", height=8)
        widths = {"c_file": 320, "c_type": 100, "c_exp": 120, "c_week": 120, "c_records": 90, "c_note": 280}
        for col, key in self._tree_cols:
            self.file_tree.heading(col, text=self._t(key))
            self.file_tree.column(col, width=widths[col], anchor="center")
        self.file_tree.grid(row=0, column=0, sticky="nsew")
        yscroll = ttk.Scrollbar(files_frame, orient="vertical", command=self.file_tree.yview)
        yscroll.grid(row=0, column=1, sticky="ns")
        self.file_tree.configure(yscrollcommand=yscroll.set)

        # ---- section 3: roles
        role_box = self._reg(ttk.LabelFrame(body, style="Card.TLabelframe", padding=10), "sec_roles")
        role_box.grid(row=3, column=0, sticky="ew", padx=16, pady=(0, 8))
        role_box.columnconfigure(0, weight=1)
        self.role_container = ttk.Frame(role_box, style="Card.TFrame")
        self.role_container.grid(row=0, column=0, sticky="ew")
        self.role_container.columnconfigure(0, weight=1)

        # ---- section 4: run
        action_frame = self._reg(ttk.LabelFrame(body, style="Card.TLabelframe", padding=12), "sec_run")
        action_frame.grid(row=4, column=0, sticky="ew", padx=16, pady=(0, 8))
        action_frame.columnconfigure(1, weight=1)
        self.run_btn = ttk.Button(action_frame, style="Accent.TButton", cursor="hand2", command=self.start_analysis)
        self._reg(self.run_btn, "btn_run")
        self.run_btn.grid(row=0, column=0, padx=(0, 12))
        self.progress = ttk.Progressbar(action_frame, mode="indeterminate")
        self.progress.grid(row=0, column=1, sticky="ew", padx=(0, 12))
        ttk.Label(action_frame, textvariable=self.status_var, style="Status.TLabel").grid(row=0, column=2, sticky="e")

        # ---- run log
        log_frame = self._reg(ttk.LabelFrame(body, style="Card.TLabelframe", padding=10), "sec_log")
        log_frame.grid(row=5, column=0, sticky="ew", padx=16, pady=(0, 16))
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)
        self.log_text = scrolledtext.ScrolledText(log_frame, wrap=tk.WORD, font=("Consolas", 10),
                                                  height=10,
                                                  bg="#0f1722", fg="#d6e2f0", insertbackground="#d6e2f0",
                                                  relief=tk.FLAT)
        self.log_text.grid(row=0, column=0, sticky="nsew")

        self.status_var.set(self._t("status_select"))
        self._apply_language()

    def toggle_intro(self) -> None:
        self.intro_expanded = not self.intro_expanded
        if self.intro_expanded:
            self.intro_body.grid(row=1, column=0, sticky="ew", pady=(8, 0))
        else:
            self.intro_body.grid_forget()
        self.intro_button_var.set(self._t("intro_hide" if self.intro_expanded else "intro_show"))

    def log(self, text: str) -> None:
        self.log_queue.put(text)

    def _bind_group_order_drag(self, listbox: tk.Listbox) -> None:
        def on_press(event):
            index = listbox.nearest(event.y)
            if 0 <= index < listbox.size():
                listbox._drag_index = index  # type: ignore[attr-defined]
                listbox.selection_clear(0, tk.END)
                listbox.selection_set(index)
            return "break"

        def on_motion(event):
            old_index = getattr(listbox, "_drag_index", None)
            if old_index is None or listbox.size() <= 1:
                return "break"
            new_index = listbox.nearest(event.y)
            if new_index == old_index or not (0 <= new_index < listbox.size()):
                return "break"
            value = listbox.get(old_index)
            listbox.delete(old_index)
            listbox.insert(new_index, value)
            listbox.selection_clear(0, tk.END)
            listbox.selection_set(new_index)
            listbox._drag_index = new_index  # type: ignore[attr-defined]
            return "break"

        listbox.bind("<Button-1>", on_press)
        listbox.bind("<B1-Motion>", on_motion)

    def _collect_group_order_map(self) -> dict[str, list[str]]:
        out: dict[str, list[str]] = {}
        for exp, listbox in self.group_order_listboxes.items():
            values = [str(listbox.get(i)).strip() for i in range(listbox.size()) if str(listbox.get(i)).strip()]
            if values:
                out[exp] = values
        return out

    def _poll_log_queue(self) -> None:
        try:
            while True:
                msg = self.log_queue.get_nowait()
                self._log_pristine = False
                self.log_text.configure(state="normal")
                timestamp = datetime.now().strftime("%H:%M:%S")
                self.log_text.insert(tk.END, f"[{timestamp}] {msg}\n")
                self.log_text.see(tk.END)
                self.log_text.configure(state="disabled")
        except queue.Empty:
            pass
        self._poll_after = self.root.after(150, self._poll_log_queue)

    def select_files(self) -> None:
        paths = filedialog.askopenfilenames(
            title="选择一个或多个 Excel 文件",
            filetypes=[("Excel 文件", "*.xlsx *.xls"), ("所有文件", "*.*")],
        )
        if not paths:
            return
        self.file_paths = list(dict.fromkeys(self.file_paths + list(paths)))
        self.refresh_file_overview()

    def _sanitize_dirname(self, name: object) -> str:
        text = str(name).strip()
        for ch in '\\/:*?"<>|':
            text = text.replace(ch, "_")
        return text.strip()

    def _auto_set_output_dir_from_experiments(self, experiments) -> None:
        """When a dataset is loaded, point the output directory at
        ``Desktop/<experiment name>`` (falls back to the default location if there
        is no Desktop). The folder is created on demand when the analysis runs.
        The user can still override it afterwards via 'choose output directory'."""
        try:
            names = [self._sanitize_dirname(e) for e in experiments if str(e).strip()]
            names = [n for n in names if n]
            if not names:
                return
            label = names[0] if len(names) == 1 else "_".join(names[:4])
            base = Path.home() / "Desktop"
            if not base.exists():
                base = self.default_output_dir()
            target = base / label
            self.output_dir = target
            self.output_var.set(str(target))
            self.log(f"输出目录已根据实验名自动设置为：{target}")
        except Exception:
            pass

    def clear_files(self) -> None:
        if self.is_running:
            messagebox.showwarning("正在运行", "分析正在进行中，暂时不能清空文件。")
            return
        self.file_paths = []
        self.role_vars = {}
        self.group_order_listboxes = {}
        self.range_map = {}
        self.true_normal_map = {}
        self.mixed_overview = None
        self.range_overview = None
        self.normal_overview = None
        for row in self.file_tree.get_children():
            self.file_tree.delete(row)
        for child in self.role_container.winfo_children():
            child.destroy()
        self.group_order_listboxes = {}
        self.status_var.set("已清空文件。")
        self.log("已清空文件列表。")

    def choose_output_dir(self) -> None:
        path = filedialog.askdirectory(title="选择结果输出目录", initialdir=self.output_var.get())
        if not path:
            return
        self.output_dir = ensure_dir(Path(path))
        self.output_var.set(str(self.output_dir))
        self.log(f"输出目录已设置为：{self.output_dir}")

    def export_input_guide(self) -> None:
        ensure_dir(Path(self.output_var.get()))
        out_path = unique_path(Path(self.output_var.get()) / "输入文件命名与填写格式说明.html")
        build_input_guide_html(out_path)
        self.log(f"已导出命名与填写说明：{out_path}")
        messagebox.showinfo("导出完成", f"命名与填写说明已保存到：\n{out_path}")
        self.open_path(out_path)

    def export_excel_template(self) -> None:
        ensure_dir(Path(self.output_var.get()))
        out_path = unique_path(Path(self.output_var.get()) / "BioEntropy_Excel_Template.xlsx")
        build_excel_template(out_path)
        self.log(f"已导出实验数据模板：{out_path}")
        messagebox.showinfo("导出完成", f"实验数据模板已保存到：\n{out_path}")
        self.open_path(out_path)

    def export_normal_sample_template(self) -> None:
        ensure_dir(Path(self.output_var.get()))
        out_path = unique_path(Path(self.output_var.get()) / "BioEntropy_NormalSamples_Template.xlsx")
        build_normal_sample_template(out_path)
        self.log(f"已导出 Normal 样本模板：{out_path}")
        messagebox.showinfo("导出完成", f"Normal 样本模板已保存到：\n{out_path}")
        self.open_path(out_path)

    def export_normal_range_template(self) -> None:
        ensure_dir(Path(self.output_var.get()))
        out_path = unique_path(Path(self.output_var.get()) / "BioEntropy_NormalRange_Template.xlsx")
        build_normal_range_template(out_path)
        self.log(f"已导出 Normal 范围模板：{out_path}")
        messagebox.showinfo("导出完成", f"Normal 范围模板已保存到：\n{out_path}")
        self.open_path(out_path)

    def show_outlier_rules(self) -> None:
        msg = (
            "异常值处理原则\n\n"
            "1. 软件默认不会直接删除整行样本。\n"
            "   当前默认流程是：稳健预处理 + 异常值标记 + 结果导出。\n\n"
            "2. winsor_q 的作用不是删样本，而是对每个指标的极端尾部做稳健截尾。\n"
            "   具体做法：在锚定周上拟合每个指标的分位数边界，再把过高/过低的极端值压回边界，"
            "减少单个极端值对均值、方差和熵的支配。\n\n"
            "3. 质量控制里的“共识异常值”来自 5 类算法联合标记：\n"
            "   - robust z-score：任一指标 |robust z| > 4\n"
            "   - IQR：任一指标落在 [Q1-1.5×IQR, Q3+1.5×IQR] 之外\n"
            "   - Isolation Forest\n"
            "   - Local Outlier Factor (LOF)\n"
            "   - Robust Mahalanobis distance（MinCovDet，阈值为卡方 0.999 分位）\n\n"
            "4. 共识规则：同一样本被 ≥2 种方法同时标记，才记为 consensus outlier。\n\n"
            "5. 主熵值 Entropy 会额外采用更严格的高置信规则：同一样本被 ≥4 种方法同时标记时，"
            "不删除整行，而是把该样本的多指标向量向本组稳健中心收缩 50%，降低单个极端样本对协方差熵的支配。"
            "原始未收缩熵会保留在 Entropy_raw 字段，便于追溯。\n\n"
            "6. 这些标记不直接解释为疗效，也不建议为了得到想要的趋势而手工删点。"
        )
        messagebox.showinfo("异常值算法说明", msg)

    def refresh_file_overview(self) -> None:
        for row in self.file_tree.get_children():
            self.file_tree.delete(row)
        for child in self.role_container.winfo_children():
            child.destroy()
        self.role_vars = {}
        self.group_order_listboxes = {}
        self.range_map = {}
        self.true_normal_map = {}
        self.mixed_overview = None
        self.range_overview = None
        self.normal_overview = None

        if not self.file_paths:
            self.status_var.set("请选择 Excel 文件。")
            return

        try:
            raw, range_refs, normal_refs, mixed_overview, data_overview, range_overview, normal_overview = load_input_bundle(self.file_paths)
            self.range_map = range_refs
            self.true_normal_map = normal_refs
            self.mixed_overview = mixed_overview
            self.range_overview = range_overview
            self.normal_overview = normal_overview
        except Exception as exc:
            self.log(f"读取文件失败：{exc}")
            messagebox.showerror("文件读取失败", str(exc))
            return

        for _, row in mixed_overview.iterrows():
            self.file_tree.insert(
                "", tk.END,
                values=(row["文件名"], row["类型"], row["实验号"], row["周数/说明"], row["记录数"], row["说明"])
            )

        experiments = get_experiments(raw)
        if not experiments:
            self.status_var.set("没有检测到有效实验。")
            return
        self._auto_set_output_dir_from_experiments(experiments)

        # Data-quality / entropy-reliability hint: the log-determinant covariance
        # entropy needs the per-group sample size n to exceed the feature count p.
        try:
            for exp in experiments:
                weeks = [k for k in raw if k[0] == exp]
                if not weeks:
                    continue
                df0 = raw[weeks[0]]
                p = len([c for c in df0.columns if c not in ("Group", "Sample ID")])
                min_n = min(int(raw[k].groupby("Group").size().min()) for k in weeks)
                ratio = (min_n / p) if p else float("nan")
                if not p:
                    continue
                if ratio < 2:
                    tag = "偏弱（n/p<2，协方差熵不稳定，建议优先参考 NMD / 参考态熵）"
                elif ratio < 5:
                    tag = "有限（2≤n/p<5）"
                else:
                    tag = "良好（n/p≥5）"
                self.log(f"数据质量 · 实验 {exp}：指标数 p={p}，最小组样本数 n={min_n}，n/p≈{ratio:.1f} → 熵估计可靠性{tag}")
        except Exception:
            pass

        ttk.Label(
            self.role_container,
            text="治疗组/比较组可手动修改；每个实验下方会自动列出组别，拖动即可调整图中显示顺序。只有当 Normal 就在实验分组里时，才需要手动改健康参考组。",
            foreground="#666666",
        ).grid(row=0, column=0, columnspan=8, sticky="w", pady=(0, 8))

        for i, exp in enumerate(experiments, start=1):
            base_row = (i - 1) * 2 + 1
            groups = get_groups(raw, exp)
            display_groups = default_group_display_order(groups)
            treat_default, comp_default = suggest_roles(groups)
            if comp_default == treat_default:
                for g in groups:
                    if g != treat_default:
                        comp_default = g
                        break
            ref_default = suggest_reference_group(groups)
            has_external = exp in range_refs and range_refs[exp] is not None and len(range_refs[exp]) > 0
            has_true_normal = exp in normal_refs and normal_refs[exp] is not None and len(normal_refs[exp]) > 0

            true_normal_names = []
            if self.normal_overview is not None and len(self.normal_overview) > 0:
                true_normal_names = self.normal_overview.loc[
                    self.normal_overview["实验号"].astype(str) == str(exp), "文件名"
                ].dropna().astype(str).tolist()

            range_names = []
            if self.range_overview is not None and len(self.range_overview) > 0:
                range_names = self.range_overview.loc[
                    self.range_overview["实验号"].astype(str) == str(exp), "文件名"
                ].dropna().astype(str).tolist()

            ref_options = []
            if has_true_normal:
                ref_options.append(UPLOADED_TRUE_NORMAL_LABEL)
            ref_options.append(NORMAL_REF_NONE)
            ref_options.extend(groups)

            if has_true_normal and has_external:
                ext_text = f"独立Normal：{' | '.join(true_normal_names)}；外部范围：{' | '.join(range_names)}"
            elif has_true_normal:
                ext_text = f"独立Normal：{' | '.join(true_normal_names)}"
            elif has_external:
                ext_text = f"外部Normal范围：{' | '.join(range_names)}"
            else:
                ext_text = "未检测到Normal参考文件"

            ttk.Label(self.role_container, text=f"实验 {exp}", font=("Microsoft YaHei UI", 10, "bold")).grid(row=base_row, column=0, sticky="w", padx=(0, 10), pady=6)
            ttk.Label(self.role_container, text="治疗组：").grid(row=base_row, column=1, sticky="e", padx=(0, 4), pady=6)
            ttk.Label(self.role_container, text="比较组：").grid(row=base_row, column=3, sticky="e", padx=(0, 4), pady=6)
            ttk.Label(self.role_container, text="健康参考组：").grid(row=base_row, column=5, sticky="e", padx=(0, 4), pady=6)

            treat_var = tk.StringVar(value=treat_default)
            comp_var = tk.StringVar(value=comp_default)
            ref_var = tk.StringVar(value=UPLOADED_TRUE_NORMAL_LABEL if has_true_normal else (ref_default if ref_default is not None else NORMAL_REF_NONE))
            ttk.Combobox(self.role_container, values=groups, textvariable=treat_var, state="readonly", width=18).grid(row=base_row, column=2, sticky="w", padx=(0, 12), pady=6)
            ttk.Combobox(self.role_container, values=groups, textvariable=comp_var, state="readonly", width=18).grid(row=base_row, column=4, sticky="w", padx=(0, 12), pady=6)
            ttk.Combobox(self.role_container, values=ref_options, textvariable=ref_var, state="readonly", width=24).grid(row=base_row, column=6, sticky="w", padx=(0, 12), pady=6)
            ttk.Label(self.role_container, text=ext_text, foreground="#666666").grid(row=base_row, column=7, sticky="w", padx=(0, 6), pady=6)
            self.role_vars[exp] = {"treatment": treat_var, "comparator": comp_var, "reference": ref_var}

            ttk.Label(self.role_container, text="图中组别顺序：").grid(row=base_row + 1, column=1, sticky="ne", padx=(0, 4), pady=(0, 10))
            order_box = tk.Listbox(
                self.role_container,
                height=min(max(len(groups), 3), 7),
                width=34,
                exportselection=False,
                activestyle="dotbox",
            )
            for group in display_groups:
                order_box.insert(tk.END, group)
            order_box.grid(row=base_row + 1, column=2, columnspan=3, sticky="ew", padx=(0, 12), pady=(0, 10))
            self._bind_group_order_drag(order_box)
            self.group_order_listboxes[exp] = order_box
            ttk.Label(
                self.role_container,
                text="按住组名上下拖动；该顺序会用于主图、辅助图和归一化展示表。",
                foreground="#666666",
            ).grid(row=base_row + 1, column=5, columnspan=3, sticky="w", padx=(0, 6), pady=(0, 10))

        summary = f"已载入 {len(self.file_paths)} 个文件，识别到 {len(experiments)} 个实验。"
        if self.normal_overview is not None and len(self.normal_overview) > 0:
            summary += f" 其中独立 Normal 样本文件 {len(self.normal_overview)} 个。"
        if range_overview is not None and len(range_overview) > 0:
            summary += f" 外部 Normal 范围文件 {len(range_overview)} 个。"
        self.status_var.set(summary)
        self.log(summary)

    def validate_before_run(self):
        if not self.file_paths:
            messagebox.showwarning("缺少文件", "请先选择 Excel 文件。")
            return None

        out_dir = Path(self.output_var.get()).expanduser()
        ensure_dir(out_dir)
        self.output_dir = out_dir

        try:
            winsor_q = float(self.winsor_var.get().strip())
        except Exception:
            messagebox.showerror("参数错误", "winsor_q 必须是数字。")
            return None

        if not (0 <= winsor_q < 0.5):
            messagebox.showerror("参数错误", "winsor_q 必须在 0 到 0.5 之间。通常建议使用较小的数值。")
            return None

        group_order_map = self._collect_group_order_map()
        role_map = {}
        for exp, vars_map in self.role_vars.items():
            treatment = vars_map["treatment"].get().strip()
            comparator = vars_map["comparator"].get().strip()
            reference = vars_map["reference"].get().strip() if "reference" in vars_map else NORMAL_REF_NONE
            if not treatment or not comparator:
                messagebox.showerror("分组设置错误", f"实验 {exp} 的治疗组或比较组为空。")
                return None
            if treatment == comparator:
                proceed = messagebox.askyesno("治疗组和比较组相同", f"实验 {exp} 的治疗组和比较组相同，通常不合理。是否仍然继续？")
                if not proceed:
                    return None
            if reference in {NORMAL_REF_NONE, UPLOADED_TRUE_NORMAL_LABEL}:
                reference = None
            elif reference in {treatment, comparator}:
                proceed = messagebox.askyesno(
                    "健康参考组提示",
                    f"实验 {exp} 的健康参考组当前设置为“{reference}”。\n\n"
                    "只有当这个分组是真正的 Normal/健康参考组时才建议这么设置。\n"
                    "如果它只是 placebo/control/模型对照，请点“否”返回修改。\n\n"
                    "是否继续？"
                )
                if not proceed:
                    return None
            role_map[exp] = {"treatment": treatment, "comparator": comparator, "reference": reference}
        return role_map, winsor_q, group_order_map

    def set_running(self, running: bool) -> None:
        self.is_running = running
        if running:
            self.run_btn.configure(state="disabled")
            if hasattr(self, "gvalues_btn"):
                self.gvalues_btn.configure(state="disabled")
            self.progress.start(12)
        else:
            self.run_btn.configure(state="normal")
            if hasattr(self, "gvalues_btn"):
                self.gvalues_btn.configure(state="normal")
            self.progress.stop()


    def _handle_existing_output_items(self) -> bool:
        """Warn when the selected output folder already contains previous BioEntropy result folders.

        Returns True if analysis may continue, False if the user cancels or path remains unsuitable.
        """
        try:
            self.output_dir = Path(self.output_var.get()).expanduser()
            ensure_dir(self.output_dir)
        except Exception as exc:
            messagebox.showerror("输出目录错误", f"无法访问输出目录：\n{exc}")
            return False

        existing = []
        for pattern in ("BioEntropy_Results*",):
            existing.extend([p for p in self.output_dir.glob(pattern) if p.exists()])
        # De-duplicate and ignore non-result working files
        existing = sorted(set(existing), key=lambda p: str(p).lower())
        if not existing:
            return True

        preview = "\n".join([f"- {p.name}" for p in existing[:8]])
        if len(existing) > 8:
            preview += f"\n... 以及 {len(existing) - 8} 个其它旧结果"

        resp = messagebox.askyesnocancel(
            "输出目录已有结果",
            "当前输出路径中已存在 BioEntropy 结果文件夹：\n\n"
            f"{preview}\n\n"
            "选择“是”：清理这些旧结果后继续。\n"
            "选择“否”：修改输出路径。\n"
            "选择“取消”：终止本次分析。"
        )
        if resp is None:
            return False
        if resp is False:
            self.choose_output_dir()
            try:
                self.output_dir = Path(self.output_var.get()).expanduser()
                ensure_dir(self.output_dir)
            except Exception as exc:
                messagebox.showerror("输出目录错误", f"无法访问新的输出目录：\n{exc}")
                return False
            still_existing = sorted(set(self.output_dir.glob("BioEntropy_Results*")), key=lambda p: str(p).lower())
            if still_existing:
                messagebox.showwarning(
                    "输出目录仍有旧结果",
                    "新的输出路径中仍存在 BioEntropy 结果文件夹。\n"
                    "请重新选择一个空目录，或选择“是”清理旧结果后继续。"
                )
                return False
            return True

        # resp is True: delete previous BioEntropy results in the selected output directory
        try:
            for p in existing:
                if p.is_dir():
                    shutil.rmtree(p)
                else:
                    p.unlink()
            self.log(f"已清理输出目录中的 {len(existing)} 个旧 BioEntropy 结果。")
            return True
        except Exception as exc:
            messagebox.showerror("清理旧结果失败", f"无法删除旧结果：\n{exc}\n\n请手动清理或更换输出目录。")
            return False


    def start_analysis(self) -> None:
        if self.is_running:
            return
        validated = self.validate_before_run()
        if validated is None:
            return
        role_map, winsor_q, group_order_map = validated

        if not self._handle_existing_output_items():
            return

        package_name = f"BioEntropy_Results_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        self.current_package_name = package_name
        explanation_options = {
            "group_order_map": group_order_map,
        }
        candidate_dir = self.output_dir / package_name

        if candidate_dir.exists():
            resp = messagebox.askyesnocancel(
                "输出目录冲突",
                "当前输出路径中已存在同名结果文件夹。\n\n"
                f"目录：\n{candidate_dir}\n\n"
                "选择“是”将覆盖同名结果；选择“否”先修改输出路径；选择“取消”终止本次分析。"
            )
            if resp is None:
                return
            if resp is False:
                self.choose_output_dir()
                candidate_dir = self.output_dir / package_name
                if candidate_dir.exists():
                    messagebox.showwarning("仍然冲突", "修改后的输出路径中仍存在同名结果，请再次调整路径或稍后重试。")
                    return
            else:
                try:
                    if candidate_dir.exists():
                        shutil.rmtree(candidate_dir)
                except Exception as exc:
                    messagebox.showerror("无法覆盖", f"删除旧结果失败：\n{exc}\n\n请手动清理后再试。")
                    return

        self.set_running(True)
        self.status_var.set("正在分析，请稍候……")
        self.log("开始分析。")
        if group_order_map:
            self.log("已应用检测到/拖拽后的组别显示顺序。")
        threading.Thread(target=self._analysis_worker, args=(package_name, role_map, winsor_q, explanation_options), daemon=True).start()

    def _analysis_worker(self, package_name, role_map, winsor_q, explanation_options=None) -> None:
        try:
            self.log("读取 Excel 文件……")
            raw, range_refs, normal_refs, mixed_overview, data_overview, range_overview, normal_overview = load_input_bundle(self.file_paths)
            self.log("开始运行疾病偏离度—熵值分析框架……")
            if normal_overview is not None and len(normal_overview) > 0:
                self.log(f"检测到 {len(normal_overview)} 个独立 Normal 样本文件，将计算 True-Normal robust distance (NMD)。")
            if range_overview is not None and len(range_overview) > 0:
                self.log(f"检测到 {len(range_overview)} 个外部 Normal 范围文件，将额外计算 NRBS / NRPS。")
            results = analyze_all(raw, role_map=role_map, winsor_q=winsor_q, external_range_map=range_refs, true_normal_map=normal_refs)
            results["tables"]["上传文件概览"] = mixed_overview
            if normal_overview is not None and len(normal_overview) > 0:
                results["tables"]["独立Normal样本文件"] = normal_overview
            if range_overview is not None and len(range_overview) > 0:
                results["tables"]["外部Normal范围文件"] = range_overview
            if (explanation_options or {}).get("group_order_map"):
                results["group_order_map"] = (explanation_options or {}).get("group_order_map")
            self.log("生成结果工作簿、报告、图片和解释报告……")
            bundle = generate_output_bundle(results, package_name=package_name, explanation_options=explanation_options)
            ensure_dir(self.output_dir)
            final_dir = self.output_dir / package_name
            if final_dir.exists():
                shutil.rmtree(final_dir)
            shutil.copytree(bundle["output_dir"], final_dir)
            self.last_result_dir = final_dir
            self.last_zip_path = None
            self.root.after(0, lambda: self._finish_success(final_dir))
        except Exception as exc:
            err = str(exc)
            tb = traceback.format_exc()
            self.log(f"分析失败：{err}")
            self.log(tb)
            # capture into a local: `exc` is unbound once the except block exits,
            # and this lambda runs later on the Tk event loop.
            self.root.after(0, lambda: self._finish_error(err))

    def _finish_success(self, result_dir: Path) -> None:
        self.set_running(False)
        self.status_var.set("分析完成。")
        self.log(f"分析完成。结果目录：{result_dir}")
        messagebox.showinfo(
            "分析完成",
            "结果已经生成。\n\n"
            "已额外生成 analysis_report.html、analysis_report.md 和 explanation.json。\n\n"
            f"结果目录：\n{result_dir}\n\n"
            "结果文件夹将自动打开。"
        )
        self.open_path(result_dir)

    def _finish_error(self, error_text: str) -> None:
        self.set_running(False)
        self.status_var.set("分析失败，请看日志。")
        messagebox.showerror(
            "分析失败",
            f"程序运行失败：\n{error_text}\n\n完整的错误信息已写入下方运行日志，可据此排查输入文件格式或分组设置。"
        )

    def start_gvalues_analysis(self) -> None:
        if self.is_running:
            return

        initial_dir = Path.cwd()
        for candidate in (
            Path.cwd() / "data" / "G值森林图数据",
            Path.cwd() / "data",
        ):
            if candidate.exists():
                initial_dir = candidate
                break

        gvalue_dir = filedialog.askdirectory(
            title="选择 Hedges'g 数据目录（可包含 1）..5） 任意子文件夹，按现有数据计算）",
            initialdir=str(initial_dir),
        )
        if not gvalue_dir:
            return

        try:
            out_dir = ensure_dir(Path(self.output_var.get()).expanduser())
            self.output_dir = out_dir
        except Exception as exc:
            messagebox.showerror("输出目录错误", f"无法访问输出目录：\n{exc}")
            return

        self.set_running(True)
        self.status_var.set("正在进行 Hedges'g 值计算……")
        self.log(f"开始 Hedges'g 值计算。输入目录：{gvalue_dir}")
        threading.Thread(target=self._gvalues_worker, args=(Path(gvalue_dir), out_dir), daemon=True).start()

    def _gvalues_worker(self, gvalue_dir: Path, out_dir: Path) -> None:
        try:
            result = analyze_gvalues(str(gvalue_dir), str(out_dir))
            final_dir = Path(result["output_dir"])
            self.last_result_dir = final_dir
            self.last_zip_path = None
            self.log(f"Hedges'g 值计算完成：{final_dir}")
            done = result.get("parts_done", [])
            if done:
                self.log("已计算的分析：" + "、".join(done))
            headline = result.get("headline", {})
            if "mahalanobis_pooled_random_g" in headline:
                self.log(f"马氏距离合并 G（随机效应）：{headline['mahalanobis_pooled_random_g']:.3f}")
            if "entropy_pooled_random_g" in headline:
                self.log(f"熵合并 G（随机效应）：{headline['entropy_pooled_random_g']:.3f}")
            self.log(f"已生成森林图 {result.get('figure_count', 0)} 张。")
            if result.get("warning_count", 0):
                self.log(f"有 {result.get('warning_count')} 条警告，详见各 report.md。")
            self.root.after(0, lambda: self._finish_gvalues_success(result))
        except Exception as exc:
            err = str(exc)
            tb = traceback.format_exc()
            self.log(f"Hedges'g 值计算失败：{err}")
            self.log(tb)
            self.root.after(0, lambda: self._finish_error(err))

    def _finish_gvalues_success(self, result: dict) -> None:
        self.set_running(False)
        result_dir = Path(result["output_dir"])
        self.status_var.set("Hedges'g 值计算完成。")
        headline = result.get("headline", {})
        done = result.get("parts_done", [])

        msg = "Hedges'g 值计算已经完成并生成森林图。\n\n"
        if done:
            msg += "已计算的分析（按现有数据自动识别）：\n" + "、".join(done) + "\n\n"
        msg += f"森林图数量：{result.get('figure_count', 0)}\n"
        if "mahalanobis_pooled_random_g" in headline:
            msg += f"马氏距离合并 G（随机效应）：{headline['mahalanobis_pooled_random_g']:.3f}\n"
        if "entropy_pooled_random_g" in headline:
            msg += f"熵合并 G（随机效应）：{headline['entropy_pooled_random_g']:.3f}\n"
        msg += f"\n结果目录：\n{result_dir}"
        msg += "\n\n结果文件夹将自动打开。"
        messagebox.showinfo("Hedges'g 值计算完成", msg)
        self.open_path(result_dir)

    def open_last_result_dir(self) -> None:
        if self.last_result_dir is None:
            messagebox.showinfo("还没有结果", "请先运行一次分析。")
            return
        self.open_path(self.last_result_dir)

    def open_path(self, path: Path) -> None:
        try:
            if sys.platform.startswith("win"):
                os.startfile(str(path))  # type: ignore[attr-defined]
            elif sys.platform == "darwin":
                subprocess.Popen(["open", str(path)])
            else:
                subprocess.Popen(["xdg-open", str(path)])
        except Exception as exc:
            messagebox.showwarning("打开失败", f"无法自动打开路径：{exc}\n路径：{path}")

    def run(self) -> None:
        self.root.mainloop()


if __name__ == "__main__":
    BioEntropyDesktopApp().run()
