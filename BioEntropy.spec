from PyInstaller.utils.hooks import collect_submodules, collect_data_files

a = Analysis(
    ["BioEntropy.py"],
    pathex=[],
    binaries=[],
    datas=[
        ("BioEntropy.ico", "."),
        ("使用说明.md", "."),
        ("profiles.json", "."),
        ("clinical_entropy_scale.json", "."),
        ("clinical_patient_ids.py", "."),
        ("template_headers.json", "."),
    ]
    + collect_data_files("scipy.stats", includes=["_sobol_direction_numbers.npz"]),
    hiddenimports=collect_submodules("sklearn.utils"),
    hookspath=[],
    hooksconfig={"matplotlib": {"backends": ["Agg"]}},
    runtime_hooks=[],
    excludes=["pytest", "IPython", "jupyter", "notebook", "statsmodels", "sympy"],
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name="BioEntropy_V1.0.2",
    version="version_info.txt",
    icon="BioEntropy.ico",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    runtime_tmpdir=None,
)
