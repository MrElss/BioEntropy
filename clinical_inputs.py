from dataclasses import dataclass
from pathlib import Path

WEEKS = (0, 6, 12, 18, 24, 31, 38, 45, 52)


class InputSelectionError(ValueError):
    pass


@dataclass(frozen=True)
class InputBundle:
    experiment: str
    folder: Path
    normal_file: Path
    source_file: Path | None = None
    visit_files: tuple = ()

    def paths(self):
        return ([self.source_file] if self.source_file else list(self.visit_files)) + [
            self.normal_file
        ]

    def record(self):
        return {
            "experiment": self.experiment,
            "folder": str(self.folder),
            "normal_file": str(self.normal_file),
            "source_file": str(self.source_file) if self.source_file else None,
            "visit_files": [str(p) for p in self.visit_files],
        }
