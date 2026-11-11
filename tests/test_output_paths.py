import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from concurrent.futures import ThreadPoolExecutor
from output_paths import available_output, reserve_output, experiment_folder_name, desktop_directory


class OutputPathTests(unittest.TestCase):
    def test_selection_only_previews_path(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "105"
            self.assertEqual(available_output(path), path)
            self.assertFalse(path.exists())

    def test_existing_files_and_empty_directories_are_preserved(self):
        with tempfile.TemporaryDirectory() as root:
            base = Path(root) / "105"
            base.mkdir()
            (base / "原结果.txt").write_text("keep", encoding="utf-8")
            base.with_name("105-2").mkdir()
            base.with_name("105-3").write_text("keep file", encoding="utf-8")
            self.assertEqual(reserve_output(base).name, "105-4")
            self.assertEqual((base / "原结果.txt").read_text("utf-8"), "keep")
            self.assertEqual(base.with_name("105-3").read_text("utf-8"), "keep file")

    def test_repeated_and_concurrent_runs_get_separate_folders(self):
        with tempfile.TemporaryDirectory() as root:
            base = Path(root) / "KK鼠"
            with ThreadPoolExecutor(max_workers=4) as pool:
                paths = list(pool.map(lambda _: reserve_output(base), range(6)))
            self.assertEqual(len(set(paths)), 6)
            self.assertEqual({p.name for p in paths}, {"KK鼠", *[f"KK鼠-{i}" for i in range(2, 7)]})
            self.assertTrue(all((p.is_dir() for p in paths)))

    def test_collision_after_preview_is_checked_again(self):
        with tempfile.TemporaryDirectory() as root:
            base = Path(root) / "地鼠"
            self.assertEqual(available_output(base), base)
            base.mkdir()
            self.assertEqual(reserve_output(base).name, "地鼠-2")

    def test_multiple_active_experiments_and_safe_names(self):
        experiments = [
            SimpleNamespace(name=n, enabled=enabled)
            for n, enabled in [("105", True), ("106", True), ("KK鼠", False)]
        ]
        self.assertEqual(experiment_folder_name(experiments), "105+106")
        self.assertEqual(
            experiment_folder_name([SimpleNamespace(name="../CON?", enabled=True)]), "_CON_"
        )
        self.assertEqual(
            experiment_folder_name([SimpleNamespace(name="CON", enabled=True)]), "实验_CON"
        )

    def test_custom_parent_and_desktop(self):
        self.assertTrue(desktop_directory().is_absolute())
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "指定位置" / "DB鼠"
            self.assertEqual(reserve_output(path), path.resolve())
            self.assertTrue(path.is_dir())


if __name__ == "__main__":
    unittest.main()
