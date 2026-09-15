import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]


class ContainerBuildTests(unittest.TestCase):
    def test_build_from_another_directory_and_push_only_when_requested(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            docker = root / "docker"
            docker.write_text(
                "#!/usr/bin/env python3\nimport json, os, sys\n"
                "with open(os.environ['DOCKER_LOG'], 'a') as log:\n"
                "    log.write(json.dumps([os.getcwd(), sys.argv[1:]]) + '\\n')\n"
            )
            docker.chmod(0o755)
            log = root / "calls.jsonl"
            env = dict(os.environ, PATH=str(root) + os.pathsep + os.environ["PATH"], DOCKER_LOG=str(log))
            script = ["bash", str(REPO / "scripts/build_dockerfiles.sh")]
            for push in (False, True):
                if log.exists():
                    log.unlink()
                subprocess.run(script + (["--push"] if push else []) +
                               ["example.org/project/rnaseq/", "v1", "umi_dup"],
                               cwd=root, env=env, check=True)
                calls = [json.loads(line) for line in log.read_text().splitlines()]
                self.assertTrue(all(cwd == str(REPO) for cwd, _ in calls))
                self.assertEqual(["build", "push", "image"] if push else ["build"],
                                 [args[0] for _, args in calls])
                self.assertEqual(["build", "--platform", "linux/amd64", "-t",
                                  "example.org/project/rnaseq/umi_dup:v1", "-f",
                                  "dockerfiles/umi_dup.Dockerfile", "."], calls[0][1])
            log.unlink()
            result = subprocess.run(script + ["example.org/rnaseq", "v1", "umi_dup", "missing"],
                                    cwd=root, env=env, capture_output=True)
            self.assertEqual(2, result.returncode)
            self.assertFalse(log.exists(), "Invalid image lists must fail before any build or push")


if __name__ == "__main__":
    unittest.main()
