import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sys
import unittest
from unittest.mock import Mock, patch
from urllib.error import HTTPError
from urllib.parse import parse_qs, urlsplit

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/gcp"))
from profile_cloud import ProfileCloud


class ProfileCloudTests(unittest.TestCase):
    def test_connections_reused_across_worker_pools_and_closed_after_read_failure(self):
        response = Mock(status=200)
        response.read.return_value = b"{}"
        response.isclosed.return_value = True
        with patch("profile_cloud.subprocess.check_output", return_value=b"token"), \
                patch("profile_cloud.getproxies", return_value={"https": "http://proxy:8080"}), \
                patch("profile_cloud.proxy_bypass", return_value=False), \
                patch("profile_cloud.http.client.HTTPSConnection") as factory:
            connection = factory.return_value
            connection.getresponse.return_value = response
            with ProfileCloud() as reader:
                for _ in range(3):
                    with ThreadPoolExecutor(max_workers=1) as pool:
                        self.assertEqual(b"{}", pool.submit(reader.request, "https://storage.googleapis.com/x").result())
                factory.assert_called_once_with("proxy", 8080, timeout=60)
                connection.set_tunnel.assert_called_once_with("storage.googleapis.com", 443, headers={})
                connection.close.assert_not_called()
                response.read.side_effect = [OSError("interrupted"), b"{}"]
                with patch("profile_cloud.time.sleep"):
                    self.assertEqual(b"{}", reader.request("https://storage.googleapis.com/x"))
                self.assertEqual(2, factory.call_count)
                connection.close.assert_called_once()
            self.assertEqual(2, connection.close.call_count)

    def test_paginated_listing_and_concurrent_reuse_preserve_all_objects(self):
        urls = []
        def open_url(request, timeout):
            urls.append(request.full_url)
            query = parse_qs(urlsplit(request.full_url).query)
            self.assertEqual(["x/rnaseq_pipeline/wf/"], query["prefix"])
            if "pageToken" not in query:
                result = {"items": [{"name": "a", "size": "7", "generation": "8"}], "nextPageToken": "next"}
            else:
                self.assertEqual(["next"], query["pageToken"])
                result = {"items": [{"name": "b", "size": "9", "generation": "10"}]}
            return io.BytesIO(json.dumps(result).encode())
        with patch("profile_cloud.subprocess.check_output", return_value=b"token") as auth, \
                patch("profile_cloud.ProfileCloud.open_url", side_effect=open_url), ProfileCloud() as reader:
            args = ("storage", "objects", "list", "gs://bucket/x/rnaseq_pipeline/wf/**", "--format=json")
            with ThreadPoolExecutor(max_workers=4) as pool:
                results = list(pool.map(lambda _: reader(*args), range(4)))
            self.assertEqual(2, len(urls))
            self.assertEqual(1, auth.call_count)
            self.assertEqual(3, reader.cache_hits)
            self.assertEqual([{"bucket": "bucket", "name": "a", "size": "7", "generation": "8"},
                              {"bucket": "bucket", "name": "b", "size": "9", "generation": "10"}],
                             json.loads(results[0]))
            self.assertTrue(all(r == results[0] for r in results))

    def test_downloads_are_generation_pinned_and_cache_is_snapshot_local(self):
        def open_url(request, timeout):
            parsed = urlsplit(request.full_url)
            self.assertTrue(parsed.path.endswith("/o/folder%2Fa%20b.txt"))
            query = parse_qs(parsed.query)
            self.assertEqual(["media"], query["alt"])
            return io.BytesIO(query["generation"][0].encode())
        with patch("profile_cloud.subprocess.check_output", return_value=b"token"), \
                patch("profile_cloud.ProfileCloud.open_url", side_effect=open_url) as request:
            with ProfileCloud() as reader:
                for generation in ("123", "124", "123"):
                    self.assertEqual(generation.encode(), reader("storage", "cat", "gs://b/folder/a b.txt#" + generation))
                self.assertEqual(2, request.call_count)
                with self.assertRaisesRegex(ValueError, "explicit generation"):
                    reader("storage", "cat", "gs://b/folder/a b.txt")
            with ProfileCloud() as reader:
                self.assertEqual(b"123", reader("storage", "cat", "gs://b/folder/a b.txt#123"))
            self.assertEqual(3, request.call_count)

    def test_transient_retry_and_expired_token_refresh_do_not_mask_permission_errors(self):
        def failure(code):
            return HTTPError("https://example.invalid", code, "failure", {}, io.BytesIO(b"detail"))
        for code, retry in ((429, True), (503, True), (401, True), (403, False), (404, False)):
            with self.subTest(code=code), patch("profile_cloud.subprocess.check_output", return_value=b"token") as auth, \
                    patch("profile_cloud.ProfileCloud.open_url", side_effect=[failure(code), io.BytesIO(b"{}")] ) as request, \
                    patch("profile_cloud.time.sleep"), ProfileCloud() as reader:
                args = ("batch", "jobs", "describe", "projects/p/locations/r/jobs/j", "--format=json")
                if retry:
                    self.assertEqual(b"{}", reader(*args))
                    self.assertEqual(2, request.call_count)
                    self.assertEqual(2 if code == 401 else 1, auth.call_count)
                else:
                    with self.assertRaisesRegex(ValueError, f"{code} " + ("not found" if code == 404 else "request failed")):
                        reader(*args)
                    self.assertEqual(1, request.call_count)

    def test_repeated_page_token_and_invalid_batch_path_fail_loudly(self):
        page = json.dumps({"nextPageToken": "same"}).encode()
        with patch("profile_cloud.subprocess.check_output", return_value=b"token"), \
                patch("profile_cloud.ProfileCloud.open_url", side_effect=lambda *args, **kwargs: io.BytesIO(page)), ProfileCloud() as reader:
            with self.assertRaisesRegex(ValueError, "repeated a page token"):
                reader("storage", "objects", "list", "gs://b/x/rnaseq_pipeline/id/**")
            with self.assertRaisesRegex(ValueError, "invalid Batch"):
                reader("batch", "jobs", "describe", "https://untrusted.invalid/job")


if __name__ == "__main__":
    unittest.main()
