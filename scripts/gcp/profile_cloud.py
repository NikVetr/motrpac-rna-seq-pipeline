"""Read-only Storage/Batch transport for one profiling snapshot, using gcloud's identity."""

from concurrent.futures import Future
from contextlib import contextmanager
import base64
import hashlib
import http.client
import io
import json
from pathlib import Path
import re
import subprocess
import tempfile
import threading
import time
from urllib.error import HTTPError, URLError
from urllib.parse import quote, unquote, urlencode, urlsplit
from urllib.request import Request, getproxies, proxy_bypass


class ProfileCloud:
    """Share reads within a capture; keep payloads on disk, never across snapshots."""

    def __init__(self):
        self._directory = tempfile.TemporaryDirectory(prefix="rnaseq-profile-reads-")
        self._lock = threading.Lock()
        self._auth_lock = threading.Lock()
        self._reads = {}
        self._idle = {}
        self._connections = set()
        self._token = None
        self._token_time = 0
        self.requests = 0
        self.cache_hits = 0

    def __enter__(self):
        return self

    def __exit__(self, *args):
        for connection in self._connections:
            connection.close()
        self._directory.cleanup()

    @contextmanager
    def open_url(self, request, timeout):
        """Reuse TLS connections across worker pools, including HTTP CONNECT proxies."""
        parsed = urlsplit(request.full_url)
        host = parsed.netloc
        with self._lock:
            idle = self._idle.setdefault(host, [])
            connection = idle.pop() if idle else None
        if connection is None:
            proxy = getproxies().get("https") if not proxy_bypass(parsed.hostname) else None
            if proxy:
                target = urlsplit(proxy)
                if target.scheme != "http":
                    raise ValueError("profiling REST reads support HTTP CONNECT proxies; use --transport gcloud for other proxies")
                connection = http.client.HTTPSConnection(target.hostname, target.port or 80, timeout=timeout)
                headers = {}
                if target.username is not None:
                    credentials = unquote(target.username) + ":" + unquote(target.password or "")
                    headers["Proxy-Authorization"] = "Basic " + base64.b64encode(credentials.encode()).decode()
                connection.set_tunnel(parsed.hostname, parsed.port or 443, headers=headers)
            else:
                connection = http.client.HTTPSConnection(parsed.hostname, parsed.port or 443, timeout=timeout)
            with self._lock:
                self._connections.add(connection)
        reusable = False
        try:
            connection.request("GET", parsed.path + ("?" + parsed.query if parsed.query else ""),
                               headers=dict(request.header_items()))
            response = connection.getresponse()
            if response.status != 200:
                data = response.read(4096)
                response.close()
                raise HTTPError(request.full_url, response.status, response.reason, response.headers, io.BytesIO(data))
            try:
                yield response
                reusable = response.isclosed()
            finally:
                response.close()
        except (OSError, http.client.HTTPException) as error:
            if isinstance(error, HTTPError):
                raise
            raise URLError(str(error)) from error
        finally:
            with self._lock:
                if reusable:
                    self._idle[host].append(connection)
                else:
                    connection.close()
                    self._connections.discard(connection)

    def token(self, expired=None):
        with self._auth_lock:
            if self._token is None or time.monotonic() - self._token_time > 3000 or (
                    expired is not None and self._token == expired):
                try:
                    self._token = subprocess.check_output(
                        ["gcloud", "auth", "print-access-token"], stderr=subprocess.PIPE).decode().strip()
                except subprocess.CalledProcessError as error:
                    raise ValueError(error.stderr.decode(errors="replace").strip()) from error
                if not self._token:
                    raise ValueError("gcloud returned an empty access token")
                self._token_time = time.monotonic()
            return self._token

    def request(self, url, params=None, media=False):
        if params:
            url += "?" + urlencode(params)
        token = self.token()
        refreshed = False
        for attempt in range(4):
            request = Request(url, headers={"Authorization": "Bearer " + token,
                                           "Accept-Encoding": "gzip" if media else "identity"})
            with self._lock:
                self.requests += 1
            try:
                with self.open_url(request, timeout=60) as response:
                    payload = response.read(64 * 1024**2 + 1) if media else response.read()
                if media and len(payload) > 64 * 1024**2:
                    raise ValueError("profiling download exceeds 64 MiB")
                return payload
            except HTTPError as error:
                code = error.code
                detail = error.read(4096).decode(errors="replace")
                error.close()
                if code == 401 and not refreshed and attempt < 3:
                    token = self.token(expired=token)
                    refreshed = True
                    continue
                if code not in (429, 500, 502, 503, 504) or attempt == 3:
                    label = "not found" if code == 404 else "request failed"
                    raise ValueError(f"{code} {label}: {url}: {detail}") from error
            except (URLError, OSError, http.client.HTTPException) as error:
                if attempt == 3:
                    raise ValueError(f"request failed: {url}: {error}") from error
            time.sleep(2 ** attempt)
        raise AssertionError("request retry loop exhausted")

    @staticmethod
    def storage_uri(uri):
        parsed = urlsplit(uri)
        if parsed.scheme != "gs" or not parsed.netloc or not parsed.path.lstrip("/") or parsed.query:
            raise ValueError(f"invalid Storage URI: {uri}")
        if parsed.fragment and not parsed.fragment.isdigit():
            raise ValueError(f"invalid object generation: {uri}")
        return parsed.netloc, parsed.path.lstrip("/"), parsed.fragment

    def fetch(self, args):
        if args[:3] == ("storage", "objects", "list"):
            if not args[3].endswith("/**"):
                raise ValueError("profiling listings require a workflow-root wildcard")
            bucket, prefix, generation = self.storage_uri(args[3][:-3])
            if generation:
                raise ValueError("a listing cannot specify a generation")
            url = "https://storage.googleapis.com/storage/v1/b/" + quote(bucket, safe="") + "/o"
            params = {"prefix": prefix + "/", "maxResults": 1000,
                      "fields": "items(name,size,generation),nextPageToken"}
            objects, tokens = [], set()
            while True:
                page = json.loads(self.request(url, params))
                objects.extend({"bucket": bucket, **obj} for obj in page.get("items", []))
                token = page.get("nextPageToken")
                if not token:
                    return json.dumps(objects).encode()
                if token in tokens:
                    raise ValueError("Storage listing repeated a page token")
                tokens.add(token)
                params["pageToken"] = token
        if args[:3] == ("storage", "objects", "describe") or args[:2] == ("storage", "cat"):
            media = args[:2] == ("storage", "cat")
            uri = args[2] if media else args[3]
            bucket, name, generation = self.storage_uri(uri)
            url = ("https://storage.googleapis.com/storage/v1/b/" + quote(bucket, safe="") +
                   "/o/" + quote(name, safe=""))
            params = {"generation": generation} if generation else {}
            if media:
                if not generation:
                    raise ValueError("profiling downloads require an explicit generation")
                params["alt"] = "media"
            else:
                params["fields"] = "size,generation"
            return self.request(url, params, media=media)
        if args[:3] == ("batch", "jobs", "describe"):
            job = args[3]
            if not re.fullmatch(r"projects/[^/]+/locations/[^/]+/jobs/[^/]+", job):
                raise ValueError(f"invalid Batch job resource name: {job}")
            return self.request("https://batch.googleapis.com/v1/" + quote(job, safe="/"))
        raise ValueError(f"unsupported profiling read: {args[:3]}")

    def __call__(self, *args):
        with self._lock:
            owner = args not in self._reads
            if owner:
                self._reads[args] = Future()
            else:
                self.cache_hits += 1
            future = self._reads[args]
        if owner:
            try:
                data = self.fetch(args)
                path = Path(self._directory.name) / hashlib.sha256(json.dumps(args).encode()).hexdigest()
                path.write_bytes(data)
                future.set_result(path)
            except BaseException as error:
                future.set_exception(error)
                raise
        return future.result().read_bytes()
