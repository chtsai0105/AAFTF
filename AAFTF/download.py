"""Download and cache AAFTF reference databases.

This module provides a one-time download command that downloads all (or a
user-selected subset) of the reference databases listed in
``AAFTF.resources`` into a persistent ``$AAFTF_DB`` folder.  Subsequent
AAFTF commands will reuse these cached files instead of re-downloading
them on first use.
"""

import os
import shutil
import sys
import urllib.request
from pathlib import Path

from AAFTF.resources import FCSADAPTOR, Contaminant_Accessions, DB_Links
from AAFTF.utility import SafeRemove, status


class _Redirect308Handler(urllib.request.HTTPRedirectHandler):
    """Extend urllib's redirect handler to also follow HTTP 308.

    Python < 3.11 does not handle 308 (Permanent Redirect).  The base
    class ``redirect_request()`` hard-codes the allowed set to
    {301,302,303,307} and raises HTTPError for anything else, so we must
    override both that method and add the http_error_308 dispatcher.
    """

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if code == 308:
            code = 307  # method-preserving permanent redirect
        return super().redirect_request(req, fp, code, msg, headers, newurl)

    def http_error_308(self, req, fp, code, msg, headers):
        return self.http_error_302(req, fp, code, msg, headers)


_opener = urllib.request.build_opener(_Redirect308Handler())


def _resolve_db_dir():
    """Resolve the AAFTF_DB directory from $AAFTF_DB, exiting with an error if it isn't set."""
    db_dir = os.environ.get("AAFTF_DB")
    if not db_dir:
        status("ERROR: No database directory specified.\n" "  Set the AAFTF_DB environment variable.\n" "  Example:\n" "    export AAFTF_DB=/path/to/aaftf_db\n" "    AAFTF download")
        sys.exit(1)
    return str(Path(db_dir).resolve())


def _human_size(num_bytes):
    """Format a byte count as a human-readable string (e.g. ``1.2 GB``)."""
    size = float(num_bytes)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024 or unit == "TB":
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"  # pragma: no cover - unreachable, but keeps a return for every path


def _expected_db_entries():
    """Return [(filename, url), ...] for every database ``download`` would fetch."""
    entries = []
    for urls in Contaminant_Accessions.values():
        for url in urls:
            entries.append((Path(url).name, url))
    for key, urls in DB_Links.items():
        for url_or_meta in urls:
            if isinstance(url_or_meta, dict):
                entries.append((url_or_meta["filename"], url_or_meta["url"]))
            else:
                entries.append((Path(url_or_meta).name, url_or_meta))
    entries.append(("run_fcsadaptor.sh", FCSADAPTOR["EXEURL"] % FCSADAPTOR["VERSION"]))
    entries.append((FCSADAPTOR["SIFLOCAL"] % FCSADAPTOR["VERSION"], str(Path(FCSADAPTOR["SIFURL"], FCSADAPTOR["VERSION"], FCSADAPTOR["SIF"]))))
    return entries


def _remote_size(url):
    """Return the remote file size in bytes via an HTTP HEAD request, or ``None`` if unavailable."""
    try:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "AAFTF/1.0"})
        with _opener.open(req, timeout=10) as response:
            length = response.headers.get("Content-Length")
            return int(length) if length is not None else None
    except Exception:
        return None


def _list_db(db_dir):
    """List every database ``download`` would fetch: local size if already downloaded, else the remote size."""
    db_path = Path(db_dir)

    status(f"Database files in {db_dir}:")
    downloaded_total = 0
    pending_total = 0
    pending_unknown = False
    for filename, url in _expected_db_entries():
        local_path = db_path / filename
        if local_path.is_file():
            size = local_path.stat().st_size
            downloaded_total += size
            print(f"  {filename:<40} {_human_size(size):>10}  (downloaded)")
        else:
            size = _remote_size(url)
            if size is None:
                pending_unknown = True
                print(f"  {filename:<40} {'unknown':>10}  (not downloaded)")
            else:
                pending_total += size
                print(f"  {filename:<40} {_human_size(size):>10}  (not downloaded)")

    # Any other files present that aren't part of the expected set (e.g. leftover reports)
    expected_names = {filename for filename, _ in _expected_db_entries()}
    if db_path.is_dir():
        extras = sorted(f for f in db_path.rglob("*") if f.is_file() and f.name not in expected_names)
        for f in extras:
            size = f.stat().st_size
            downloaded_total += size
            print(f"  {str(f.relative_to(db_path)):<40} {_human_size(size):>10}  (other)")

    print("-" * 68)
    print(f"  {'Downloaded':<40} {_human_size(downloaded_total):>10}")
    pending_label = _human_size(pending_total) + ("+" if pending_unknown else "")
    print(f"  {'Not yet downloaded (estimated)':<40} {pending_label:>10}")


def run(force=False, skip_core=False, skip_sourmash=False, skip_fcs=False, sourdb_type="all", list_db=False, **kwargs):
    """Execute the ``download`` subcommand.

    Downloads reference databases to the ``AAFTF_DB`` directory so that
    later AAFTF commands do not need to fetch them on-the-fly. If
    ``list_db`` is set, lists the existing database files (and their sizes)
    instead of downloading anything.
    """
    db_dir = _resolve_db_dir()

    if list_db:
        _list_db(db_dir)
        return

    Path(db_dir).mkdir(parents=True, exist_ok=True)
    status(f"AAFTF database directory: {db_dir}")

    errors = []

    # Core databases (always downloaded unless --skip-core)
    if not skip_core:
        try:
            _download_contaminants(db_dir, force=force)
        except Exception as e:
            errors.append(f"Contaminant accessions: {e}")
        try:
            _download_db_links(db_dir, force=force)
        except Exception as e:
            errors.append(f"DB_Links: {e}")

    # Sourmash taxonomy databases (optional)
    if not skip_sourmash:
        try:
            _download_sourmash(db_dir, sourdb_type=sourdb_type, force=force)
        except Exception as e:
            errors.append(f"Sourmash DB: {e}")

    # NCBI FCS-adaptor resources (optional)
    if not skip_fcs:
        try:
            _download_fcs(db_dir, force=force)
        except Exception as e:
            errors.append(f"FCS resources: {e}")

    if errors:
        status("\nSome downloads failed:")
        for err in errors:
            status(f"  - {err}")
        status("\nYou can re-run 'AAFTF download' later; already-downloaded " "files will be skipped unless you pass --force.")
        sys.exit(1)

    status("Setup complete. Future AAFTF runs will use cached files from " f"{db_dir}")


def _download(url, dest, force=False):
    """Download ``url`` to ``dest`` if it does not already exist.

    Writes to a temporary file first and renames atomically on success so
    that an interrupted download never leaves a partial file that would be
    mistaken for a complete one on the next run.

    Args:
        url: Remote URL to download.
        dest: Local file path to write.
        force: If True, re-download even if ``dest`` exists.

    Returns:
        The absolute path to the downloaded file.
    """
    if Path(dest).exists() and not force:
        status(f"  Already present: {dest}")
        return dest

    status(f"  Downloading {Path(dest).name} ...")
    Path(dest).parent.mkdir(parents=True, exist_ok=True)

    tmp = dest + ".tmp"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "AAFTF/1.0"})
        with _opener.open(req, timeout=300) as response:
            final_url = response.geturl()
            if final_url != url:
                status(f"  Redirected to {final_url}")
            with open(tmp, "wb") as outfh:
                shutil.copyfileobj(response, outfh)
        os.rename(tmp, dest)
    except Exception as e:
        status(f"  ERROR downloading {url}: {e}")
        if Path(tmp).exists():
            SafeRemove(tmp)
        raise

    status(f"  Saved {dest}")
    return dest


def _download_db_links(db_dir, keys=None, force=False):
    """Download files referenced in ``DB_Links``.

    Args:
        db_dir: Target directory.
        keys: Iterable of ``DB_Links`` keys to download, or ``None``
            for the subset used by *filter* / *vecscreen*.
        force: If True, overwrite existing files.
    """
    if keys is None:
        # Default set needed by filter / vecscreen
        keys = ("UniVec", "CONTAM_EUKS", "CONTAM_PROKS", "MITO")

    for key in keys:
        if key not in DB_Links:
            status(f"  WARNING: unknown DB_Links key '{key}', skipping")
            continue
        status(f"Downloading {key} ...")
        for url_or_meta in DB_Links[key]:
            if isinstance(url_or_meta, dict):
                url = url_or_meta["url"]
                filename = url_or_meta["filename"]
            else:
                url = url_or_meta
                filename = Path(url).name
            dest = str(Path(db_dir, filename))
            _download(url, dest, force=force)


def _download_contaminants(db_dir, force=False):
    """Download contaminant accessions (e.g. PhiX).

    Args:
        db_dir: Target directory.
        force: If True, overwrite existing files.
    """
    status("Downloading contaminant accessions ...")
    for name, urls in Contaminant_Accessions.items():
        status(f"  {name} ...")
        for url in urls:
            filename = Path(url).name
            dest = str(Path(db_dir, filename))
            _download(url, dest, force=force)


def _download_sourmash(db_dir, sourdb_type="gbk", force=False):
    """Download sourmash LCA taxonomy databases.

    Args:
        db_dir: Target directory.
        sourdb_type: One of ``gbk``, ``gtdb``, ``gtdbrep``, or ``all``.
        force: If True, overwrite existing files.
    """
    type_map = {
        "gbk": "sourmash_gbk",
        "gtdb": "sourmash_gtdb",
        "gtdbrep": "sourmash_gtdbrep",
    }

    if sourdb_type == "all":
        indices = list(type_map.values())
    else:
        if sourdb_type not in type_map:
            status(f"  ERROR: unknown sourdb_type '{sourdb_type}'. " f"Choose from {list(type_map.keys()) + ['all']}")
            return
        indices = [type_map[sourdb_type]]

    for idx in indices:
        status(f"Downloading sourmash database ({idx}) ...")
        for entry in DB_Links[idx]:
            dest = str(Path(db_dir, entry["filename"]))
            _download(entry["url"], dest, force=force)


def _download_fcs(db_dir, force=False):
    """Download NCBI FCS-adaptor script and SIF image.

    Args:
        db_dir: Target directory.
        force: If True, overwrite existing files.
    """
    status("Downloading NCBI FCS-adaptor resources ...")

    # Wrapper script
    script_url = FCSADAPTOR["EXEURL"] % FCSADAPTOR["VERSION"]
    script_dest = str(Path(db_dir, "run_fcsadaptor.sh"))
    if Path(script_dest).exists() and not force:
        status(f"  Already present: {script_dest}")
    else:
        _download(script_url, script_dest, force=force)
        os.chmod(script_dest, 0o555)

    # Singularity image
    image_name = FCSADAPTOR["SIFLOCAL"] % FCSADAPTOR["VERSION"]
    image_dest = str(Path(db_dir, image_name))
    if Path(image_dest).exists() and not force:
        status(f"  Already present: {image_dest}")
    else:
        image_url = str(Path(FCSADAPTOR["SIFURL"], FCSADAPTOR["VERSION"], FCSADAPTOR["SIF"]))
        _download(image_url, image_dest, force=force)
