#!/usr/bin/env python3
"""
cluster_sync.py — upload/download files to/from an HPC cluster via rsync over ssh.

UPLOAD
------
Copies a local folder to a *new*, timestamped folder on the cluster:

    ./cluster_sync.py upload user@cluster:/path/where/to/copy

This creates `/path/where/to/copy/jobcraft_folder_YYYYMMDD_HHMMSS/` on the
remote host and copies everything from the local folder (default: current
directory) into it. The resulting remote path is remembered locally in a
hidden state file (`.jobcraft_sync.json`) so you don't have to type it again.

DOWNLOAD
--------
Syncs the local folder with the remote folder recorded during the last
upload — no path needed:

    ./cluster_sync.py download

This uses rsync, so only new/changed files are transferred and nothing is
duplicated. Re-running it after a job produces more output just pulls down
the new files.

Requires: `rsync` and `ssh` on PATH, and key-based (passwordless) SSH access
to the cluster is strongly recommended.

Examples
--------
    # first time: push a job folder to the cluster
    cd my_job
    ../cluster_sync.py upload hpc-login:/scratch/vondrakmar/jobs

    # later: pull results back into the same local folder
    cd my_job
    ../cluster_sync.py download

    # pull results back and remove local files no longer present remotely
    ../cluster_sync.py download --delete

    # check what's on record for this folder
    ../cluster_sync.py status
"""

import argparse
import datetime
import json
import os
import shlex
import subprocess
import sys

STATE_FILENAME = ".jobcraft_sync.json"


def state_path(local_dir):
    return os.path.join(local_dir, STATE_FILENAME)


def load_state(local_dir):
    p = state_path(local_dir)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def save_state(local_dir, data):
    with open(state_path(local_dir), "w") as f:
        json.dump(data, f, indent=2)


def run(cmd, dry_run=False):
    print("+ " + " ".join(shlex.quote(c) for c in cmd))
    if dry_run:
        return 0
    return subprocess.run(cmd).returncode


def parse_remote(remote):
    if ":" not in remote:
        sys.exit("Remote must be in the form host:/path or user@host:/path")
    host, path = remote.split(":", 1)
    if not path.startswith("/"):
        sys.exit("Remote path must be absolute (start with /)")
    return host, path


def rsync_base_args(exclude, ssh_opts):
    args = ["-avz", "--progress"]
    for pattern in exclude:
        args += ["--exclude", pattern]
    if ssh_opts:
        args += ["-e", f"ssh {ssh_opts}"]
    return args


def cmd_upload(args):
    local_dir = os.path.abspath(args.local)
    if not os.path.isdir(local_dir):
        sys.exit(f"Local path does not exist: {local_dir}")

    host, base_path = parse_remote(args.remote)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    folder_name = f"jobcraft_folder_{timestamp}"
    remote_full_path = base_path.rstrip("/") + "/" + folder_name

    ssh_cmd = ["ssh"] + shlex.split(args.ssh_opts) + [
        host,
        f"mkdir -p {shlex.quote(remote_full_path)}",
    ]
    rc = run(ssh_cmd, dry_run=args.dry_run)
    if rc != 0:
        sys.exit("Failed to create remote directory")

    rsync_args = ["rsync"] + rsync_base_args(args.exclude, args.ssh_opts) + [
        local_dir + "/",
        f"{host}:{remote_full_path}/",
    ]
    rc = run(rsync_args, dry_run=args.dry_run)
    if rc != 0:
        sys.exit(f"rsync failed with exit code {rc}")

    if not args.dry_run:
        save_state(local_dir, {
            "remote_host": host,
            "remote_path": remote_full_path,
            "local_path": local_dir,
            "ssh_opts": args.ssh_opts,
            "last_upload": datetime.datetime.now().isoformat(timespec="seconds"),
        })
        print(f"\nUploaded to {host}:{remote_full_path}")
        print(f"Saved sync target in {state_path(local_dir)}")


def cmd_download(args):
    local_dir = os.path.abspath(args.local)
    state = load_state(local_dir)

    if state is None and (args.remote_host is None or args.remote_path is None):
        sys.exit(
            f"No sync record found in {local_dir} ({STATE_FILENAME} missing).\n"
            "Run an upload first, or pass both --remote-host and --remote-path explicitly."
        )

    host = args.remote_host or state["remote_host"]
    remote_path = args.remote_path or state["remote_path"]
    ssh_opts = args.ssh_opts if args.ssh_opts is not None else (state or {}).get("ssh_opts", "")

    os.makedirs(local_dir, exist_ok=True)

    rsync_args = ["rsync"] + rsync_base_args(args.exclude, ssh_opts)
    if args.delete:
        rsync_args.append("--delete")
    rsync_args += [f"{host}:{remote_path}/", local_dir + "/"]

    rc = run(rsync_args, dry_run=args.dry_run)
    if rc != 0:
        sys.exit(f"rsync failed with exit code {rc}")

    if not args.dry_run:
        new_state = state.copy() if state else {}
        new_state.update({
            "remote_host": host,
            "remote_path": remote_path,
            "local_path": local_dir,
            "ssh_opts": ssh_opts,
            "last_download": datetime.datetime.now().isoformat(timespec="seconds"),
        })
        save_state(local_dir, new_state)
        print(f"\nSynced from {host}:{remote_path}")


def cmd_status(args):
    local_dir = os.path.abspath(args.local)
    state = load_state(local_dir)
    if state is None:
        print(f"No sync record in {local_dir}.")
        return
    print(json.dumps(state, indent=2))


def main():
    parser = argparse.ArgumentParser(
        description="Upload/download files to/from an HPC cluster via rsync over ssh."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_up = sub.add_parser(
        "upload",
        help="Create a new timestamped folder on the cluster and copy local files there.",
    )
    p_up.add_argument(
        "remote",
        help="cluster:/path/where/to/copy (e.g. user@login.hpc.edu:/scratch/user/jobs)",
    )
    p_up.add_argument("--local", default=".", help="Local folder to upload (default: current directory)")
    p_up.add_argument("--exclude", action="append", default=[], help="rsync --exclude pattern (repeatable)")
    p_up.add_argument("--ssh-opts", default="", help="Extra options passed to ssh, e.g. '-p 2222 -i ~/.ssh/id_rsa'")
    p_up.add_argument("--dry-run", action="store_true", help="Show what would happen without doing it")
    p_up.set_defaults(func=cmd_upload)
    p_down = sub.add_parser(
        "download",
        help="Sync local folder with the remote folder recorded during the last upload.",
    )
    p_down.add_argument("--local", default=".", help="Local folder to sync into (default: current directory)")
    p_down.add_argument("--remote-host", default=None, help="Override remote host (default: from saved state)")
    p_down.add_argument("--remote-path", default=None, help="Override remote path (default: from saved state)")
    p_down.add_argument("--ssh-opts", default=None, help="Extra options passed to ssh (default: from saved state)")
    p_down.add_argument("--exclude", action="append", default=[], help="rsync --exclude pattern (repeatable)")
    p_down.add_argument("--delete", action="store_true", help="Delete local files that no longer exist on remote")
    p_down.add_argument("--dry-run", action="store_true", help="Show what would happen without doing it")
    p_down.set_defaults(func=cmd_download)

    p_status = sub.add_parser("status", help="Show the saved sync record for a local folder.")
    p_status.add_argument("--local", default=".", help="Local folder to inspect (default: current directory)")
    p_status.set_defaults(func=cmd_status)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
