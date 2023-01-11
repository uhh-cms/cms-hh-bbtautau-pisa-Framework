#!/usr/bin/env python3
# coding: utf-8
# flake8: noqa

import os
import subprocess

import yaml
from collections import OrderedDict


thisdir = os.path.dirname(os.path.abspath(__file__))


def main(config_file):
    # read the config
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    # set global settings first
    remote_bases = {}
    for name, data in config.items():
        if name == "GLOBAL":
            if "remoteBases" in data:
                remote_bases = dict(data["remoteBases"])
            break

    timestamps = OrderedDict()

    # loop through datasets
    for i, (name, data) in enumerate(config.items()):
        if name == "GLOBAL":
            continue

        # validate
        if "remoteBase" not in data:
            print(f"skip {name}, no remoteBase entry")
            continue

        # read settings
        dataset_key = data["miniAOD"]
        remote_base = data["remoteBase"]
        remote_base = remote_bases.get(remote_base, remote_base)

        # read the directory
        remote_dir = os.path.join(remote_base, dataset_key.split("/")[1], "crab_" + name)
        cmd = f"gfal-ls {remote_dir}"
        print(f"reading {name}")
        p, out = run_command(cmd)
        timestamps[name] = " or ".join(f"\"{t}\"" for t in out.split("\n")) if p.returncode == 0 else "null"

    # print timestamps
    print("\n" + "-" * 80 + "\n")
    for name, timestamp in timestamps.items():
        print(f"{name}:\n    timestamp: {timestamp}\n")


def run_command(cmd, timeout=None, attempts=1, _attempt=1, **kwargs):
    kwargs["preexec_fn"] = os.setsid
    p = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, **kwargs)

    try:
        out, _ = p.communicate(timeout=timeout)
        return p, out.decode("utf-8").strip()

    except KeyboardInterrupt:
        p.terminate()
        try:
            pgid = os.getpgid(p.pid)
            if p.poll() is None:
                os.killpg(pgid, signal.SIGTERM)
        except:
            pass
        raise

    except Exception:
        if attempts > _attempt:
            return run_command(cmd, timeout=timeout, attempts=attempts, _attempt=_attempt + 1)

        print(f"command failed after {attempts} attempt(s):\n{cmd}")
        raise


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="list timestamps of crab submissions")
    parser.add_argument(
        "--config",
        "-c",
        default="samples_2017_uhh.yaml",
        help="the sample config to read",
    )

    args = parser.parse_args()

    main(
        os.path.join(thisdir, args.config),
    )
