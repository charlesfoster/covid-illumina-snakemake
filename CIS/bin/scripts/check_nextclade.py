#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import shlex

def check_nextclade(thisdir):
    snakefile = os.path.join(thisdir, "bin", "rules", "routine.smk")
    print(
        "Checking whether nextclade needs to be updated..."
    )
    update_file = os.path.join(thisdir, "bin", "scripts", "get_latest_tag.sh")
    cmd = f'grep "docker://nextstrain" {snakefile} | uniq | tr -d " " | sed "s|docker://nextstrain/nextclade:||" | tr -d "\n"'
    current_version = subprocess.check_output(cmd, shell=True).decode('utf-8').replace('"','')
    cmd = f"bash {update_file} --source dockerhub --repo nextstrain --image nextclade --no-alphabet"
    newest_version = subprocess.check_output(cmd, shell=True).decode('utf-8')
    if current_version != newest_version:
        cmd = f'sed -i "s|docker://nextstrain/nextclade:{current_version}|docker://nextstrain/nextclade:{newest_version}|g" {snakefile}'
        update_proc = subprocess.Popen(
            shlex.split(cmd),
            shell=False,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            )
        print(
            f"Updated nextclade version dependency to {newest_version}"
        )
    else:
        print(
            f"Nextclade version dependency is already at the newest version: {current_version}"
        )
