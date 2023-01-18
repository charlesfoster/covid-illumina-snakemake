#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import shlex

def check_freyja(thisdir,workflow):
    snakefile = os.path.join(thisdir, "bin", "rules", f"{workflow}.smk")
    print(
        "Checking whether freyja needs to be updated..."
    )
    print(f"Working with: {snakefile}")
    update_file = os.path.join(thisdir, "bin", "scripts", "get_latest_tag.sh")
    cmd = f'grep "docker://staphb" {snakefile} | uniq | tr -d " " | sed "s|docker://staphb/freyja:||" | tr -d "\n"'
    current_version = subprocess.check_output(cmd, shell=True).decode('utf-8').replace('"','')
    cmd = f"bash {update_file} --source dockerhub --repo staphb --image freyja"
    newest_version = subprocess.check_output(cmd, shell=True).decode('utf-8')
    if current_version != newest_version:
        cmd = f'sed -i "s|docker://staphb/freyja:{current_version}|docker://staphb/freyja:{newest_version}|g" {snakefile}'
        update_proc = subprocess.Popen(
            shlex.split(cmd),
            shell=False,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            )
        print(
            f"Updated freyja version dependency to {newest_version}"
        )
    else:
        print(
            f"freyja version dependency is already at the newest version: {current_version}"
        )
