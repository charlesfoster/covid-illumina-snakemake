#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import shlex
from cmp_version import cmp_version

def check_freyja(thisdir,workflow):
    snakefile = os.path.join(thisdir, "bin", "rules", f"{workflow}.smk")
    update_file = os.path.join(thisdir, "bin", "scripts", "get_latest_tag.sh")
    print(
        "Checking whether freyja needs to be updated..."
    )
    print(f"Working with: {snakefile}")
    print(f"Checking StaPH-B...")
    cmd = f"bash {update_file} --source dockerhub --repo staphb --image freyja"
    newest_staphb = subprocess.check_output(cmd, shell=True).decode('utf-8')
    print(f"Checking biocontainers...")
    cmd = f'bash {update_file} --source quay --repo biocontainers --image freyja'
    newest_biocontainers = subprocess.check_output(cmd, shell=True).decode('utf-8')
    comparison = cmp_version(newest_staphb,newest_biocontainers)
    if comparison == 0:
        print(f"StaPH-B and biocontainers are equivalent --> using StaPH-B")
        newest_version = newest_staphb
        choice = 'staphb'
    elif comparison == -1:
        print(f"Using biocontainers container")
        newest_version = newest_biocontainers
        choice = 'bioc'
    elif comparison == 1:
        print(f"Using StaPH-B container")
        newest_version = newest_staphb
        choice = 'staphb'
    cmd = f'grep "docker://.*freyja" {snakefile} | uniq | tr -d " " | sed "s|docker://.*:||" | tr -d "\n"'
    current_version = subprocess.check_output(cmd, shell=True).decode('utf-8').replace('"','')
    comparison = cmp_version(current_version,newest_version)
    if comparison == 0:
        print(
            f"freyja version dependency is already at the newest version: {current_version}"
        )
    else:
        if choice == 'staphb':
            cmd = f'sed -i "s|docker://.*/freyja:{current_version}|docker://staphb/freyja:{newest_version}|g" {snakefile}'
        elif choice == 'bioc':
            cmd = f'sed -i "s|docker://.*/freyja:{current_version}|docker://quay.io/biocontainers/freyja:{newest_version}|g" {snakefile}'
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
