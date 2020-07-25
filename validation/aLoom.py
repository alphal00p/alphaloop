'''
Avoid OOM state due to aL generation during compilation
'''

from datetime import datetime
from time import sleep

import psutil
from psutil import NoSuchProcess, Process

VMEM_LIMIT_PERC = 95


def check_mem(processes):
    vmem = psutil.virtual_memory()
    if vmem.percent > VMEM_LIMIT_PERC:
        print(datetime.now())
        for process in processes:
            cmdline = process.cmdline()
            command_short = "{} [...] {}".format(
                ' '.join(cmdline[:2]), cmdline[-1])
            print("\t", command_short)
            process.kill()
        return True
    return False


def get_aL_pids():
    aL_pids = set()
    pids = psutil.pids()
    for pid in pids:
        try:
            p = Process(pid)
        except NoSuchProcess:
            pass
        if '--mode=alphaloop' in p.cmdline():
            aL_pids.add(pid)
    return aL_pids


while True:
    aL_pids = get_aL_pids()
    aL_processes = []
    for pid in aL_pids:
        try:
            aL_processes.append(Process(pid))
        except NoSuchProcess:
            pass
    base_vmem = sum(p.memory_info().vms for p in aL_processes)
    for p in aL_processes:
        aL_vmem = base_vmem
        childrens = set(filter(lambda cp: not cp.pid in aL_pids,
                               p.children(recursive=True)))
        for c in childrens:
            aL_vmem += c.memory_info().vms
        if check_mem(childrens):
            sleep(1)
            print("   >> Killed aL childrens when reached {:.2f}G, freed {:.2f}%"
                  .format(aL_vmem/1024**3,
                          VMEM_LIMIT_PERC - psutil.virtual_memory().percent))
            break
    sleep(1)
