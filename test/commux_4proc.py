#!/usr/bin/env python3
"""Run a real four-process commux/UCX communication check."""

import os
import socket
import sys

import commux
import torch
import torch.distributed as dist
import torch.multiprocessing as mp


def _free_local_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


def _pin_to_core(rank):
    try:
        cpus = sorted(os.sched_getaffinity(0))
        if len(cpus) >= 4:
            os.sched_setaffinity(0, {cpus[rank]})
            return cpus[rank]
    except (AttributeError, OSError):
        pass
    return None


def _worker(rank, world_size, port):
    core = _pin_to_core(rank)
    commux.register()
    dist.init_process_group(
        backend="ucx",
        init_method=f"tcp://127.0.0.1:{port}",
        rank=rank,
        world_size=world_size,
    )

    value = torch.tensor([rank + 1], dtype=torch.int64)
    dist.all_reduce(value, op=dist.ReduceOp.SUM)
    if value.item() != 10:
        raise RuntimeError(f"rank {rank}: all_reduce returned {value.item()}")

    bcast = torch.tensor([42 if rank == 2 else 0], dtype=torch.int64)
    dist.broadcast(bcast, src=2)
    if bcast.item() != 42:
        raise RuntimeError(f"rank {rank}: broadcast returned {bcast.item()}")

    if rank == 0:
        dist.send(torch.tensor([123], dtype=torch.int64), dst=1, tag=77)
    elif rank == 1:
        recv = torch.zeros(1, dtype=torch.int64)
        dist.recv(recv, src=0, tag=77)
        if recv.item() != 123:
            raise RuntimeError(f"rank {rank}: recv returned {recv.item()}")

    dist.barrier()
    dist.destroy_process_group()
    print(f"rank {rank}/4 ok" + (f" on cpu {core}" if core is not None else ""))


def main():
    world_size = 4
    if hasattr(os, "sched_getaffinity") and len(os.sched_getaffinity(0)) < world_size:
        print("SKIP: fewer than 4 CPU cores are available", file=sys.stderr)
        return 77
    mp.spawn(_worker, args=(world_size, _free_local_port()), nprocs=world_size)
    print("PASS: commux UCX 4-process test")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
