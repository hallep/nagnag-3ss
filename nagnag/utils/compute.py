''' Parallel computing functions '''

from collections.abc import Callable, Iterable
import multiprocessing as mp

def single_map(func:Callable, inputs:Iterable, n_procs:int=24):
    return mp.Pool(processes=n_procs).map(func, inputs)

def multi_map(func:Callable, inputs:Iterable[Iterable], n_procs:int=24):
    return mp.Pool(processes=n_procs).starmap(func, inputs)
