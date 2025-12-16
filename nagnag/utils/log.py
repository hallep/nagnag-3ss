import sys
from datetime import datetime

def TIME():
    return f"[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}]"

def log_script(name:str):
    print(f"{TIME()} ---------- {name} ----------", file=sys.stderr)

def log_fn(desc:str, sub:bool=False):
    print(f"{TIME()}{"\t" if sub else " "}{desc}", file=sys.stderr)
