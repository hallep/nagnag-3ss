import argparse, subprocess
import pandas as pd
from pathlib import Path
import maskpass
from sqlalchemy import URL, create_engine

parser = argparse.ArgumentParser(prog="export_hgmd.py", add_help=False,
                                 description="Export the 'splice' and 'hg38_coords' tables from HGMD Pro through MySQL")
parser.add_argument("--help", action="help", help="Show this help message and exit")
parser.add_argument("-u", "--user", type=str, default="root", help="mysql user [root]")
parser.add_argument("-h", "--host", type=str, default="localhost", help="mysql host [localhost]")
args = parser.parse_args()

# create directory
pro_dir = f"{Path(__file__).parent}/src/HGMD_pro"
subprocess.run(["mkdir", "-p", pro_dir])

# get password
password = maskpass.askpass(mask="")

# connect to database
engine = create_engine(URL.create(
    "mysql+mysqlconnector",
    username=args.user,
    password=password,
    host=args.host,
    database="hgmd_pro"
))

# export tables
for table in ["splice", "hg38_coords"]:
    df = pd.read_sql(sql=f"SELECT * FROM {table}", con=engine)
    df.to_csv(f"{pro_dir}/{table}.txt", sep="\t", index=False)
