import os
import dask_mpi as dm
import dask.distributed as dd

# Initialise Dask cluster and store worker files in current work directory
dm.initialize(local_directory=os.getcwd())

# Define two simple test functions
def inc(x):
 return x + 1

def add(x, y):
 return x + y

client = dd.Client()

# Submit chain of computations using futures
a = client.submit(inc, 1)
b = client.submit(inc, 2)
c = client.submit(add, a, b)

# Expect the same answer
print("Dask result:", c.result())
print("Local result:", add(inc(1), inc(2)))
