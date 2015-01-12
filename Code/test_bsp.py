from Scientific.BSP import ParSequence, ParFunction, ParRootFunction
from sys import argv
from random import randint

def sumWat(ws): return len(p[w for w in ws if sum(w)==k])
global_sumWat = ParFunction(sumWat)
def asynchronous(result): print( float(result) / len(wat) )
global_asynchronous = ParRootFunction(output)
wat = [[randint(1,6) for _ in range(n)] for _ in range(m)]
global_wat = ParSequence(wat)

# the parallel calc:
err_argument = global_sumWat(global_wat)
# collecting the results in processor 0:
err_argument_gs = results.reduce(lambda x, y: x + y, 0)

global_output(err_argument_gs)
