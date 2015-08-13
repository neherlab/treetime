import numpy as np

def opt(n, a, b, c):
    return a+b+c

def a(callable, args):
    n = 0
    print (opt(n, *args))


if __name__=='__main__':

    a(0, (1,2))
