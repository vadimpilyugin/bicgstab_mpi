import os
import re
from os import path as path
import remap as rm

class Experiment:
  def __init__(self, n, op, time, ntr, nproc, supercomp, speedup=None, test_no=None):
    self.n = n
    self.op = op
    self.time = time
    self.ntr = ntr
    self.nproc = nproc
    self.supercomp = supercomp
    self.speedup = speedup
    self.test_no = test_no

  def __repr__(self):
    return f"FROM: {self.supercomp}: N={self.n} {self.op}\ttime={self.time}s NTR={self.ntr} nproc={self.nproc}"

def get_test_times(basename):
  md = re.match(rf'test_\d+_(\d+)', basename)
  return int(md[1])

def read_tests(test_no, supercomp="bluegene"):
  tests_all = []
  for basename in filter(lambda fn: re.match(rf'test_{test_no}_(\d+)', fn), os.listdir('.')):
    test_times = get_test_times(basename)
    tests = []
    for fn in os.listdir(basename):
      f = path.join(basename, fn)
      with open(f, 'r') as file:
        content = file.read()
        matches = re.findall(r"([a-zA-Z]+)\s+time=\s*([\d.]+)s\s+NTR=(\d+)\s+nproc=(\d+)\s+N=(\d+)", content)
        if len(matches) < 4:
          print(f"{f}: Error!\n")
        for m in matches:
          tests.append(Experiment(
            int(m[4]),
            m[0], 
            float(m[1]), 
            int(m[2]), 
            int(m[3]), 
            supercomp,
            None,
            test_times,
          ))
    tests_all = tests_all + tests
  return tests_all

n_p = None
def cond_fix_n_p(e):
  global n_p
  return e.n // e.nproc == n_p and e.op == "solver"

def series_corr_to_n(n):
  def foo(e):
    return e.n == n
  return foo

tests_1 = read_tests(1)

series = rm.divide(tests_1, 'test_no', 'n', 'op')
for ks, s_k in rm.iter_series(series):
  e_start = min(s_k, key=lambda e: e.nproc)
  for e in s_k:
    e.speedup = e_start.time / e.time

with open("exp2.csv", 'w') as f:
  t = rm.table(tests_1, 'nproc', 'op', 'speedup', series_corr_to_n(10**5))
  rm.print_table(t)
  f.write(rm.print_table(t, csv=True))
  f.write("\n")
  t = rm.table(tests_1, 'nproc', 'op', 'speedup', series_corr_to_n(10**6))
  rm.print_table(t)
  f.write(rm.print_table(t, csv=True))
  f.write("\n")
  t = rm.table(tests_1, 'nproc', 'op', 'speedup', series_corr_to_n(10**7))
  rm.print_table(t)
  f.write(rm.print_table(t, csv=True))
  f.write("\n")

test_2 = read_tests(2)
newer = list(filter(lambda e: e.test_no >= 0, test_2))

with open("exp3.csv", 'w') as f:
  n_p = 8192
  t = rm.table(newer, 'nproc', 'n', 'time', cond_fix_n_p)
  rm.print_table(t)
  f.write(rm.print_table(t, csv=True))

  n_p = 8192*2**3
  t = rm.table(newer, 'nproc', 'n', 'time', cond_fix_n_p)
  rm.print_table(t)
  f.write(rm.print_table(t, csv=True))

  n_p = 8192*2**6
  t = rm.table(newer, 'nproc', 'n', 'time', cond_fix_n_p)
  rm.print_table(t)
  f.write(rm.print_table(t, csv=True))

test_3 = read_tests(3)
nproc_1 = list(filter(lambda e: e.nproc == 1, test_3))
series = rm.divide(nproc_1, 'test_no', 'op')
for ks, s_k in rm.iter_series(series):
  e_start = min(s_k, key=lambda e: e.ntr)
  for e in s_k:
    e.speedup = e_start.time / e.time

with open("exp1.csv", 'w') as f:
  t = rm.table(nproc_1, 'ntr', 'op', 'speedup')
  rm.print_table(t)
  f.write(rm.print_table(t, csv=True))