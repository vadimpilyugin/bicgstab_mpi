import os
import math
import sys
import time

command = "mpisubmit.bg -n %d -m %s -w %s --stdout %s --stderr %s.err main %s 1e-6 50 %d %s"
# command = "echo %d %s %s; touch %s; touch %s.err; echo %s %d %s"

filename = "bluegene_%s_%d_%d_%d-%d"

n_to_t = {
  1: " 00:15:00",
  2: " 00:15:00",
  4: " 00:15:00",
  8: " 00:15:00",
  16: " 00:15:00",
  32: " 00:15:00",
  64: " 00:15:00",
  128: " 00:15:00",
  256: " 00:10:00",
  512: " 00:05:00",
  1024: " 00:03:00",
}

def do(cmd):
  print cmd
  os.system(cmd)

n_times = 35

def test_N_fixed_nproc_changes(skip=[0,0,0]):
  deg = [5,6,7]
  Ns = [10**d for d in deg]
  max_values = [
    "100 100 10",
    "100 100 100",
    "100 200 500",
  ]
  launch_idx = [i for i in range(11)]
  nproc = [2**i for i in launch_idx]
  mode = "smp"
  tests = "2000 500 1000 10"
  name = "test1"
  n_threads = 1
  print 
  for i in range(len(Ns)):
    if skip[i] == 1:
      continue
    for j in launch_idx:
      fn = filename % (name, 10, deg[i], nproc[j], n_times)
      cmd = command % (nproc[j], mode, n_to_t[nproc[j]], fn, fn, max_values[i], n_threads, tests)
      do(cmd)

# ar = [5,4,4,]

# j = 1
# for i in range(11):
#   print ar
#   ar[j % 3] += 1
#   j += 1

def log2(x):
  return int(math.log(x)/math.log(2))

def xyz_from_N(N):
  i = log2(N)
  ar = [i/3, i/3, i/3]
  for k in range(i % 3):
    ar[k] += 1
  for k in range(len(ar)):
    ar[k] = 2 ** ar[k]
  return ar

def test_N_fixed_nproc_and_n_threads_change():
  name = "test2"
  hsh = {
    8192*2**0: [2**i for i in range(0, 11)],
    8192*2**3: [2**i for i in range(0, 11)],
    8192*2**6: [2**i for i in range(0, 11)],
  }
  for k,v in hsh.items():
    print
    for nproc in v:
      N = k * nproc
      max_values = xyz_from_N(N)
      max_values = "%d %d %d" % (max_values[0], max_values[1], max_values[2])
      print "N=%d nx,ny,nz=%s nproc=%d" % (N, max_values, nproc)
      fn = filename % (name, 2, log2(N), nproc, n_times)
      cmd = command % (nproc, "smp", n_to_t[nproc], fn, fn, max_values, 1, "2000 500 1000 10")
      do(cmd)
    sl()

def test_N_and_nproc_fixed_n_threads_changes():
  N = 10**6
  name = "test3"
  max_values = "100 100 100"
  n_tests = "2000 500 1000 10"
  mode = "smp"
  for i in range(11):
    nproc = 2**i
    fn = filename % (name, 10, 6, nproc, n_times)
    cmd = command % (nproc, mode, n_to_t[nproc], fn, fn, max_values, 12345, n_tests)
    do(cmd)

def sl():
  am = 300
  time.sleep(am)

was_nt = []

def mv():
  for test_no in [1,2,3]:
    for nt in was_nt:
      folder = "test_%d_%d" % (test_no, nt)
      cmd = "mkdir -p %s" % folder
      do(cmd)
      cmd = "mv bluegene_test%d_*-%d* %s" % (test_no, nt, folder)
      do(cmd)

while True:
  print("\n\n------------- New cycle! [%d] ---------------" % n_times)
  was_nt.append(n_times)
  # test_N_fixed_nproc_changes(skip=[0,1,1])
  # sl()
  # test_N_fixed_nproc_changes(skip=[1,0,1])
  # sl()
  # test_N_fixed_nproc_changes(skip=[1,1,0])
  # sl()
  # mv()
  test_N_fixed_nproc_and_n_threads_change()
  # sl()
  mv()
  # test_N_and_nproc_fixed_n_threads_changes()
  # sl()
  # mv()
  n_times += 1