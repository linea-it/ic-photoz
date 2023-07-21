from itertools import islice

with open("rail-condor.log") as f:
    for line in islice(f, 100):
        print(line)