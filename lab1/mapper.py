
#!/usr/bin/python

import sys
import time
from dateutil.parser import parse

for line in sys.stdin:
    line = line.strip()
    items = line.split('\t')


    d = parse(items[0])
    hour = d.hour

    if hour > 9 and hour < 18:
        if items[3] == "BUY":
            sys.stdout.write("{}\t{}\n".format(items[1]))



