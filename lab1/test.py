

import sys
import time
from dateutil.parser import parse


file = open("dannie.txt")




for line in file:
    line = line.strip()
    items = line.split('  ')
    items = list(filter(None, items))
    # tm = time.strptime(items[0], "%a %b %d %H:%M:%S %Y")
    # tm = time.strptime(items[0], "%Y-%M")

    d = parse(items[0])
    hour = d.hour

    if hour > 9 and hour < 18:
        if item[3] == "BUY":
            pass

    print(items)
    # print(tm)



