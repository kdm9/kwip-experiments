import numpy as np
import matplotlib.pyplot as plt

arr = np.zeros((96, 96))
for x in range(2):
    for y in range(2):
        for j in range(x * 48, (x+1) * 48):
            for k in range(y * 48, (y+1) * 48):
                if x == y:
                    arr[j, k] = 2
                else:
                    arr[j, k] = 4

for i in range(16):
    st = i * 6
    sp = st + 6
    for j in range(st, sp):
        for k in range(st, sp):
            if j == k:
                arr[j, k] = 0
            else:
                arr[j, k] = 1

plt.imshow(arr, interpolation='none')
plt.show()
