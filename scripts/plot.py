import matplotlib.pyplot as plt
import csv

p = []
speedup = []

plt.xlabel('Number of threads (p)')
plt.ylabel('Speedup S(p)')

# 1024 matrix, b = 0
with open('../data/matrixMult.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    for row in lines:
        if (row[0] == '1024' and row[1] == '0'):
            p.append(int(row[2]))
            speedup.append(float(row[3]))

plt.plot(p, speedup, color='b', linestyle='dotted', marker='o', label='1024, b=0')

# Ideal speedup
plt.plot(p, p, color='b', linestyle='--', label='Ideal speedup (1024, b=0)')

plt.grid()
plt.legend()
plt.xscale('log')
plt.yscale('linear')
plt.show()


p = []
speedup = []

plt.xlabel('Number of threads (p)')
plt.ylabel('Speedup S(p)')

# 1024 matrix, b = 64
with open('../data/matrixMult.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    for row in lines:
        if (row[0] == '1024' and row[1] == '64'):
            p.append(int(row[2]))
            speedup.append(float(row[3]))

plt.plot(p, speedup, color='r', linestyle='dotted', marker='o', label='1024, b=64')

plt.plot(p, p, color='r', linestyle='--', label='Ideal speedup (1024, b=64)')

plt.grid()
plt.legend()
plt.xscale('log')
plt.yscale('linear')
plt.show()


p = []
speedup = []

plt.xlabel('Number of threads (p)')
plt.ylabel('Speedup S(p)')

# 2048 matrix, b = 0
with open('../data/matrixMult.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    for row in lines:
        if (row[0] == '2048' and row[1] == '0'):
            p.append(int(row[2]))
            speedup.append(float(row[3]))

plt.plot(p, speedup, color='g', linestyle='dotted', marker='o', label='2048, b=0')

# Ideal speedup
plt.plot(p, p, color='g', linestyle='--', label='Ideal speedup (2048, b=0)')

plt.grid()
plt.legend()
plt.xscale('log')
plt.yscale('linear')
plt.show()


p = []
speedup = []

plt.xlabel('Number of threads (p)')
plt.ylabel('Speedup S(p)')

# 2048 matrix, b = 64
with open('../data/matrixMult.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    for row in lines:
        if (row[0] == '2048' and row[1] == '64'):
            p.append(int(row[2]))
            speedup.append(float(row[3]))

plt.plot(p, speedup, color='m', linestyle='dotted', marker='o', label='2048, b=64')

# Ideal speedup
plt.plot(p, p, color='m', linestyle='--', label='Ideal speedup (2048, b=64)')

plt.grid()
plt.legend()
plt.xscale('log')
plt.yscale('linear')
plt.show()
