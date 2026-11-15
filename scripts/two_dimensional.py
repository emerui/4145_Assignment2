import matplotlib.pyplot as plt
import csv

p = []
speedup = []

plt.xlabel('Number of threads (p)')
plt.ylabel('Speedup S(p)')

# 1024 matrix, 1D partition, b = 64
with open('../data/matrixMult.csv', 'r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    for row in lines:
        if row[0] == '1024' and row[1] == '64':
            p.append(int(row[2]))
            speedup.append(float(row[3]))

plt.plot(p, speedup, marker='o', label='1D partition (b=64)')

# 1024 matrix, 2D square grid, b = 64
p2 = [2, 4, 8, 12, 18, 36]
speedup2 = [1.918, 3.828, 7.155, 10.827, 13.855, 16.391]

plt.plot(p2, speedup2, marker='s', label='2D square grid (b=64)')

# Ideal speedup
plt.plot(p, p, linestyle='--', label='Ideal speedup S(p)=p')

plt.xticks(p)
plt.grid()
plt.legend()
plt.xscale('log')
plt.yscale('linear')
plt.title('1024 x 1024 Matrix')
plt.show()


p = []
speedup = []

plt.xlabel('Number of threads (p)')
plt.ylabel('Speedup S(p)')

# 2048 matrix, 1D partition, b = 64
with open('../data/matrixMult.csv', 'r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    for row in lines:
        if row[0] == '2048' and row[1] == '64':
            p.append(int(row[2]))
            speedup.append(float(row[3]))

plt.plot(p, speedup, marker='o', label='1D partition (b=64)')

# 2048 matrix, 2D square grid, b = 64
p2 = [2, 4, 8, 12, 18, 36]
speedup2 = [2.968, 5.896, 11.718, 17.523, 24.104, 39.259]

plt.plot(p2, speedup2, marker='s', label='2D square grid (b=64)')

# 2048 matrix, 2D column grid, 1 x 36
plt.scatter([36], [19.966], marker='D',
            label='2D column grid (1x36)')

# Ideal speedup
plt.plot(p, p, linestyle='--', label='Ideal speedup S(p)=p')

plt.xticks(p)
plt.grid()
plt.legend()
plt.xscale('log')
plt.yscale('linear')
plt.title('2048 x 2048 Matrix')
plt.show()
