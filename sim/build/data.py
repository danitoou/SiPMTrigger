file = open('/home/danitoou/geant4/sipmtrigger/build/data.txt')

channels = [[] for i in range(8)]
while True:

    line = file.readline().rstrip('\n').split(',')

    if line[0] == '':
        break


    for i in range(8):
        if line[i] != '0':
            channels[i].append(int(line[i]))


for i in range(8):
    channels[i].sort()

two_ch = [0 for i in range(68)]
three_ch = [0 for i in range(568)]
four_ch = [0 for i in range(4568)]

# 01 02 03 04 05 06 07
# 12 13 14 15 16 17
# 23 24 25 26 27
# 34 35 36 37
# 45 46 47
# 56 57
# 67

# 012 013 014 015 016 017 023 024 025 026 027 034 035 036 037 045 046 047 056 057 067
# 123 124 125 126 127 134 135 136 137 145 146 147 156 157 167
# 234 235 236 237 245 246 247 256 257 267
# 345 346 347 356 357 367
# 456 457 467
# 567


for i in range(8):
    for event in channels[i]:
        for j in range(i+1, 8):
            if event in channels[j]:
                two_ch[10*i + j] += 1
                for k in range(j+1, 8):
                    if event in channels[k]:
                        three_ch[100*i + 10*j + k] += 1
                        for l in range(k+1, 8):
                            if event in channels[l]:
                                four_ch[1000*i + 100*j + 10*k + l] += 1

# two_ch = [i for i in two_ch if i != 0]
# three_ch = [i for i in three_ch if i != 0]
# four_ch = [i for i in four_ch if i != 0]

# print(two_ch)
# print(three_ch)
# print(four_ch)
                                
for i in range(8):
    print("%d - %d" % (i, len(channels[i])))

for i in range(len(two_ch)):
    num = two_ch[i]
    if num != 0:
        print('%dx%d - %d' % (i/10, i%10, num))

for i in range(len(three_ch)):
    num = three_ch[i]
    if num != 0:
        print('%dx%dx%d - %d' % (i/100, (i/10)%10, i%10, num))

for i in range(len(four_ch)):
    num = four_ch[i]
    if num != 0:
        print('%dx%dx%dx%d - %d' % (i/1000, (i/100)%10, (i/10)%10, i%10, num))



file.close()