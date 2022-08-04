import os

print("Compilando ...")
os.system('gcc -fopenmp lcs.c -O3')

print("Executando ...")

for i in range(20):

    print("Executando " + str(i) + "...")
 
    stream = os.popen('time ./a.out')
    output = stream.read()

    with open('results-serial.txt', 'a') as f:
        f.write(output)

print("Compilando ...")
os.system('gcc -fopenmp paralela.c -O3 -lm')

print("Executando ...")

for i in range(20):

    print("Executando " + str(i) + "...")
 
    stream = os.popen('time ./a.out')
    output = stream.read()

    with open('results-paralela.txt', 'a') as f:
        f.write(output)