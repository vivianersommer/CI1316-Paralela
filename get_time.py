import os

print("Compilando ...")
os.system('gcc -fopenmp paralela.c')

print("Executando ...")

for i in range(5):

    print("Executando " + str(i) + "...")
 
    stream = os.popen('time ./a.out')
    output = stream.read()

    with open('results.txt', 'a') as f:
        f.write('Teste ' + str(i) + ' -------------------')
        f.write(output)
        f.write('-------------------------  ')
        f.write('\n')
