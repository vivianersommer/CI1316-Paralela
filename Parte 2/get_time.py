import os

# SERIAL ----------------------------------------------------------
print("Compilando ...")
os.system('gcc lcs.c -o lcs -O3')

print("Executando ...")

for i in range(1):

    print("Executando " + str(i + 1) + "...")
 
    stream = os.popen('time ./lcs')
    output = stream.read()

    with open('results-serial.txt', 'a') as f:
        f.write(output)

# ----------------------------------------------------------------

# PARALELO -------------------------------------------------------

print("Compilando ...")
os.system('mpicc paralela.c -o paralela -O3')

print("Executando ...")

for i in range(20):

    print("Executando " + str(i + 1) + "...")
 
    stream = os.popen('mpirun -np 8 ./paralela')
    output = stream.read()

    with open('results-paralela.txt', 'a') as f:
        f.write(output)

# ----------------------------------------------------------------
