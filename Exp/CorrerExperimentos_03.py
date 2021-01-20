import math, subprocess
import psutil
import numpy as np
from tqdm import tqdm
import csv
import os

def correr_experimento(algoritmo, secuencia1, secuencia2):
    # Crear proceso para ejecutar el codigo.
    process = subprocess.Popen(["../cli", "-a", algoritmo, "-s", secuencia1, "-t", secuencia2], stderr=subprocess.PIPE, stdout=subprocess.PIPE, stdin=subprocess.PIPE, universal_newlines = True)
    psutil.Process(process.pid).nice(-20)

    exit_code = process.wait()

    # Verificar que el proceso no fallo.
    assert exit_code == 0, F"Hubo un error en la experimentacion para el algoritmo: {algoritmo} con la secuencia {secuencia1} y {secuencia2}."
    # Leer salida de STDERR con los tiempos de ejecucion de cada metodo.
    outputs = process.stderr.readline().split();

    tiempo = outputs[0]
    celdas = outputs[1]
    score = outputs[2]
    longitud1 = outputs[3]
    longitud2 = outputs[4]

    process.stdin.close();
    process.stdout.close();
    process.stderr.close();

    return tiempo, celdas, score, longitud1, longitud2

experimentos = [];

# Experimentos SW & NW
randoms = [f"../sequences/random/{name}" for name in os.listdir("../sequences/random/")]

for i in range(len(randoms)):
    random1 = randoms[i]
    random1_len = random1.split('_')[1]
    for j in range(i+1,len(randoms)):
        random2 = randoms[j]
        random2_len = random2.split('_')[1]
        if random1_len != random2_len:
            continue
        for alineamiento in ["SW","NW"]:
            for lenguaje in ["C","ASM"]:
                for tecnologia in ["LIN","SSE","AVX","AVX512"]:
                    if not (lenguaje == "ASM" and tecnologia == "LIN"):
                        algoritmo = f"{alineamiento}_{lenguaje}_{tecnologia}"
                        experimentos.append([f"{alineamiento}_random", algoritmo, random1, random2])

# Correr experimentos
csv_file = open("results/results_03.csv", "w")
writer = csv.writer(csv_file, delimiter=',')

columnas = [
    "experimento",
    "algoritmo",
    "secuencia1",
    "secuencia2",
    "longitud1",
    "longitud2",
    "tiempo",
    "celdas",
    "score"
];

writer.writerow(columnas)

T = 1 # Numero de veces que se ejecuta cada experimento
for experimento in tqdm(experimentos):
    exp = experimento[0]
    algoritmo = experimento[1]
    secuencia1 = experimento[2]
    secuencia2 = experimento[3]

    for i in range(0, T):
        tiempo, celdas, score, longitud1, longitud2 = correr_experimento(algoritmo, secuencia1, secuencia2)

        fila = [
            exp,
            algoritmo,
            secuencia1,
            secuencia2,
            longitud1,
            longitud2,
            tiempo,
            celdas,
            score
        ]

        writer.writerow(fila)
