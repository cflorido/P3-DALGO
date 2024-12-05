import math
from collections import defaultdict

def construir_grafo(celulas, d):
    """
    Construye un grafo donde las células están conectadas si están a una distancia <= d 
    y comparten al menos un péptido.
    """
    n = len(celulas)
    grafo = defaultdict(list)
    
    for i in range(n):
        id1, x1, y1, peptidos1 = celulas[i]
        for j in range(i + 1, n):
            id2, x2, y2, peptidos2 = celulas[j]
            
            # Calcular distancia entre células
            distancia = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
            
            # Si están dentro de la distancia y comparten péptidos, conectar en el grafo
            if distancia <= d and peptidos1 & peptidos2:
                grafo[id1].append(id2)
                grafo[id2].append(id1)
    
    return grafo

def greedy_componentes_clique(grafo):
    """
    Asigna greedy las células a cliques minimizando el número de cliques necesarios.
    """
    cliques = {}
    
    for celula in sorted(grafo.keys()):
        # Encuentra los cliques vecinos ya asignados
        cliques_vecinos = {cliques[vecino] for vecino in grafo[celula] if vecino in cliques}
        
        # Asigna el menor número de clique disponible
        clique_asignado = 1
        while clique_asignado in cliques_vecinos:
            clique_asignado += 1
        cliques[celula] = clique_asignado
    
    return cliques

def resolver_caso(caso):
    """
    Procesa un caso específico.
    """
    n, d = caso["n"], caso["d"]
    celulas = []
    for entrada in caso["celulas"]:
        id_celula = entrada[0]
        x, y = entrada[1], entrada[2]
        peptidos = set(entrada[3:])
        celulas.append((id_celula, x, y, peptidos))

    # Construir grafo basado en las restricciones
    grafo = construir_grafo(celulas, d)

    # Resolver el problema con un enfoque greedy
    resultado = greedy_componentes_clique(grafo)

    return resultado

def main():
    casos = [
        {
            "n": 7,
            "d": 1,
            "celulas": [
                [1, 0, 0, "AETQT", "DFTYA", "PHLYT"],
                [2, 0, 2, "DSQTS", "IYHLK", "LHGPS", "LTLLS"],
                [3, 1, 0, "AETQT", "DFTYA", "HGCYS", "LSVGG", "SRFNH"],
                [4, 1, 1, "DFTYA", "HGCYS", "IYHLK", "SRFNH"],
                [5, 1, 2, "DSQTS", "IYHLK", "LSVGG", "LTLLS", "TTVTG"],
                [6, 2, 1, "AETQT", "HGCYS", "IYHLK", "LSVGG", "LTLLS"],
                [7, 2, 2, "HGCYS", "SRFNH", "TTVTG"],
            ]
        },
        {
            "n": 4,
            "d": 1,
            "celulas": [
                [1, 0, 0, "AETQT", "DFTYA"],
                [2, 0, 1, "AETQT", "HGCYS"],
                [3, 1, 0, "DFTYA", "IYHLK"],
                [4, 1, 1, "HGCYS", "IYHLK"],
            ]
        }
    ]

    for i, caso in enumerate(casos, start=1):
        resultado = resolver_caso(caso)
        print(f"Caso {i}:")
        for id_celula in sorted(resultado):
            print(id_celula, resultado[id_celula])

if __name__ == "__main__":
    main()
