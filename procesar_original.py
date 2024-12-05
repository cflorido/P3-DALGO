import re

def procesar_archivo(entrada, salida):
    with open(entrada, 'r') as f:
        lineas = f.readlines()

    resultado = []
    i = 0

    while i < len(lineas):
        linea = lineas[i].strip()
        
        # Si la línea contiene dos números al inicio, conservarla
        if re.match(r'^\d+ \d+', linea):
            resultado.append(linea)
            numeros = list(map(int, linea.split()))
            num_lineas = numeros[0]  
            

            for j in range(i + 1, i + 1 + num_lineas):
                if j < len(lineas):  
                    elementos = lineas[j].strip().split()
                    
                    # Eliminar el cuarto número si existe
                    if len(elementos) > 3:
                        elementos.pop(3)
                    
                    resultado.append(" ".join(elementos))
            

            i += num_lineas
        else:
         
            resultado.append(linea)
        
        i += 1

  
    with open(salida, 'w') as f:
        f.write("\n".join(resultado))


procesar_archivo("P2_cases.in", "P2_cases_Def.in")
