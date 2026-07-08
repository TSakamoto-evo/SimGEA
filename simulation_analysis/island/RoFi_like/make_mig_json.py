output = open("mig.json", mode="w")
print("{", file=output)

with open("snpsfile.txt", mode="r") as f:
  lineno = 0
  
  for line in f:
    line = line.rstrip("\n")
    parts = line.split()

    n_pop = int(len(parts) / 2)

    if lineno == 0:
      print("\t\"N_patch\": " + str(n_pop) + ",", file=output)

      dij = [0.0 for i in range(int(n_pop * (n_pop - 1) / 2))]
    
    idx = 0
    for i in range(n_pop - 1):
      for j in range(i + 1, n_pop):
        a1 = int(parts[2 * i])
        a0 = int(parts[2 * i + 1])

        b1 = int(parts[2 * j])
        b0 = int(parts[2 * j + 1])

        c0 = (a0 * (a0 - 1) / 2.0) / ((a0 + a1) * (a0 + a1 - 1) / 2.0)
        c1 = (a0 * a1) / ((a0 + a1) * (a0 + a1 - 1) / 2.0)
        c2 = (a1 * (a1 - 1) / 2.0) / ((a0 + a1) * (a0 + a1 - 1) / 2.0)

        d0 = (b0 * (b0 - 1) / 2.0) / ((b0 + b1) * (b0 + b1 - 1) / 2.0)
        d1 = (b0 * b1) / ((b0 + b1) * (b0 + b1 - 1) / 2.0)
        d2 = (b1 * (b1 - 1) / 2.0) / ((b0 + b1) * (b0 + b1 - 1) / 2.0)

        dij[idx] += 1 * c1 * d0 + 4 * c2 * d0 + 1 * c0 * d1 + 1 * c2 * d1 + 4 * c0 * d2 + 1 * c1 * d2

        idx += 1

    lineno += 1

print("\t\"Dij\": [" + str(dij[0] / lineno), end="", file=output)
for i in range(1, len(dij)):
  print(", " + str(dij[i] / lineno), end="", file=output)
print("] ", file=output)
print("}", file=output)

output.close()
