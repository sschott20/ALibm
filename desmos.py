# read in file poly.txt and print the polynomial where each line a coefficient
# output of the form a0 + a1x + a2x^2 + ... + anx^n

infile = open("../dump/FPoly.txt", "r")
lines = infile.readlines()
combined = ""
for line in lines:
    combined += line.strip() 
    combined += "x^" + str(lines.index(line)) + " + "

print(combined.strip())
