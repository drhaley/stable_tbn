def main():
    for size in range(2,11):
        for top_amount in range(2,11):
            produce(size, top_amount, include_catalytic_region=False)
        produce(size, "inf")

def produce(size, top_amount, include_catalytic_region=True):
    if type(top_amount) is int:
        top_amount_str = f"{top_amount:03}"
    else:
        top_amount_str = top_amount
    autocat_string = "_autocat" if include_catalytic_region else ""
    filename = f"tbn_gg{autocat_string}_size_{size:03}_amount_{top_amount_str}.txt"
    with open(filename, 'w') as outFile:
        for i in range(size):
            H_domains = ' '.join([f"x{i}{j}" for j in range(size)])
            outFile.write(f"{top_amount}[{H_domains} >H{i}]\n")
        for j in range(size):
            V_list = [f"x{i}{j}" for i in range(size)]
            if include_catalytic_region:
                V_list.extend([f"x{i}{j}" for i in range(j, size)])
            V_domains = ' '.join(V_list)
            outFile.write(f"{top_amount}[{V_domains} >V{j}]\n")

        G_domains = ' '.join([f"x{i}{j}*" for i in range(size) for j in range(size)])
        outFile.write(f"2[{G_domains} >G]\n")

if __name__ == "__main__":
    main()
