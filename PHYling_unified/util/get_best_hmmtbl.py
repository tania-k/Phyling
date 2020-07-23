#!/usr/bin/env python
import argparse, sys

parser = argparse.ArgumentParser(description="Get top hit from HMM search",
                                 add_help=True)
parser.add_argument('-c','--cutoff',default='1e-20',type=float,
                    help="HMM parse cutoff")
parser.add_argument('-m','--maxhits',default='1',type=int,
                    help="Maximum best hits to return")
parser.add_argument('-i','--input',help="Input domtbl file")
args = parser.parse_args(sys.argv[1:])

#print(args)
#print(args.cutoff)
#print(args.input)


with open(args.input,"r") as fh:
    seenbest = {}
    gene_matches  = {}
    ortho_matches = {}
    for line in fh:
        if line[0] != "#":
            line = line.strip("\n")
            row = line.split()
            gene = row[0]
            orthohmm = row[3]
            evalue = float(row[6])
            if (float(evalue) <= args.cutoff):
                if gene not in gene_matches:
                    gene_matches[gene] = []
                gene_matches[gene].append([orthohmm,evalue])

                if orthohmm not in ortho_matches:
                    ortho_matches[orthohmm] = {gene: evalue}
                elif ( gene not in ortho_matches[orthohmm] or
                       ortho_matches[orthohmm][gene] > evalue):
                    ortho_matches[orthohmm][gene] = evalue

    for orthohmm in ortho_matches.keys():
        print("orthohmm='%s'"%(orthohmm))
        lastevalue = ''
        for gene in sorted(ortho_matches[orthohmm],
                            key=lambda x: ortho_matches[orthohmm][x], reverse = False):
            evalue = ortho_matches[orthohmm][gene]
            print("  gene is %s evalue is %s"%(gene,evalue))

            if ( len(gene_matches[gene]) == 1 ):
                print("only one match for %s - %s (evalue=%s)" %(gene,orthohmm,evalue))
            else:
                print("more than one match for %s - %s (evalue=%s)" %(gene,orthohmm,evalue))
                for hits in gene_matches[gene]:
                    print("\tgene=%s orthohmm=%s evalue=%s"%(gene,hits[0],hits[1]))


#                (not q in seenbest)):    # and not already seen
#                seenbest[q] = [q, t, evalue]

#    for s in seenbest:
#        print("\t".join(seenbest[s]))
