def main():
    pass

if __name__ == '__main__':
    chr1_file = open('./E116_15_coreMarks_chr1_statebyline.txt')
    # Read in all the lines at once
    contents = chr1_file.read()
    # Dump the first 2 lines since we don't need them
    track = contents.split('\n')[2:]
    bedfile = open('./wgEncodeBroadHistoneGm12878H3k27acStdSig.bedGraph')
    enhancers = []
    curr = []
    for line in bedfile:
        tokens = line.strip().split('\t')
        if track[int(tokens[1]) / 200] == '7':
            curr.append(float(tokens[3]))
        else:
            if curr:
                enhancers.append(curr)
                curr = []

    print enhancers[:10]
