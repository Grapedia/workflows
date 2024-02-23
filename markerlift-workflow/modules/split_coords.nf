process SPLIT_COORDS {

    input:
    val(coords_old)

    output:
        val "$coords_old" + "_START", emit: start_coords
        val "$coords_old" + "_END", emit: end_coords
    """
    #!/usr/bin/python3
    
    f1 = open("$coords_old" + "_START", "w")
    f2 = open("$coords_old" + "_END", "w")

    with open("$coords_old") as f:
        for l in f:
            s = l.strip().split('\t')
            f1.write('\t'.join([s[0], s[1], s[1], s[3]]) + '\\n')
            f2.write('\t'.join([s[0], s[2], s[2], s[3]]) + '\\n')

    f1.close()
    f2.close()
    """
}