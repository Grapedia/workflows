process SPLIT_COORDS {

    input:
    path(coords_old)

    output:
        val "START_$coords_old", emit: start_coords
        val "END_$coords_old", emit: end_coords
    """
    #!/usr/bin/python3
    
    f1 = open("$projectDir/START_$coords_old", "w")
    f2 = open("$projectDir/END_$coords_old", "w")

    with open("$projectDir/$coords_old") as f:
        for l in f:
            s = l.strip().split('\t')
            f1.write('\t'.join([s[0], s[1], s[1], s[3]]) + '\\n')
            f2.write('\t'.join([s[0], s[2], s[2], s[3]]) + '\\n')

    f1.close()
    f2.close()
    """
}