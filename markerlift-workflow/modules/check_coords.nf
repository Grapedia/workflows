process CHECK_COORDS {

    input:
    path(coords_old)

    output:
        stdout
    """
    #!/usr/bin/python3

    single = True
    with open("$projectDir/$coords_old") as f:
      for l in f:
        s = l.strip().split('\t')
        if len(s) < 4 or int(s[1]) != int(s[2]):
            single = False
            break
    print(int(single))
    """
}