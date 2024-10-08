process CHECK_GENOME_INDEX {

    input:
    path(filename)

    output:
        stdout
    """
    #!/usr/bin/python3
    
    import os

    print(int(os.path.isfile(os.path.join("$projectDir", "$filename" + ".bwt"))))

    """
}