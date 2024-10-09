process CHECK_GENOME_INDEX {

    input:
    path(filename)
    val(working_dir)

    output:
        stdout
    """
    #!/usr/bin/python3
    
    import os

    print(int(os.path.isfile(os.path.join("$working_dir", "$filename" + ".bwt"))))

    """
}