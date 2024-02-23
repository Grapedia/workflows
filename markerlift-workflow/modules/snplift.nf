process SNPLIFT {
    scratch true

    input:
    val(x)
    val(config)
    val(new_filename)
    val(working_dir)

    output:
    val(config), emit: config_file

    script:
    """
    time /snplift/snplift "$working_dir/$config" > $new_filename
    """
}
