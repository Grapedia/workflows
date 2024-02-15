process SNPLIFT {
    scratch true

    input:
    val(x)
    path(config)
    val(new_filename)

    output:
    path(config), emit: config_file

    script:
    """
    time /snplift/snplift "$projectDir/$config" > $projectDir/output_$new_filename
    """
}
