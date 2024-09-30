process MERGE_RESULT {

    input:
    val(x)
    path(start_coords)
    path(end_coords)
    val(coords_new)

    output:
        val coords_new

    """
    #!/usr/bin/python3

    import os
    
    d_start = {}
    d_end = {}
    with open(os.path.join(os.path.dirname("$coords_new"), "$start_coords")) as f:
      for l in f:
        s = l.strip().split('\t')
        d_start[s[3]] = s

    with open(os.path.join(os.path.dirname("$coords_new"), "$end_coords")) as f:
      for l in f:
        s = l.strip().split('\t')
        d_end[s[3]] = s

    with open("$coords_new", "w") as f:
      for k, s in d_start.items():
        if k not in d_end:
          continue
        e = d_end[k]
        if s[0] != e[0]:
          continue
        f.write('\t'.join([s[0], s[1], e[1], k]) + '\\n')
    
    """
}
