rule start:
    input:
        "summary/disk_space_start.txt" 
    output: 
        touch("flags/00_start.done")

rule check_disk_space:
    output:
        "summary/disk_space_start.txt"
    params:
        # required disk space in Megabytes
        required_space=100_000
    run:
        import shutil
        _,_,free = shutil.disk_usage(".")
        free = free // (1024**2)
        if free < params.required_space:
            raise ValueError(f"Insufficient disk space: {free} Mb available, {params.required_space} Mb required.")
        with open(output[0], 'w') as f:
            f.write(f"Available disk space: {free} Mb\n")
            f.write(f"Required disk space: {params.required_space} Mb\n")
