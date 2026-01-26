version 1.0

task rsem_reference {
    input {
        File reference_fasta
        File annotation_gtf
        String prefix
        Int memory
        Int disk_space
        Int ncpu
        Int preemptible
        String docker
    }

    command <<<
        mkdir ~{prefix}
        cd ~{prefix} || exit 126
        rsem-prepare-reference --gtf ~{annotation_gtf} --num-threads ~{ncpu} ~{reference_fasta} rsem_reference
        cd .. && tar -cvzf ~{prefix}.tar.gz ~{prefix}
    >>>

    output {
        File rsem_reference = "~{prefix}.tar.gz"
    }

    runtime {
        cpu: ncpu
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        docker: docker
        preemptible: preemptible
    }

    meta {
        author: "Archana Raja"
    }
}

workflow rsem_reference_workflow {
    input {
        File reference_fasta
        File annotation_gtf
        String prefix
        Int memory
        Int disk_space
        Int ncpu
        Int preemptible
        String docker
    }

    call rsem_reference {
        input:
            reference_fasta = reference_fasta,
            annotation_gtf = annotation_gtf,
            prefix = prefix,
            memory = memory,
            disk_space = disk_space,
            ncpu = ncpu,
            preemptible = preemptible,
            docker = docker
    }

    output {
        File rsem_reference_tar = rsem_reference.rsem_reference
    }
}
