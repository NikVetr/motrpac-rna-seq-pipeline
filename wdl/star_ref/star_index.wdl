version 1.0

  task star_index {
      input {
          File reference_fasta
          File annotation_gtf
          String prefix
          Int overhang

          Int memory
          Int disk_space
          Int ncpu
          Int preemptible
          String docker
      }

      command <<<
          mkdir ~{prefix}
          STAR \
          --runMode genomeGenerate \
          --genomeDir ~{prefix} \
          --genomeFastaFiles ~{reference_fasta} \
          --sjdbGTFfile ~{annotation_gtf} \
          --sjdbOverhang ~{overhang} \
          --runThreadN ~{ncpu}
          tar -cvzf ~{prefix}.tar.gz ~{prefix}
      >>>

      output {
          File star_index = "~{prefix}.tar.gz"
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

  workflow star_index_workflow {
      input {
          File reference_fasta
          File annotation_gtf
          String prefix
          Int overhang
          Int memory
          Int disk_space
          Int ncpu
          Int preemptible
          String docker
      }

      call star_index {
          input:
              reference_fasta = reference_fasta,
              annotation_gtf = annotation_gtf,
              prefix = prefix,
              overhang = overhang,
              memory = memory,
              disk_space = disk_space,
              ncpu = ncpu,
              preemptible = preemptible,
              docker = docker
      }

      output {
          File star_index_tar = star_index.star_index
      }
  }