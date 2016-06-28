#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "pcawg-merge-annotate"
label: "PCAWG pre-merge annotation"
cwlVersion: cwl:draft-3
description: |
    A Docker container for the uniform annotation as part of the meging process.
    ```
    Usage:
    # fetch CWL
    $> dockstore cwl --entry quay.io/ljdursi/pcawg-merge-annotate:1.0.0 > Dockstore.cwl
    # make a runtime JSON template and edit it (or use the content of sample_configs.json in this git repo)
    $> dockstore convert cwl2json --cwl Dockstore.cwl > Dockstore.json
    # run it locally with the Dockstore CLI
    $> dockstore launch --entry quay.io/ljdursi/pcawg-merge-annotate:1.0.0 \
        --json Dockstore.json
    ```

dct:creator:
  "@id": "http://orcid.org/0000-0002-4697-798X"
  foaf:name: Jonathan Dursi
  foaf:mbox: "mailto:jonathan@dursi.ca"

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/ljdursi/pcawg-merge-annotate:1.0.0"

#hints:
#  - class: ResourceRequirement
#    coresMin: 1
#    ramMin: 4092
#    outdirMin: 512000
#    description: "the process requires at least 4G of RAM"

inputs:
  - id: "#variant_type"
    type: string
    default: "SNV"
    description: "Annotate SNV (SNV) or indel (indel) variants"
    inputBinding:
      position: 1

  - id: "#input_vcf"
    type: File
    description: "The VCF file to be annotated"
    format: "http://edamontology.org/format_3016"
    inputBinding:
      position: 2

  - id: "#normal_bam"
    type: File
    description: "The normal BAM file"
    secondaryFiles: ".bai"
    format: "http://edamontology.org/format_2572"
    inputBinding:
      position: 3

  - id: "#tumour_bam"
    type: File
    secondaryFiles: ".bai"
    description: "The tumour BAM file"
    format: "http://edamontology.org/format_2572"
    inputBinding:
      position: 4

  - id: "#output"
    type: string

stdout: $(inputs.output)

outputs:
  - id: "#annotated_vcf"
    type: File
    format: "http://edamontology.org/format_3016"
    description: "The annotated VCF"
    outputBinding:
      glob: $(inputs.output)
    description: "A zip file that contains the HTML report and various graphics."

baseCommand: []
