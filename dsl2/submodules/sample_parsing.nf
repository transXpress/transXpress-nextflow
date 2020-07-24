nextflow.preview.dsl=2

process sampleTSVtoTrimmedRelative {
input:
    path toRelative
output:
    path "relative_samples.tsv"
script:
"""
#!/usr/bin/env python
import re
import os
import os.path

output_handle = open("relative_samples.tsv", "w")

with open("${toRelative}", "r") as input_handle:
      for line in input_handle:
        row = re.split("[\t ]+", line)
        newRow = [row[0],row[1]]
        forwardReads = os.path.basename(row[2]).strip()+".R1-P.qtrim.fastq.gz"
        reverseReads = os.path.basename(row[3]).strip()+".R2-P.qtrim.fastq.gz"
        newRow.append(forwardReads)
        newRow.append(reverseReads)
        newRowString = "\t".join(newRow)+os.linesep
        output_handle.write(newRowString)
      output_handle.close()
"""
}

process tsvToYAML {
    input:
      path rnaSPAdes_sample
    output:
      path "relative_samples.yaml"
script:
    """
    #!/usr/bin/env python
    import re
    import os
    import os.path
    import pprint
    sample_list = []
    with open("${rnaSPAdes_sample}", "r") as input_handle:
      for line in input_handle:
        row = re.split("[\t ]+", line)
        if (len(row) > 3): # paired reads
          paired_dict = {}
          paired_dict['orientation'] = 'fr'
          paired_dict['type'] = 'paired-end'
          f_reads = row[2].strip()
          r_reads = row[3].strip()
          paired_dict['left reads'] = [f_reads]
          paired_dict['right reads'] = [r_reads]
          sample_list.append(paired_dict)
        if (len(row) == 3): # unpaired reads
          unpaired_dict = {}
          unpaired_dict['type'] = 'single'
          u_reads = row[2].strip()
          assert os.path.isfile(u_reads)
          unpaired_dict['single reads'] = [u_reads]
          sample_list.append(unpaired_dict)
    with open("relative_samples.yaml", "w") as output_handle:
      output_handle.write(pprint.pformat(sample_list))
    """
}

workflow parseSamplesTSV {
take: sampleTSV
main:
sampleTSV_ch = Channel.fromPath( sampleTSV , checkIfExists:true)

//Note for map{} below. file() is correct. Channel.fromPath(), and path(), doesn't work.
 sampleTSV_ch.splitCsv(sep:'\t',header:false)
     .map{ row -> 
      return tuple(file(row[2]), file(row[3])) }.set{ reads_ch }

sampleTSVtoTrimmedRelative(sampleTSV_ch)
tsvToYAML(sampleTSVtoTrimmedRelative.out)

emit:
 reads=reads_ch;tsv=sampleTSVtoTrimmedRelative.out;yaml=tsvToYAML.out 
}

workflow {
main:
 println(params.samples)
 parseSamplesTSV(params.samples)
}
