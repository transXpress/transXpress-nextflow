# transXpress-nextflow
transXpress: a Nextflow pipeline for rapid de novo transcriptome assembly and annotation

Also see our sister project: [transXpress-snakemake](https://github.com/transXpress/transXpress-snakemake)

## Intro

## Dependencies

Requires
* NextFlow 19.02.0+
* BioPython
* samtools
* R
* infernal
* seqkit
* basic linux utitilies: wget, split

## Installation


## Usage
<!--
  ~ Copyright 2013-2019, Centre for Genomic Regulation (CRG)
  ~
  ~ Licensed under the Apache License, Version 2.0 (the "License");
  ~ you may not use this file except in compliance with the License.
  ~ You may obtain a copy of the License at
  ~
  ~     http://www.apache.org/licenses/LICENSE-2.0
  ~
  ~ Unless required by applicable law or agreed to in writing, software
  ~ distributed under the License is distributed on an "AS IS" BASIS,
  ~ WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ~ See the License for the specific language governing permissions and
  ~ limitations under the License.
  -->

<html>

    <head>
        <title>Nextflow Cytoscape.js with Dagre</title>

        <meta name="viewport" content="width=device-width, user-scalable=no, initial-scale=1, maximum-scale=1">

        <script type="text/javascript">
        var prot = (("https:" == document.location.protocol) ? "https://" : "http://");
        document.write(unescape("%3Cscript src='" + prot + "code.jquery.com/jquery-2.0.3.min.js' type='text/javascript' %3E%3C/script%3E"));
        document.write(unescape("%3Cscript src='" + prot + "cdnjs.cloudflare.com/ajax/libs/cytoscape/2.6.12/cytoscape.min.js' type='text/javascript' %3E%3C/script%3E"));
        document.write(unescape("%3Cscript src='" + prot + "cdn.rawgit.com/cpettitt/dagre/v0.7.4/dist/dagre.min.js' type='text/javascript' %3E%3C/script%3E"));
        document.write(unescape("%3Cscript src='" + prot + "cdn.rawgit.com/cytoscape/cytoscape.js-dagre/1.1.2/cytoscape-dagre.js' type='text/javascript' %3E%3C/script%3E"));
        </script>

        <style>
            body {
                font-family: helvetica;
                font-size: 14px;
            }

            #cy {
                width: 100%;
                height: 100%;
                position: absolute;
                left: 0;
                top: 0;
                z-index: 999;
            }

            h1 {
                opacity: 0.5;
                font-size: 1em;
            }
        </style>

        <script>
            $(function(){
                var cy = window.cy = cytoscape({
                    container: document.getElementById('cy'),
                    boxSelectionEnabled: false,
                    autounselectify: true,

                    layout: {
                        name: 'dagre'
                    },

                    style: cytoscape.stylesheet()
                        .selector( 'node')
                            .css({
                                'width': 10,
                                'height': 10,
                                'content': 'data(label)',
                                'text-valign': 'center',
                                'text-halign': 'center',
                                'text-opacity': 0.5,
                            })
                        .selector('node.PROCESS')
                            .css({
                                'width': 100,
                                'height': 50,
                                'text-opacity': 0.9,
                                'background-color': '#009911'
                            })
                        .selector('node.OPERATOR')
                            .css({
                                'background-color': '#11479e',
                                'text-halign': 'right',
                            })
                        .selector('node.ORIGIN')
                            .css({
                                'background-color': '#999999',
                                'text-halign': 'right',
                            })
                        .selector('node.TERMINATION')
                            .css({
                                'background-color': '#999999',
                                'text-halign': 'right',
                            })
                        .selector('edge')
                            .css({
                                'content': 'data(label)',
                                'text-opacity': 0.5,
                                'width': 4,
                                'target-arrow-shape': 'triangle',
                                'line-color': '#9dbaea',
                                'target-arrow-color': '#9dbaea'
                            }),

elements: {
nodes: [
{ data: { id: 'p0', label: 'Channel.fromPath'}, classes: 'ORIGIN' },
{ data: { id: 'p1', label: 'splitCsv'}, classes: 'OPERATOR' },
{ data: { id: 'p2', label: 'map'}, classes: 'OPERATOR' },
{ data: { id: 'p3', label: 'trimmomatic'}, classes: 'PROCESS' },
{ data: { id: 'p4'}, classes: 'NODE' },
{ data: { id: 'p5', label: 'into'}, classes: 'OPERATOR' },
{ data: { id: 'p6'}, classes: 'NODE' },
{ data: { id: 'p7', label: 'into'}, classes: 'OPERATOR' },
{ data: { id: 'p8'}, classes: 'NODE' },
{ data: { id: 'p9'}, classes: 'ORIGIN' },
{ data: { id: 'p10', label: 'convertSamplesToRelative'}, classes: 'PROCESS' },
{ data: { id: 'p11', label: 'into'}, classes: 'OPERATOR' },
{ data: { id: 'p12'}, classes: 'NODE' },
{ data: { id: 'p13', label: 'collect'}, classes: 'OPERATOR' },
{ data: { id: 'p14', label: 'collect'}, classes: 'OPERATOR' },
{ data: { id: 'p15', label: 'collect'}, classes: 'OPERATOR' },
{ data: { id: 'p16', label: 'convertReadsToYAML'}, classes: 'PROCESS' },
{ data: { id: 'p17'}, classes: 'NODE' },
{ data: { id: 'p18', label: 'collect'}, classes: 'OPERATOR' },
{ data: { id: 'p19', label: 'collect'}, classes: 'OPERATOR' },
{ data: { id: 'p20', label: 'trinityInchwormChrysalis'}, classes: 'PROCESS' },
{ data: { id: 'p21', label: 'splitText'}, classes: 'OPERATOR' },
{ data: { id: 'p22', label: 'trinityButterflyParallel'}, classes: 'PROCESS' },
{ data: { id: 'p23', label: 'collectFile'}, classes: 'OPERATOR' },
{ data: { id: 'p24', label: 'collect'}, classes: 'OPERATOR' },
{ data: { id: 'p25', label: 'collect'}, classes: 'OPERATOR' },
{ data: { id: 'p26', label: 'trinityFinish'}, classes: 'PROCESS' },
{ data: { id: 'p27'}, classes: 'ORIGIN' },
{ data: { id: 'p28', label: 'renameTrinityAssembly'}, classes: 'PROCESS' },
{ data: { id: 'p29', label: 'transdecoderLongOrfs'}, classes: 'PROCESS' },
{ data: { id: 'p30', label: 'collect'}, classes: 'OPERATOR' },
{ data: { id: 'p31', label: 'collect'}, classes: 'OPERATOR' },
{ data: { id: 'p32', label: 'kallisto'}, classes: 'PROCESS' },
{ data: { id: 'p33', label: 'concat'}, classes: 'OPERATOR' },
{ data: { id: 'p34', label: 'first'}, classes: 'OPERATOR' },
{ data: { id: 'p35', label: 'into'}, classes: 'OPERATOR' },
{ data: { id: 'p36', label: 'trinityStats'}, classes: 'PROCESS' },
{ data: { id: 'p37', label: 'splitFasta'}, classes: 'OPERATOR' },
{ data: { id: 'p38', label: 'into'}, classes: 'OPERATOR' },
{ data: { id: 'p39', label: 'splitFasta'}, classes: 'OPERATOR' },
{ data: { id: 'p40', label: 'into'}, classes: 'OPERATOR' },
{ data: { id: 'p41', label: 'downloadPfam'}, classes: 'PROCESS' },
{ data: { id: 'p42', label: 'downloadVirusesUniref50'}, classes: 'PROCESS' },
{ data: { id: 'p43'}, classes: 'NODE' },
{ data: { id: 'p44', label: 'downloadRfam'}, classes: 'PROCESS' },
{ data: { id: 'p45', label: 'downloadSprot'}, classes: 'PROCESS' },
{ data: { id: 'p46', label: 'sprotBlastxParallel'}, classes: 'PROCESS' },
{ data: { id: 'p47', label: 'sprotBlastpParallel'}, classes: 'PROCESS' },
{ data: { id: 'p48', label: 'collectFile'}, classes: 'OPERATOR' },
{ data: { id: 'p49', label: 'collectFile'}, classes: 'OPERATOR' },
{ data: { id: 'p50', label: 'into'}, classes: 'OPERATOR' },
{ data: { id: 'p51', label: 'pfamParallel'}, classes: 'PROCESS' },
{ data: { id: 'p52', label: 'collectFile'}, classes: 'OPERATOR' },
{ data: { id: 'p53', label: 'collectFile'}, classes: 'OPERATOR' },
{ data: { id: 'p54', label: 'into'}, classes: 'OPERATOR' },
{ data: { id: 'p55', label: 'rfamParallel'}, classes: 'PROCESS' },
{ data: { id: 'p56', label: 'collectFile'}, classes: 'OPERATOR' },
{ data: { id: 'p57', label: 'publishRfamResults'}, classes: 'PROCESS' },
{ data: { id: 'p58'}, classes: 'NODE' },
{ data: { id: 'p59', label: 'transdecoderPredict'}, classes: 'PROCESS' },
{ data: { id: 'p60', label: 'splitFasta'}, classes: 'OPERATOR' },
{ data: { id: 'p61', label: 'into'}, classes: 'OPERATOR' },
{ data: { id: 'p62', label: 'deeplocParallel'}, classes: 'PROCESS' },
{ data: { id: 'p63', label: 'collectFile'}, classes: 'OPERATOR' },
{ data: { id: 'p64', label: 'tmhmmParallel'}, classes: 'PROCESS' },
{ data: { id: 'p65', label: 'collectFile'}, classes: 'OPERATOR' },
{ data: { id: 'p66', label: 'annotatedFasta'}, classes: 'PROCESS' },
{ data: { id: 'p67', label: 'into'}, classes: 'OPERATOR' },
{ data: { id: 'p68'}, classes: 'NODE' },
{ data: { id: 'p69', label: 'into'}, classes: 'OPERATOR' },
{ data: { id: 'p70'}, classes: 'NODE' },
{ data: { id: 'p71', label: 'final_checksum'}, classes: 'PROCESS' },
],
edges: [
{ data: { source: 'p0', target: 'p1'} },
{ data: { source: 'p1', target: 'p2'} },
{ data: { source: 'p2', target: 'p3', label: 'readPairs_ch' } },
{ data: { source: 'p3', target: 'p5', label: 'filteredForwardReads_ch' } },
{ data: { source: 'p3', target: 'p7', label: 'filteredReverseReads_ch' } },
{ data: { source: 'p3', target: 'p15', label: 'filteredSingleReads_ch1' } },
{ data: { source: 'p3', target: 'p4', label: 'filteredSingleReads_ch2' } },
{ data: { source: 'p5', target: 'p13', label: 'filteredForwardReads_ch1' } },
{ data: { source: 'p5', target: 'p18', label: 'filteredForwardReads_ch2' } },
{ data: { source: 'p5', target: 'p24', label: 'filteredForwardReads_ch3' } },
{ data: { source: 'p5', target: 'p6', label: 'filteredForwardReads_ch4' } },
{ data: { source: 'p5', target: 'p30', label: 'filteredForwardReads_ch5' } },
{ data: { source: 'p7', target: 'p14', label: 'filteredReverseReads_ch1' } },
{ data: { source: 'p7', target: 'p19', label: 'filteredReverseReads_ch2' } },
{ data: { source: 'p7', target: 'p25', label: 'filteredReverseReads_ch3' } },
{ data: { source: 'p7', target: 'p8', label: 'filteredReverseReads_ch4' } },
{ data: { source: 'p7', target: 'p31', label: 'filteredReverseReads_ch5' } },
{ data: { source: 'p9', target: 'p10', label: 'samples.txt' } },
{ data: { source: 'p10', target: 'p11', label: 'relative_samples_txt_ch' } },
{ data: { source: 'p11', target: 'p20', label: 'relative_samples_txt_ch1' } },
{ data: { source: 'p11', target: 'p26', label: 'relative_samples_txt_ch2' } },
{ data: { source: 'p11', target: 'p12', label: 'relative_samples_txt_ch3' } },
{ data: { source: 'p11', target: 'p32', label: 'relative_samples_txt_ch4' } },
{ data: { source: 'p13', target: 'p16'} },
{ data: { source: 'p14', target: 'p16'} },
{ data: { source: 'p15', target: 'p16'} },
{ data: { source: 'p16', target: 'p17', label: 'datasets_YAML_ch' } },
{ data: { source: 'p18', target: 'p20'} },
{ data: { source: 'p19', target: 'p20'} },
{ data: { source: 'p20', target: 'p26', label: 'trinityWorkDir' } },
{ data: { source: 'p20', target: 'p21', label: 'trinityCmds' } },
{ data: { source: 'p21', target: 'p22'} },
{ data: { source: 'p22', target: 'p23', label: 'trinityFinishedCmds' } },
{ data: { source: 'p23', target: 'p26'} },
{ data: { source: 'p24', target: 'p26'} },
{ data: { source: 'p25', target: 'p26'} },
{ data: { source: 'p26', target: 'p28', label: 'originalGeneTransMap' } },
{ data: { source: 'p26', target: 'p28', label: 'Trinity_fasta_ch' } },
{ data: { source: 'p27', target: 'p28', label: 'species.txt' } },
{ data: { source: 'p28', target: 'p32', label: 'transcriptomeKallisto' } },
{ data: { source: 'p28', target: 'p29', label: 'transcriptomeTransdecoder' } },
{ data: { source: 'p28', target: 'p59', label: 'transcriptomeTransdecoderPredict' } },
{ data: { source: 'p28', target: 'p36', label: 'transcriptomeStats' } },
{ data: { source: 'p28', target: 'p37', label: 'transcriptomeSplit' } },
{ data: { source: 'p28', target: 'p66', label: 'transcriptomeAnnotation' } },
{ data: { source: 'p28', target: 'p32', label: 'geneTransMap' } },
{ data: { source: 'p29', target: 'p39', label: 'longOrfsProteomeSplit' } },
{ data: { source: 'p29', target: 'p59', label: 'transdecoderWorkDir' } },
{ data: { source: 'p30', target: 'p32'} },
{ data: { source: 'p31', target: 'p32'} },
{ data: { source: 'p32', target: 'p33', label: 'rawKallistoTable' } },
{ data: { source: 'p32', target: 'p33', label: 'normalizedKallistoTable' } },
{ data: { source: 'p33', target: 'p34'} },
{ data: { source: 'p34', target: 'p35'} },
{ data: { source: 'p35', target: 'p66', label: 'transcriptExpression' } },
{ data: { source: 'p35', target: 'p36', label: 'expressionStats' } },
{ data: { source: 'p37', target: 'p38'} },
{ data: { source: 'p38', target: 'p46', label: 'sprotBlastxChunks' } },
{ data: { source: 'p38', target: 'p55', label: 'rfamChunks' } },
{ data: { source: 'p39', target: 'p40'} },
{ data: { source: 'p40', target: 'p47', label: 'sprotBlastpChunks' } },
{ data: { source: 'p40', target: 'p51', label: 'pfamChunks' } },
{ data: { source: 'p41', target: 'p51', label: 'pfamDb' } },
{ data: { source: 'p42', target: 'p43', label: 'virusDb' } },
{ data: { source: 'p44', target: 'p55', label: 'rfamDb' } },
{ data: { source: 'p45', target: 'p46', label: 'sprotDb' } },
{ data: { source: 'p46', target: 'p48', label: 'sprotBlastxResults' } },
{ data: { source: 'p45', target: 'p47', label: 'sprotDb' } },
{ data: { source: 'p47', target: 'p49', label: 'sprotBlastpResults' } },
{ data: { source: 'p48', target: 'p66', label: 'blastxResult' } },
{ data: { source: 'p49', target: 'p50'} },
{ data: { source: 'p50', target: 'p59', label: 'blastpForTransdecoder' } },
{ data: { source: 'p50', target: 'p66', label: 'blastpResult' } },
{ data: { source: 'p51', target: 'p52', label: 'pfamResults' } },
{ data: { source: 'p51', target: 'p53', label: 'pfamDomResults' } },
{ data: { source: 'p52', target: 'p66', label: 'pfamResult' } },
{ data: { source: 'p53', target: 'p54'} },
{ data: { source: 'p54', target: 'p66', label: 'pfamDomResult' } },
{ data: { source: 'p54', target: 'p59', label: 'pfamForTransdecoder' } },
{ data: { source: 'p55', target: 'p56', label: 'rfamResults' } },
{ data: { source: 'p56', target: 'p57'} },
{ data: { source: 'p57', target: 'p58', label: 'rfamResultPub' } },
{ data: { source: 'p59', target: 'p66', label: 'predictProteome' } },
{ data: { source: 'p59', target: 'p60', label: 'predictProteomeSplit' } },
{ data: { source: 'p60', target: 'p61'} },
{ data: { source: 'p61', target: 'p62', label: 'deeplocChunks' } },
{ data: { source: 'p61', target: 'p64', label: 'tmhmmChunks' } },
{ data: { source: 'p62', target: 'p63', label: 'deeplocResults' } },
{ data: { source: 'p63', target: 'p66', label: 'deeplocResult' } },
{ data: { source: 'p64', target: 'p65', label: 'tmhmmResults' } },
{ data: { source: 'p65', target: 'p66', label: 'tmhmmResult' } },
{ data: { source: 'p66', target: 'p67', label: 'transcriptome_annotated_fasta_ch' } },
{ data: { source: 'p66', target: 'p69', label: 'transcriptome_annotated_pep_ch' } },
{ data: { source: 'p67', target: 'p71', label: 'transcriptome_annotated_fasta_ch1' } },
{ data: { source: 'p67', target: 'p68', label: 'transcriptome_annotated_fasta_ch2' } },
{ data: { source: 'p69', target: 'p71', label: 'transcriptome_annotated_pep_ch1' } },
{ data: { source: 'p69', target: 'p70', label: 'transcriptome_annotated_pep_ch2' } },
],
},

                });

            });
        </script>
    </head>

    <body>
        <h1>Nextflow Cytoscape.js with Dagre</h1>
        <div id="cy"></div>
    </body>

</html>
