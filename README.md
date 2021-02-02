# RTNano
Real-Time Nanopore sequencing monitor (RTNano)

RTNano is a tool for real-time analysis of Nanopore sequencing data to detect pathogen-positive samples and call variants. The analysis is ultra-fast, without consuming too much PC resourse. It is used in NIRVANA (Nanopore sequencing of Isothermal Rapid Viral Amplification for Near real-time Analysis) and useful in the SARS-CoV-2 related Nanopore sequencing.

More detail: [coming...]

## INSTALLATION
### Prerequisites
Anaconda (https://www.anaconda.com/distribution/) or Miniconda (https://conda.io/miniconda.html).
For example, users may download and install the following Anaconda3 package:
```
wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh
```

### Download the RTNano package
```
git clone https://github.com/milesjor/RTNano.git
cd ./RTNano/
```

### Install all required modules
```
conda env create --name rtnano --file ./conda_env.yaml
conda activate rtnano
python -m pip install pandas
```

## USAGE

### The parameters of RTNano
Users can use the following command to print out usage:
```
python ./rt_nano.py -h
```

A list of full usage:
```
usage: rt_nano.py [-h] -p PATH [-s SAVE_PATH] [-r REFER_SEQ] [-t THREAD]
                  [-T INTERVAL_TIME] [-g GUPPY_BARCODER] [-k BARCODE_KITS]
                  [--run_time RUN_TIME] [--align_identity ALIGN_IDENTITY]
                  [--covered_pcent COVERED_PCENT]
                  [--housekeep_gene HOUSEKEEP_GENE] [--ntc NTC] [--resume]
                  [--put_back] [--call_variant] [--alle_freq ALLE_FREQ] [-v]

Real-Time analysis of Nanopore data for Covid-19 sequencing.

optional arguments:
  -h, --help            show this help message and exit
  -p PATH, --path PATH  path/to/nanopore_result_folder
  -s SAVE_PATH, --save_path SAVE_PATH
                        path/to/saved_folder Default: rtnano_result in -p PATH
                        folder
  -r REFER_SEQ, --refer_seq REFER_SEQ
                        path/to/reference_genome.fa, default is using amplicon.fa
                        in program folder
  -t THREAD, --thread THREAD
                        working thread [1]
  -T INTERVAL_TIME, --interval_time INTERVAL_TIME
                        interval time for scanning minknow folder in minutes
                        [1]
  -g GUPPY_BARCODER, --guppy_barcoder GUPPY_BARCODER
                        Optional: path/to/guppy_barcoder, when offering this
                        parameter, it will do additional demultiplexing using
                        guppy_barcoder --require_barcodes_both_ends
                        --trim_barcodes
  -k BARCODE_KITS, --barcode_kits BARCODE_KITS
                        barcode kits used, e.g. "EXP-NBD114 EXP-NBD104" it is
                        required when providing -g/--guppy_barcoder
  --run_time RUN_TIME   total run time in hours [48]
  --align_identity ALIGN_IDENTITY
                        minimum [0.89] alignment identity when considering an
                        effective alignment
  --covered_pcent COVERED_PCENT
                        minimum covering [0.96] of amplicon length when
                        considering an effective alignment
  --housekeep_gene HOUSEKEEP_GENE
                        the amplicon name of housekeeping gene in the
                        sequencing, default: [ACTB_263bp]
  --ntc NTC             barcode number of no template control in the
                        sequencing, e.g. barcode96
  --resume              resume the unexpectedly interrupted analysis. Please
                        use the same [-p] [-s] as before
  --put_back            return fastq file to their original fastq_pass folder.
                        Please use the same [-p] [-s] [-g] as you generated
                        the result, together with --put_back
  --call_variant        call variants using samtools and filter by alleic
                        frequency (>=0.5). Please use the same [-p] [-s] as
                        you generated the result. It uses fastq file in
                        analyzed_achieve/accumulated_reads folder. If you want
                        to use it for your own data, please put fastq file in
                        this folder, one sample one fastq file
  --alle_freq ALLE_FREQ
                        filter SNVs by allelic frequency [0.5]
  -v, --version         show the current version
```

## EXAMPLE

### Real-time analysis during sequencing run
```
python ./rt_nano.py -p ./example/
```

### Results
```
./example/rtnano_result/
├── 20200607_22.20.29_result.txt
├── 20200607_22.20.29_rt_nano.log
├── analyzed_achieve
│   ├── 20200607_22.20.29
│   ├── accumulated_reads
│   └── pooled_result.txt
└── analyzing
```

### One can also open a new shell and check variants using all of current reads by running:
```
conda activate rtnano
python ./rt_nano.py -p ./example/ --call_variant
```

### Results
```
./example/rtnano_result/snv/
└── 20200607_22.26.20
    ├── barcode15
    ├── barcode16
    ├── barcode17
    ├── barcode18
    ├── barcode19
    ├── barcode20
    └── unclassified
```

## TIPS

### Provide guppy_barcoder
It is strongly recommended to download guppy from Nanopore community, and provide guppy_barcoder for RTNano to ensure a confident demultiplexing.
Run the command like below:
```
python ./rt_nano.py -p ./example/ -g ~/Downloads/ont-guppy-cpu/bin/guppy_barcoder -k "EXP-NBD114"
```
